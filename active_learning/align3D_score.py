# ----------------------------------------------------------------------------
# Authors: Bill Tatsis, Matt Seddon, Dan Mason, and Dan O'Donovan
# Company: Healx
# Date: August 17, 2023
#
# Description: This script performs a 3D overlay of a probe and reference molecule
#              and returns similarity scores of the two compounds
# ----------------------------------------------------------------------------

import collections
import copy
import os
from typing import Optional

import espsim
from rdkit import Chem, RDConfig, RDLogger
from rdkit.Chem import AllChem, Mol, rdMolAlign, rdMolDescriptors, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps

# Silence non-critical RDKit warnings to minimize unnecessary outputs
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# Set up features to use in FeatureMap
fdef_path = os.path.join(RDConfig.RDDataDir, "BaseFeatures.fdef")
FEATURE_DEFINITIONS = AllChem.BuildFeatureFactory(fdef_path)

FEATURE_MAP_PARAMS = {}
for k in FEATURE_DEFINITIONS.GetFeatureFamilies():
    feature_params = FeatMaps.FeatMapParams()
    FEATURE_MAP_PARAMS[k] = feature_params

MOL_FEATURES = (
    "Donor",
    "Acceptor",
    "NegIonizable",
    "PosIonizable",
    "ZnBinder",
    "Aromatic",
    "Hydrophobe",
    "LumpedHydrophobe",
)


AlignmentScore = collections.namedtuple(
    "AlignmentScore", ["shape_score", "esp_score", "rdkit_score", "best_mol"]
)


def _align(
    prb_mol: Mol,
    ref_mol: Mol,
    prb_crippen: Optional[list] = None,
    ref_crippen: Optional[list] = None,
    i: int = -1,
    j: int = -1,
):
    """
    Alignment helper function.

    Parameters
    ----------
    prb_mol: Prb Mol object
    ref_mol: Reference Mol object
    prb_crippen: Optional list of Crippen contributions, default is None
    ref_crippen: Optional list of Crippen contributions, default is None
    i: i index
    j: j index

    Returns
    -------
    None, note that this modifies the alignment object!
    """
    if prb_crippen is None:
        prb_crippen = rdMolDescriptors._CalcCrippenContribs(prb_mol)
    if ref_crippen is None:
        ref_crippen = rdMolDescriptors._CalcCrippenContribs(ref_mol)
    alignment = rdMolAlign.GetCrippenO3A(
        prb_mol, ref_mol, prb_crippen, ref_crippen, i, j
    )
    alignment.Align()


def score_alignment(
    prb_mol: Mol, ref_mol: Mol, partial_charges: str = "gasteiger"
) -> AlignmentScore:
    """
    Score the 3D overlay of two molecules (probe against reference).

    Parameters
    ----------
    prb_mol: RDKit Mol object to score
    ref_mol: RDKit Mol object reference
    partial_charges: str defining partial charges

    Returns
    -------
    AlignmentScore named tuple with scores and Mol object with best score
    """
    # Copy mols, otherwise the alignment will change the positions
    # of the molecule object permanently
    prb_mol = copy.deepcopy(prb_mol)
    ref_mol = copy.deepcopy(ref_mol)

    prb_crippen = rdMolDescriptors._CalcCrippenContribs(prb_mol)
    ref_crippen = rdMolDescriptors._CalcCrippenContribs(ref_mol)

    # Get best conformer
    best_score = 0
    best_i = 0

    for i in range(prb_mol.GetNumConformers()):
        _align(prb_mol, ref_mol, prb_crippen, ref_crippen, i, 0)
        shape_score = espsim.GetShapeSim(prb_mol, ref_mol, i, 0)
        if shape_score > best_score:
            best_score = shape_score
            best_i = i

    # Use best conformer for scoring
    _align(prb_mol, ref_mol, prb_crippen, ref_crippen, best_i, 0)

    shape_score = espsim.GetShapeSim(prb_mol, ref_mol, best_i, 0)
    esp_score = espsim.GetEspSim(
        prb_mol, ref_mol, best_i, 0, partialCharges=partial_charges, metric="tanimoto"
    )

    sc_rdkit_score = rdkit_score(prb_mol, ref_mol, best_i, 0)

    return AlignmentScore(
        shape_score, esp_score, sc_rdkit_score, Chem.Mol(prb_mol, confId=best_i)
    )


def rdkit_score(query_mol: Mol, ref_mol: Mol, query_id: int, ref_id: int) -> float:
    """
    Function to calculate the feature map score

    Parameters
    ----------
    query_mol: RDKit mol object of the query molecule
    ref_mol: RDKit mol object of the reference molecule
    query_id: id of query mol
    ref_id: id of reference mol

    Returns
    -------
    shape and color similarity score
    """

    query_mol = Chem.Mol(query_mol, confId=int(query_id))
    ref_mol = Chem.Mol(ref_mol, confId=int(ref_id))

    feature_lists = []

    for m in [query_mol, ref_mol]:
        rawFeats = FEATURE_DEFINITIONS.GetFeaturesForMol(m)
        # filter that list down to only include the ones we're interested in
        feature_lists.append([f for f in rawFeats if f.GetFamily() in MOL_FEATURES])
    fms = [
        FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=FEATURE_MAP_PARAMS)
        for x in feature_lists
    ]

    fms[0].scoreMode = FeatMaps.FeatMapScoreMode.Best

    fm_score = fms[0].ScoreFeats(feature_lists[1]) / min(
        fms[0].GetNumFeatures(), len(feature_lists[1])
    )

    protrude_dist = rdShapeHelpers.ShapeProtrudeDist(
        query_mol, ref_mol, allowReordering=False
    )
    rdkit_score = 0.5 * fm_score + 0.5 * (1 - protrude_dist)

    return rdkit_score
