# ----------------------------------------------------------------------------
# Authors: Bill Tatsis, Matt Seddon, Dan Mason, and Dan O'Donovan
# Company: Healx
# Date: August 17, 2023
#
# Description: This script generates 3D conformers for a given Mol.
# ----------------------------------------------------------------------------

from typing import Tuple

from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem, Mol


def generate_conformers(
    molecule: Mol,
    num_conformers: int,
    prune_rms_threshold: float = 0.5,
    add_hydrogens: bool = True,
    optimize: bool = False,
) -> Mol:
    """
    Generate RDKit conformers for a Mol object.

    Parameters
    ----------
    molecule: RDKit Mol object
        The molecule for which conformers need to be generated.
    num_conformers: int
        Number of conformers to generate.
    prune_rms_threshold: float, optional
        RMS threshold for pruning, default 0.5.
    add_hydrogens: bool, optional
        Whether to add Hydrogen atoms to the Mol, default True.
    optimize: bool, optional
        Whether to optimize conformer generation, default False.

    Returns
    -------
    Mol
        Mutated Mol object with conformers generated.
    """

    molecule.RemoveAllConformers()

    if add_hydrogens:
        molecule = Chem.AddHs(molecule)

    AllChem.EmbedMultipleConfs(
        molecule,
        numConfs=num_conformers,
        maxAttempts=100,
        pruneRmsThresh=prune_rms_threshold,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True,
        useRandomCoords=True,
        randomSeed=42,
        enforceChirality=True,
        numThreads=0,
    )

    if optimize:
        AllChem.MMFFOptimizeMoleculeConfs(molecule, mmffVariant="MMFF94s")

    # Trying to catch some edge cases (mainly with bridged ring systems)
    if molecule.GetNumConformers() > 0:
        return molecule
    else:
        logger.error(
            "Failed to generate 3D confs for "
            + Chem.MolToSmiles(Chem.RemoveHs(molecule)),
        )


def conf_gen(molecule: Mol, mol_id: str, num_conformers: int) -> Tuple[str, Mol]:
    """
    Generate conformers for an RDKit Mol.

    Parameters
    ----------
    molecule: RDKit Mol object
        The molecule for which conformers need to be generated.
    mol_id: str
        A string representing the molecule ID.
    num_conformers: int
        Number of conformers to be generated.

    Returns
    -------
    Tuple[str, Mol]
        A tuple containing the molecule ID and the Mol object mutated to have conformers
        computed.
    """
    if isinstance(molecule, Mol):
        h_mol = generate_conformers(
            molecule, num_conformers, prune_rms_threshold=0.5, add_hydrogens=True
        )

        return mol_id, h_mol


# Helper function
# We can also use functools.partial to create a partial function of conf_gen
# map_conf_gen = partial(conf_gen)


def map_conf_gen(args):
    return conf_gen(*args)
