import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from chemographykit.gtm import GTM
from chemographykit.utils.molecules import calculate_latent_coords
import torch


def clean_and_scale_descriptors(desc):
    """
    Cleans and scales an array-like of molecular descriptor vectors.

    Parameters
    ----------
    desc : array-like of array-like or np.ndarray
        Nested list or 2D numpy array of descriptors, may contain NaNs.

    Returns
    -------
    scaled_desc : np.ndarray
        2D numpy array, columns with NaNs removed, values standardized (zero mean, unit variance).
    kept_indices : list
        Indices of columns that were kept after NaN filtering.
    """
    # Convert to 2D numpy array
    desc = np.asarray(desc)
    if desc.ndim != 2:
        raise ValueError("Input descriptors must be a 2D array or nested list.")

    # Identify columns without any NaN/infinite values
    valid_mask = ~np.any(np.isnan(desc) | np.isinf(desc), axis=0)
    if not np.any(valid_mask):
        raise ValueError("All descriptor columns contain NaN or infinite values.")
    clean_desc = desc[:, valid_mask]

    # Standardize features (zero mean, unit variance)
    scaler = StandardScaler()
    scaled_desc = scaler.fit_transform(clean_desc)

    kept_indices = np.where(valid_mask)[0].tolist()
    return scaled_desc, kept_indices


def calc_gtm(
    desc,
    device=None,
    num_nodes=36**2,
    num_basis_functions=31**2,
    basis_width=5.726809,
    reg_coeff=543.61223,
    max_iter=300,
    tolerance=0.001,
    standardize=False,
    seed=1234,
    pca_scale=True,
):
    """
    Fits a GTM (Generative Topographic Mapping) model to the given descriptors.

    Parameters
    ----------
    desc : np.ndarray
        Descriptor matrix (n_samples x n_features).
    device : str or torch.device, optional
        Device to use for computation.
    num_nodes : int, optional
        Number of grid nodes in GTM map.
    num_basis_functions : int, optional
        Number of basis functions for GTM.
    basis_width : float, optional
        Basis function width.
    reg_coeff : float, optional
        Regularization coefficient for GTM.
    max_iter : int, optional
        Maximum number of EM iterations.
    tolerance : float, optional
        EM convergence tolerance.
    standardize : bool, optional
        Whether to internally standardize descriptors.
    seed : int, optional
        Random seed for reproducibility.
    pca_scale : bool, optional
        Whether to scale using PCA during GTM initialization.

    Returns
    -------
    gtm : GTM object
        Trained GTM object.
    crds_2d : np.ndarray
        2D GTM projections for each sample.
    resps : np.ndarray
        GTM responsibilities for each sample (n_samples x n_nodes).
    llhs : np.ndarray
        Log-likelihoods for each sample.
    """

    # Defensive device resolution
    if device is None:
        device = torch.device("cpu")

    gtm = GTM(
        num_nodes=num_nodes,
        num_basis_functions=num_basis_functions,
        basis_width=basis_width,
        reg_coeff=reg_coeff,
        max_iter=max_iter,
        tolerance=tolerance,
        standardize=standardize,
        seed=seed,
        device=device,
        pca_scale=pca_scale,
    )
    if not isinstance(desc, torch.Tensor):
        desc = torch.tensor(desc, dtype=torch.float64, device=device)
    gtm.fit_transform(desc)
    resps, llhs = gtm.project(desc)
    resps = resps.detach().cpu().numpy()
    llhs = llhs.detach().cpu().numpy()
    crds_2d = calculate_latent_coords(resps, correction=True, return_node=True)

    return gtm, crds_2d, resps, llhs


def plot_gtm_class_map(resp, gtm_coords, labels, class_name="Active"):
    """
    Plots a GTM class map showing the spatial distribution of a class probability
    and overlays molecule points, highlighting misclassifications by dynamically sizing markers.

    Parameters
    ----------
    resp : np.ndarray
        The Responsibility Matrix (n_molecules x n_nodes).
    gtm_coords : np.ndarray
        GTM node coordinates (n_molecules x 2).
    labels : array-like
        Class labels (same length as number of molecules).
    class_name : str, optional
        Name of the class to show in the plot title/colorbar.

    Returns
    -------
    prob_map : np.ndarray
        1D array of class probabilities for each node.
    """

    n_molecules, n_nodes = resp.shape
    grid_dim = int(np.sqrt(n_nodes))

    # Robust guard for non-square grids
    if grid_dim * grid_dim != n_nodes:
        grid_dim = int(np.round(np.sqrt(n_nodes)))
        if grid_dim * grid_dim != n_nodes:
            print(f"Warning: nodes ({n_nodes}) not a perfect square. Adjusting plot.")

    # Normalize/scale coordinates to imshow grid
    x_span = gtm_coords[:, 0].max() - gtm_coords[:, 0].min()
    y_span = gtm_coords[:, 1].max() - gtm_coords[:, 1].min()
    x_scaled = (
        (gtm_coords[:, 0] - gtm_coords[:, 0].min())
        / (x_span if x_span != 0 else 1)
        * (grid_dim - 1)
    )
    y_scaled = (
        (gtm_coords[:, 1] - gtm_coords[:, 1].min())
        / (y_span if y_span != 0 else 1)
        * (grid_dim - 1)
    )

    # Prepare label array
    y = np.asarray(labels).flatten().astype(float)

    # Weighted sum and normalization
    class_sum_per_node = np.dot(y, resp)
    total_weight_per_node = np.sum(resp, axis=0)

    # Avoid divide-by-zero warnings and invalid values
    small_val = 1e-10
    prob_map = np.divide(
        class_sum_per_node,
        total_weight_per_node,
        out=np.full_like(class_sum_per_node, np.nan),
        where=total_weight_per_node > small_val,
    )

    # Reshape for 2D plotting
    if grid_dim * grid_dim == n_nodes:
        prob_map_2d = prob_map.reshape((grid_dim, grid_dim)).T
    else:
        prob_map_2d = prob_map.reshape((1, n_nodes))

    # Plotting
    plt.figure(figsize=(8, 7))
    interpolation_method = "nearest"
    im = plt.imshow(
        prob_map_2d,
        origin="lower",
        cmap="plasma",
        interpolation=interpolation_method,
        aspect="auto" if prob_map_2d.shape[0] == 1 else "equal",
    )
    # Contour for separating probability regions, if 2D
    if prob_map_2d.shape[0] > 1 and np.all(np.isfinite(prob_map_2d)):
        try:
            plt.contour(
                prob_map_2d, levels=[0.5], colors="white", linewidths=1, alpha=0.5
            )
        except Exception:
            pass  # Ignore if contour isn't possible (e.g. constant array)

    # Map each compound to its most likely node index
    node_assignments = np.argmax(resp, axis=1)
    node_probs = prob_map.flatten()[node_assignments]
    diff = np.abs(y - node_probs)

    # Marker size reflects misclassification
    base_size = 15
    max_size_extra = 45
    sizes = base_size + (diff * max_size_extra)

    # Improved edge color for visibility
    labels_arr = np.asarray(labels).flatten()
    edges = np.where(labels_arr == 0, "white", "black")

    scatter = plt.scatter(
        x_scaled,
        y_scaled,
        c=labels_arr,
        #        cmap='Set1',
        edgecolors=edges,
        linewidths=0.8,
        s=sizes,
        alpha=0.8,
        zorder=10,
    )

    plt.colorbar(im, label=f"Probability of {class_name}")
    plt.title(f"GTM Class Map: {class_name} ({n_nodes} Nodes)")
    plt.axis("off")
    plt.tight_layout()
    plt.show()

    return prob_map
