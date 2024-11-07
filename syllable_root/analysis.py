# analysis.py

import numpy as np
import glob
import pickle
from scipy import stats
from scipy.stats import spearmanr
from skbio.stats.distance import mantel
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict
import os
import matplotlib.pyplot as plt


def load_behavior_distances(pickle_pattern):
    """
    Load all behavior distance matrices matching the given glob pattern.

    Parameters
    ----------
    pickle_pattern : str
        Glob pattern to match pickle files.

    Returns
    -------
    dict
        Dictionary mapping filename to loaded data. Each entry contains 'distances' and 'syllable_ixs'.
    """
    pickle_files = glob.glob(pickle_pattern)
    behavior_data = {}
    print(f"Loading behavior distance matrices with pattern: {pickle_pattern}")
    print(f"Number of files found: {len(pickle_files)}")
    for file in pickle_files:
        print('you found something!')
        with open(file, 'rb') as f:
            try:
                data = pickle.load(f)
                behavior_data[file] = data  # data should contain 'distances' and 'syllable_ixs'
                print(f"Loaded file: {file}")
            except Exception as e:
                print(f"Failed to load {file}: {e}")
    return behavior_data


def compute_photometry_distance_matrix(syllable_snippets, signal_type='DS', top_syllables=None, min_snippets=10):
    """
    Compute the photometry distance matrix based on mean signatures per syllable.

    Parameters
    ----------
    syllable_snippets : defaultdict
        Nested dictionary containing DS_470 or VS_470 signal snippets.
    signal_type : str, optional
        Type of photometry signal to use ('DS' or 'VS'). Default is 'DS'.
    top_syllables : list, optional
        List of top syllables to include in the analysis. If None, include all.
    min_snippets : int, optional
        Minimum number of snippets required to include a syllable. Default is 10.

    Returns
    -------
    tuple
        photometry_distances (ndarray): Photometry distance matrix.
        syllables_photometry (list): List of syllables included in the distance matrix.
    """
    photometry_signatures = {}

    for genotype in syllable_snippets:
        for injection_phase in syllable_snippets[genotype]:
            for syllable, snippets in syllable_snippets[genotype][injection_phase].items():
                # Include only top syllables and those with at least 'min_snippets' snippets
                if top_syllables is not None and syllable not in top_syllables:
                    continue
                if len(snippets) < min_snippets:
                    continue  # Exclude syllables with fewer than 'min_snippets' snippets

                if syllable not in photometry_signatures:
                    photometry_signatures[syllable] = []
                photometry_signatures[syllable].extend(snippets)

    # Compute mean signature per syllable
    syllable_means = {s: np.mean(snippets, axis=0) for s, snippets in photometry_signatures.items() if snippets}

    # Sort syllables
    syllables_photometry = sorted(syllable_means.keys())

    if not syllables_photometry:
        print("No syllables meet the criteria for photometry distance matrix computation.")
        return None, []

    # Create feature matrix
    photometry_matrix = np.array([syllable_means[s] for s in syllables_photometry])

    # Standardize
    scaler = StandardScaler()
    photometry_matrix_std = scaler.fit_transform(photometry_matrix)

    # Compute distance matrix using cosine distance
    photometry_distances = squareform(pdist(photometry_matrix_std, metric='cosine'))

    print(f"Computed {signal_type} photometry distance matrix with shape: {photometry_distances.shape}")
    return photometry_distances, syllables_photometry


def compare_distance_matrices(behavior_data, photometry_distances, syllables_photometry, top_syllables=None, min_snippets=10, output_file=None):
    """
    Compare each behavior distance matrix with the photometry distance matrix.

    Parameters
    ----------
    behavior_data : dict
        Dictionary mapping filenames to behavior distance data.
    photometry_distances : ndarray
        Photometry distance matrix.
    syllables_photometry : list
        List of syllable indices for photometry data.
    top_syllables : list, optional
        List of top syllables to include in the comparison. If None, include all.
    min_snippets : int, optional
        Minimum number of snippets required to include a syllable. Default is 10.
    output_file : file object, optional
        File object to write the results. Default is None.

    Returns
    -------
    None
    """
    if photometry_distances is None or not syllables_photometry:
        print("Photometry distance matrix is empty. Skipping comparison.")
        return

    for file, data in behavior_data.items():
        behavior_distances = data.get('distances')
        syllables_behavior = list(data.get('syllable_ixs', []))

        if not isinstance(behavior_distances, np.ndarray):
            print(f"Warning: 'distances' in {file} is not a numpy array. Skipping.")
            continue

        # Include only top syllables
        if top_syllables is not None:
            syllables_behavior = [s for s in syllables_behavior if s in top_syllables]

        # Find common syllables
        common_syllables = sorted(set(syllables_behavior).intersection(syllables_photometry))
        if len(common_syllables) < 2:
            message = f"Not enough common syllables between {file} and photometry data. Skipping.\n"
            print(message.strip())
            if output_file:
                output_file.write(message)
            continue

        # Get indices for common syllables
        behavior_indices = [syllables_behavior.index(s) for s in common_syllables]
        photometry_indices = [syllables_photometry.index(s) for s in common_syllables]

        # Extract sub-distance matrices
        behavior_sub = behavior_distances[np.ix_(behavior_indices, behavior_indices)]
        photometry_sub = photometry_distances[np.ix_(photometry_indices, photometry_indices)]

        # Flatten upper triangles
        triu = np.triu_indices(len(common_syllables), k=1)
        behavior_vec = behavior_sub[triu]
        photometry_vec = photometry_sub[triu]

        # Compute Spearman's rank correlation
        corr, p = spearmanr(behavior_vec, photometry_vec)

        # Compute Mantel Test
        try:
            mantel_corr, mantel_p, _ = mantel(
                behavior_sub, photometry_sub, method='spearman', permutations=999
            )
        except Exception as e:
            mantel_corr, mantel_p = np.nan, np.nan
            print(f"  Mantel test failed for {file}: {e}")

        # Prepare result strings
        result_str = (
            f"\nComparison for {os.path.basename(file)}:\n"
            f"  Spearman correlation coefficient: {corr:.3f}\n"
            f"  P-value: {p:.3e}\n"
            f"  Mantel test correlation coefficient: {mantel_corr:.3f}\n"
            f"  Mantel test P-value: {mantel_p:.3e}\n"
        )

        # Print and write the results
        print(result_str)
        if output_file:
            output_file.write(result_str)

        # Optional: Plotting
        plt.figure(figsize=(6, 6))
        plt.scatter(behavior_vec, photometry_vec, alpha=0.6)
        plt.xlabel('Behavior Distances')
        plt.ylabel('Photometry Distances')
        plt.title(f'{os.path.basename(file)}\nSpearman r={corr:.2f}, p={p:.3f}')
        plt.grid(True)
        plt.show()
