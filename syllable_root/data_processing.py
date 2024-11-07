# data_processing.py

import numpy as np
import pandas as pd
from collections import defaultdict
from scipy import stats


def compute_syllable_counts(merged_dfs):
    """
    Compute syllable initiation counts and time spent proportions across all DataFrames.

    Parameters
    ----------
    merged_dfs : dict
        Dictionary where keys are Mouse IDs and values are merged pandas DataFrames.

    Returns
    -------
    tuple
        initiations_proportions (dict): Proportion of initiations per syllable.
        time_spent_proportions (dict): Proportion of time spent per syllable.
    """
    initiations_proportions = defaultdict(list)
    time_spent_proportions = defaultdict(list)
    all_syllables = set()

    # Gather all unique syllables
    for df in merged_dfs.values():
        all_syllables.update(df['syllable_0'].dropna().unique())

    all_syllables = sorted(all_syllables)

    for mouse_id, df in merged_dfs.items():
        # Ensure syllable_0 is numeric
        df['syllable_0'] = pd.to_numeric(df['syllable_0'], errors='coerce')

        # Compute initiations
        syllable_series = df['syllable_0']
        syllable_shift = syllable_series.shift(1)
        initiations = (syllable_series != syllable_shift) & (~syllable_series.isna())

        # Count initiations per syllable
        num_initiations = syllable_series[initiations].value_counts()
        total_initiations = num_initiations.sum()

        # Calculate initiation proportions
        proportions_initiations = num_initiations / total_initiations if total_initiations > 0 else pd.Series()

        # Count time spent per syllable
        time_spent = syllable_series.value_counts()
        total_time = time_spent.sum()

        # Calculate time spent proportions
        proportions_time_spent = time_spent / total_time if total_time > 0 else pd.Series()

        # Collect proportions for each syllable
        for syllable in all_syllables:
            prop_init = proportions_initiations.get(syllable, 0)
            prop_time = proportions_time_spent.get(syllable, 0)
            initiations_proportions[syllable].append(prop_init)
            time_spent_proportions[syllable].append(prop_time)

    return initiations_proportions, time_spent_proportions


def rank_syllables(initiations_proportions, time_spent_proportions):
    """
    Rank syllables using the Borda count method based on initiation and time spent proportions.

    Parameters
    ----------
    initiations_proportions : dict
        Proportion of initiations per syllable.
    time_spent_proportions : dict
        Proportion of time spent per syllable.

    Returns
    -------
    tuple
        avg_initiations_proportions (dict): Average initiation proportions per syllable.
        avg_time_spent_proportions (dict): Average time spent proportions per syllable.
        final_ranking (list): Final ranked list of syllables based on Borda count.
    """
    avg_initiations_proportions = {}
    avg_time_spent_proportions = {}
    borda_scores = defaultdict(int)

    syllables = sorted(initiations_proportions.keys())
    num_syllables = len(syllables)

    # Compute average proportions
    for syllable in syllables:
        avg_init = np.mean(initiations_proportions[syllable]) if initiations_proportions[syllable] else 0
        avg_time = np.mean(time_spent_proportions[syllable]) if time_spent_proportions[syllable] else 0
        avg_initiations_proportions[syllable] = avg_init
        avg_time_spent_proportions[syllable] = avg_time

    # Sort syllables based on initiation proportions
    sorted_initiations = sorted(avg_initiations_proportions.items(), key=lambda x: x[1], reverse=True)
    # Sort syllables based on time spent proportions
    sorted_time_spent = sorted(avg_time_spent_proportions.items(), key=lambda x: x[1], reverse=True)

    # Assign Borda points
    for rank, (syllable, _) in enumerate(sorted_initiations):
        points = num_syllables - rank
        borda_scores[syllable] += points

    for rank, (syllable, _) in enumerate(sorted_time_spent):
        points = num_syllables - rank
        borda_scores[syllable] += points

    # Final ranking based on Borda scores
    final_ranking = sorted(borda_scores.items(), key=lambda x: x[1], reverse=True)

    return avg_initiations_proportions, avg_time_spent_proportions, final_ranking


def extract_signal_snippets(merged_dfs, parameters):
    """
    Extract and normalize signal snippets around syllable initiations.

    Parameters
    ----------
    merged_dfs : dict
        Dictionary where keys are Mouse IDs and values are merged pandas DataFrames.
    parameters : dict
        Dictionary containing parameters for snippet extraction:
            - m: int, number of frames before initiation
            - n: int, number of frames after initiation
            - normalization_frame: int, frame index corresponding to initiation
            - window_size: int, rolling mean window size
            - min_snippets_required: int, minimum number of snippets required per syllable

    Returns
    -------
    tuple
        syllable_snippets_DS (defaultdict): Normalized DS_470 signal snippets.
        syllable_snippets_VS (defaultdict): Normalized VS_470 signal snippets.
    """
    m = parameters.get('m', 150)
    n = parameters.get('n', 150)
    normalization_frame = parameters.get('normalization_frame', 0)
    window_size = parameters.get('window_size', 15)
    min_snippets_required = parameters.get('min_snippets_required', 10)

    syllable_snippets_DS = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    syllable_snippets_VS = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    time_axis = np.arange(-m, n + 1)

    for mouse_id, df in merged_dfs.items():
        # Get genotype
        genotype = str(df['Genotype_syllable'].iloc[0]) if 'Genotype_syllable' in df.columns else 'Unknown'

        df = df.reset_index(drop=True)

        # Ensure required columns are present
        required_columns = {'DS_470', 'VS_470', 'syllable_0', 'SecFromInjection_fiber'}
        if not required_columns.issubset(df.columns):
            print(f"Warning: Required columns missing in DataFrame for Mouse ID {mouse_id}. Skipping.")
            continue

        # Convert 'syllable_0' and 'SecFromInjection_fiber' to numeric
        df['syllable_0'] = pd.to_numeric(df['syllable_0'], errors='coerce')
        df['SecFromInjection_fiber'] = pd.to_numeric(df['SecFromInjection_fiber'], errors='coerce')

        # Identify syllable initiations
        syllable_shift = df['syllable_0'].shift(1)
        initiations = (df['syllable_0'] != syllable_shift) & (~df['syllable_0'].isna())

        # Define time windows in minutes
        pre_injection_window_start = -25
        pre_injection_window_end = -5
        post_injection_window_start = 10
        post_injection_window_end = 40

        # Iterate over top syllables (0 to 29)
        top_syllables = list(range(40))
        for syllable in top_syllables:
            # Find initiation indices for the current syllable
            syllable_initiations = initiations & (df['syllable_0'] == syllable)
            initiation_indices = syllable_initiations[syllable_initiations].index

            for idx in initiation_indices:
                # Get the 'SecFromInjection_fiber' value at initiation
                sec_from_injection = df.at[idx, 'SecFromInjection_fiber']

                # Determine injection phase
                if pre_injection_window_start <= sec_from_injection <= pre_injection_window_end:
                    injection_phase = 'pre'
                elif post_injection_window_start <= sec_from_injection <= post_injection_window_end:
                    injection_phase = 'post'
                else:
                    continue  # Skip syllable initiations outside the time windows

                start_idx = idx - m
                end_idx = idx + n

                # Check if the window is within DataFrame bounds
                if start_idx >= 0 and end_idx < len(df):
                    # Extract signal snippets
                    snippet_DS = df['DS_470'].iloc[start_idx:end_idx + 1].values
                    snippet_VS = df['VS_470'].iloc[start_idx:end_idx + 1].values

                    # Apply rolling mean
                    snippet_DS_rolled = pd.Series(snippet_DS).rolling(window=window_size, min_periods=1).mean().values
                    snippet_VS_rolled = pd.Series(snippet_VS).rolling(window=window_size, min_periods=1).mean().values

                    # Normalize the signals
                    normalization_value_DS = snippet_DS_rolled[normalization_frame]
                    normalization_value_VS = snippet_VS_rolled[normalization_frame]

                    snippet_DS_normalized = snippet_DS_rolled - normalization_value_DS
                    snippet_VS_normalized = snippet_VS_rolled - normalization_value_VS

                    # Append the normalized snippets to the respective lists
                    syllable_snippets_DS[genotype][injection_phase][syllable].append(snippet_DS_normalized)
                    syllable_snippets_VS[genotype][injection_phase][syllable].append(snippet_VS_normalized)
                else:
                    continue  # Skip if the window is incomplete

    # After collecting all the snippets:
    for genotype in syllable_snippets_DS:
        for injection_phase in syllable_snippets_DS[genotype]:
            for syllable in syllable_snippets_DS[genotype][injection_phase]:
                syllable_snippets_DS[genotype][injection_phase][syllable] \
                    = np.vstack(syllable_snippets_DS[genotype][injection_phase][syllable])
                syllable_snippets_VS[genotype][injection_phase][syllable] \
                    = np.vstack(syllable_snippets_VS[genotype][injection_phase][syllable])

    return syllable_snippets_DS, syllable_snippets_VS


def process_snippets(data):
    """
    Process signal snippets to calculate mean, confidence interval, and sample size.

    Parameters
    ----------
    snippets : list of np.ndarray
        List of signal snippets.

    Returns
    -------
    tuple
        mean (np.ndarray): Mean signal across snippets.
        ci (np.ndarray): Confidence interval for the mean.
        n (int): Number of snippets.
    """
    if data.ndim == 2 and data.size > 0:
        mean = np.mean(data, axis=0)
        sem = stats.sem(data, axis=0, nan_policy='omit')
        t_value = stats.t.ppf(0.975, df=len(data) - 1) if len(data) > 1 else 0
        ci = sem * t_value
        n = len(data)
        return mean, ci, n
    else:
        return None, None, 0
