# data_loader.py

import numpy as np
import pandas as pd
import h5py
import re
from collections import defaultdict


def load_syllable_data(file_path, include_latent_state=True):
    """
    Load syllable data from an HDF5 file.

    Parameters
    ----------
    file_path : str
        Path to the HDF5 file containing syllable data.
    include_latent_state : bool, optional
        Whether to include 'latent_state' datasets. Default is True.

    Returns
    -------
    dict
        A dictionary where keys are Mouse IDs (as strings) and values are pandas DataFrames containing the syllable data.
    """
    syllable_dfs = {}
    mouse_id_pattern = re.compile(r'^T(\d+)_')

    with h5py.File(file_path, 'r') as file:
        for group_name in file.keys():
            # Extract Mouse ID using regex
            match = mouse_id_pattern.match(group_name)
            if match:
                mouse_id = match.group(1)
            else:
                print(f"Warning: Mouse ID not found in group name '{group_name}'. Skipping this group.")
                continue  # Skip this group

            # Check for duplicate Mouse IDs
            if mouse_id in syllable_dfs:
                raise ValueError(f"Error: Multiple groups found for Mouse ID '{mouse_id}' in the file. Only one group per Mouse ID is expected.")

            group = file[group_name]
            df_list = []

            for dataset_name in group.keys():
                # Exclude 'latent_state' datasets if the flag is False
                if not include_latent_state and dataset_name.startswith('latent_state'):
                    continue

                dataset = group[dataset_name][:]
                
                # Convert 1D datasets to 2D
                if dataset.ndim == 1:
                    dataset = dataset.reshape(-1, 1)
                
                num_cols = dataset.shape[1]
                columns = [f"{dataset_name}_{i}" for i in range(num_cols)]
                
                df = pd.DataFrame(dataset, columns=columns)
                df_list.append(df)

            if df_list:
                # Concatenate all DataFrames horizontally
                concatenated_df = pd.concat(df_list, axis=1)
                syllable_dfs[mouse_id] = concatenated_df
            else:
                print(f"Warning: No datasets found in group '{group_name}'.")

    return syllable_dfs


def load_tracking_data(mouse_ids, tracking_prefix_func):
    """
    Load tracking data from CSV files for each mouse ID.

    Parameters
    ----------
    mouse_ids : list of int
        List of mouse IDs to load tracking data for.
    tracking_prefix_func : function
        Function that takes a mouse ID and returns the corresponding CSV file path.

    Returns
    -------
    dict
        A dictionary where keys are Mouse IDs (as strings) and values are pandas DataFrames containing the tracking data.
    """
    tracking_dfs = {}

    for mouse_id in mouse_ids:
        file_path = tracking_prefix_func(mouse_id)
        try:
            df = pd.read_csv(file_path)
            # Select the first column and the last three columns
            cols_to_keep = df.columns[[0] + list(range(-3, 0))]
            tracking_dfs[str(mouse_id)] = df[cols_to_keep].reset_index(drop=True)
        except FileNotFoundError:
            print(f"Warning: Tracking file for Mouse ID {mouse_id} not found at '{file_path}'. Skipping.")
        except Exception as e:
            print(f"Error loading tracking data for Mouse ID {mouse_id}: {e}")

    return tracking_dfs


def load_fiber_data(mouse_ids, fiber_prefix_func):
    """
    Load fiber photometry data from CSV files for each mouse ID.

    Parameters
    ----------
    mouse_ids : list of int
        List of mouse IDs to load fiber data for.
    fiber_prefix_func : function
        Function that takes a mouse ID and returns the corresponding CSV file path.

    Returns
    -------
    dict
        A dictionary where keys are Mouse IDs (as strings) and values are pandas DataFrames containing the fiber data.
    """
    fiber_dfs = {}

    for mouse_id in mouse_ids:
        file_path = fiber_prefix_func(mouse_id)
        try:
            df = pd.read_csv(file_path)
            # Keep columns that don't start with 'fiber' or 'Fiber' and the last two columns
            columns_to_keep = [col for col in df.columns if not col.lower().startswith('fiber')] + list(df.columns[-2:])
            fiber_dfs[str(mouse_id)] = df[columns_to_keep].reset_index(drop=True)
        except FileNotFoundError:
            print(f"Warning: Fiber file for Mouse ID {mouse_id} not found at '{file_path}'. Skipping.")
        except Exception as e:
            print(f"Error loading fiber data for Mouse ID {mouse_id}: {e}")

    return fiber_dfs


def load_fiber_data_guppy(mouse_ids, fiber_prefix_func):
    """
    Load fiber photometry data from CSV files for each mouse ID.

    Parameters
    ----------
    mouse_ids : list of int
        List of mouse IDs to load fiber data for.
    fiber_prefix_func : function
        Function that takes a mouse ID and returns the corresponding CSV file path.

    Returns
    -------
    dict
        A dictionary where keys are Mouse IDs (as strings) and values are pandas DataFrames containing the fiber data.
    """
    fiber_dfs = {}

    for mouse_id in mouse_ids:
        file_path = fiber_prefix_func(mouse_id)
        try:
            df = pd.read_csv(file_path)
            # Keep columns that don't start with 'fiber' or 'Fiber' and the last two columns
            columns_to_keep = [col for col in df.columns if not col.lower().startswith('fiber')] + list(df.columns[-2:])
            fiber_dfs[str(mouse_id)] = df[columns_to_keep].reset_index(drop=True)
        except FileNotFoundError:
            print(f"Warning: Fiber file for Mouse ID {mouse_id} not found at '{file_path}'. Skipping.")
        except Exception as e:
            print(f"Error loading fiber data for Mouse ID {mouse_id}: {e}")

    return fiber_dfs


def align_and_merge_data(syllable_dfs, tracking_dfs, fiber_dfs):
    """
    Align the lengths of syllable and tracking dataframes and merge them with fiber data.

    Parameters
    ----------
    syllable_dfs : dict
        Dictionary of syllable DataFrames with Mouse IDs as keys.
    tracking_dfs : dict
        Dictionary of tracking DataFrames with Mouse IDs as keys.
    fiber_dfs : dict
        Dictionary of fiber DataFrames with Mouse IDs as keys.

    Returns
    -------
    dict
        A dictionary where keys are Mouse IDs and values are merged pandas DataFrames.
    """
    merged_dfs = {}

    for mouse_id in syllable_dfs.keys():
        if mouse_id not in tracking_dfs:
            print(f"Warning: Tracking data for Mouse ID {mouse_id} not found. Skipping merge for this mouse.")
            continue
        if mouse_id not in fiber_dfs:
            print(f"Warning: Fiber data for Mouse ID {mouse_id} not found. Skipping merge for this mouse.")
            continue

        syllable_df = syllable_dfs[mouse_id]
        tracking_df = tracking_dfs[mouse_id]
        fiber_df = fiber_dfs[mouse_id]

        # Align lengths
        len1 = syllable_df.shape[0]
        len2 = tracking_df.shape[0]
        min_len = min(len1, len2)

        syllable_df_aligned = syllable_df.iloc[:min_len].reset_index(drop=True)
        tracking_df_aligned = tracking_df.iloc[:min_len].reset_index(drop=True)

        # Concatenate tracking data to syllable data
        combined_df = pd.concat([syllable_df_aligned, tracking_df_aligned], axis=1)

        # Merge with fiber data using merge_asof on timestamps
        if 'BonsaiTrackingTimestamp' not in combined_df.columns or 'BonsaiFlyTimestamp' not in fiber_df.columns:
            print(f"Warning: Required timestamp columns missing for Mouse ID {mouse_id}. Skipping merge.")
            continue

        try:
            merged_df = pd.merge_asof(
                combined_df.sort_values('BonsaiTrackingTimestamp'),
                fiber_df.sort_values('BonsaiFlyTimestamp'),
                left_on='BonsaiTrackingTimestamp',
                right_on='BonsaiFlyTimestamp',
                direction='nearest',
                suffixes=('_syllable', '_fiber')
            )

            # Rename timestamp columns for clarity
            merged_df = merged_df.rename(columns={
                'BonsaiTrackingTimestamp': 'timestamp_syllable',
                'BonsaiFlyTimestamp': 'timestamp_fiber'
            })

            # Drop specific columns
            columns_to_drop = ['MetaFlyTimestamp', 'ExtraVar1', 'FrameNo']
            for col in columns_to_drop:
                if col in merged_df.columns:
                    merged_df = merged_df.drop(columns=col)

            # Rename fiber columns
            rename_columns = {}
            for col in merged_df.columns:
                if col == 'Fiber1_ZdFF_scaled':
                    rename_columns[col] = 'VS_470'
                elif col == 'Fiber3_ZdFF_scaled':
                    rename_columns[col] = 'DS_470'
            merged_df = merged_df.rename(columns=rename_columns)

            merged_dfs[mouse_id] = merged_df
        except Exception as e:
            print(f"Error merging data for Mouse ID {mouse_id}: {e}")

    return merged_dfs
