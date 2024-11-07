import os
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import matplotlib.ticker as ticker

def list_folders(directory):
    # Get a list of all entries in the directory
    entries = os.listdir(directory)
    
    # Filter the list to include only directories
    folders = [entry for entry in entries if os.path.isdir(os.path.join(directory, entry))]
    folders = sorted(folders)
    return folders

def load_hdf5_file(file_path):
    """
    Load all datasets from an HDF5 file into a single DataFrame.
    
    Parameters:
    - file_path: str, path to the HDF5 file
    
    Returns:
    - df: pandas DataFrame containing all datasets as columns
    """
    data_dict = {}
    
    def extract_datasets(name, obj):
        if isinstance(obj, h5py.Dataset):
            key = name.decode() if isinstance(name, bytes) else name
            data_dict[key] = np.array(obj).flatten() if len(obj.shape) > 1 else np.array(obj)
    
    with h5py.File(file_path, 'r') as file:
        file.visititems(extract_datasets)
    
    # Create DataFrame from dictionary, padding with NaNs for unequal lengths
    max_len = max(len(v) for v in data_dict.values())
    df = pd.DataFrame({k: pd.Series(v).reindex(range(max_len), fill_value=np.nan) for k, v in data_dict.items()})
    
    return df

def load_h5_file(file_path, key):
    """
    Load a specific dataset from an H5 file into a pandas DataFrame.
    
    Parameters:
    - file_path: str, path to the H5 file
    - key: str, key of the dataset to load
    
    Returns:
    - df: pandas DataFrame containing the dataset
    """
    df = pd.read_hdf(file_path, key=key)
    return df


def list_files_and_load_data(directory_path):
    """
    List all HDF5 and H5 files in a directory, load their data into dictionaries and DataFrames,
    and save all of it in an "all_data" dictionary.
    
    Parameters:
    - directory_path: str, path to the directory containing HDF5 and H5 files
    
    Returns:
    - all_data: dict, a dictionary containing data from all files
    """
    all_data = {}

    for filename in os.listdir(directory_path):
        file_path = os.path.join(directory_path, filename)
        
        if filename.endswith(".hdf5"):
            #print(f"Processing HDF5 file: {filename}")
            data_dict = load_hdf5_file(file_path)
            all_data[filename] = data_dict
        
        elif filename.endswith(".h5"):
            #print(f"Processing H5 file: {filename}")
            with h5py.File(file_path, 'r') as file:
                keys = list(file.keys())
                for key in keys:
                    df = load_h5_file(file_path, key)
                    all_data[f"{filename}_{key}"] = df
    
    return all_data


def create_all_z_score_dfs(guppy_path):
    all_z_score_dfs = {}
    all_timestamps = {}

    for folder in list_folders(guppy_path):
        dir_guppy_name = f'{folder}/{folder}_output_1'
        directory_path = guppy_path + '/' + dir_guppy_name
        all_data = list_files_and_load_data(directory_path)
        for k, v in all_data.items():
            if k.startswith('z_score'):
                z_score_fname = k.split(sep='.')[0].split(sep='_')[-1]
                #match = re.match(r"(\d+)(dls|vs|dms)", z_score_fname)
                match = re.match(r"(\d+)(LH|mPFC)", z_score_fname)
                if match:
                    mouse_id, brain_reg = match.groups()
                    all_z_score_dfs[(mouse_id, brain_reg)] = v
            elif isinstance(v, pd.DataFrame) and 'timestamps' in v.columns:
                mouse_id = k.split('_')[1]
                all_timestamps[mouse_id] = v['timestamps']
        #z_score_dfs = {k: v for k, v in all_data.items() if k.startswith('z_score')}

    for k, v in all_z_score_dfs.items():
        mouse_id = k[0]
        v['timestamps'] = all_timestamps[mouse_id]

    for k, v in all_z_score_dfs.items():
    # Remove rows where the 'data' column is NaN
        all_z_score_dfs[k] = v.dropna(subset=['data'])
    
    return all_z_score_dfs