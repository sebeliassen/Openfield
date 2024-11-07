# all_mouse_intervals = {}
# all_mouse_intervals_timestamps = {}
from data.guppy_loading import list_folders

import os
import pandas as pd

def find_non_nan_intervals(centroid_coords):
    nan_to_vals = []
    val_to_nans = []

    # Use pandas' vectorized operations for efficiency
    is_nan = centroid_coords.isna()
    transitions = is_nan != is_nan.shift()
    transition_indices = transitions[transitions].index

    for idx in transition_indices:
        if not is_nan[idx]:
            nan_to_vals.append(idx)
        else:
            val_to_nans.append(idx)

    # Handle the start of the series
    if not is_nan.iloc[0]:
        nan_to_vals.insert(0, 0)

    # Handle the end of the series
    if not is_nan.iloc[-1]:
        val_to_nans.append(len(centroid_coords))

    # Pair nan_to_vals with val_to_nans to create intervals
    val_intervals = list(zip(val_to_nans, nan_to_vals))
    val_intervals = [(a, b) for a, b in val_intervals if (b - a) > 4 * 60]
    if len(val_intervals) != 3:
        print(f'warning: length of val_intervals is {len(val_intervals)} and not 3, None returned')
        return None
    else:
        return val_intervals

    
def process_mouse_data_from_path(raw_path):
        all_mouse_intervals = {}
        all_mouse_intervals_timestamps = {}
        #return non_nan_triplet is not None
        for dir_name in list_folders(raw_path):
            directory_path = os.path.join(raw_path, dir_name)
            first_csv_file = next(f for f in os.listdir(directory_path) 
                                if f.endswith('.csv') and f.startswith('CenterTopCam_TrackData'))
            
            cam_df = pd.read_csv(os.path.join(directory_path, first_csv_file))
            trial, mice_ids = dir_name.split('_')
            mice_ids = mice_ids.split('.')

            arena_columns = ['CentroidCoords.Arena1.X', 'CentroidCoords.Arena2.X']
            
            for mouse_id, arena_column in zip(mice_ids, arena_columns):
                if mouse_id != 'e':
                    nan_mask = cam_df[arena_column].isna()
                    non_nan_triplet = find_non_nan_intervals(cam_df[arena_column])
                    
                    if non_nan_triplet:
                        all_mouse_intervals_timestamps[mouse_id] = [
                            (cam_df['Timestamp.FP3002_System'].iloc[a], 
                             cam_df['Timestamp.FP3002_System'].iloc[b])
                            for a, b in non_nan_triplet
                        ]
                        all_mouse_intervals[mouse_id] = non_nan_triplet
        return all_mouse_intervals_timestamps


def nan_interval_to_segments(nan_intervals):
    if len(nan_intervals) == 3 and all(len(x) == 2 for x in nan_intervals):
        # 15 minutes before c to c
        # d to 30 minutes after d
        # f to 30 minutes after f
        (a, b), (c, d), (e, f) = nan_intervals
        return (c - 15*60, c), (d, d + 30*60), (f, f + 30*60)
    else:
        print('nan intervals do not fit proper pattern')
        return None       