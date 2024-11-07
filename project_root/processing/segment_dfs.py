import numpy as np

def create_segment_dfs(segment_names, all_z_score_dfs, all_mouse_segments):
    """
    Create segment DataFrames based on provided segment names.

    Args:
    segment_names (list): List of segment names (e.g., ['baseline', 'nicotine', 'antagonist'])
    all_z_score_dfs (dict): Dictionary of z-score DataFrames
    all_mouse_segments (dict): Dictionary of mouse segments

    Returns:
    dict: Dictionary of segment DataFrames
    """
    # Initialize a dictionary to store segment DataFrames
    segment_dfs = {name: {} for name in segment_names}

    for (mouse_id, brain_region), df in all_z_score_dfs.items():
        intervals = all_mouse_segments.get(mouse_id)
        if intervals is None:
            continue

        # Extract start and end times
        start_times = [interval[0] for interval in intervals]
        end_times = [interval[1] for interval in intervals]

        # Use searchsorted to find the indices
        start_indices = np.searchsorted(df['timestamps'], start_times)
        end_indices = np.searchsorted(df['timestamps'], end_times)

        # Split the dataframe and assign to corresponding segment
        for name, start, end in zip(segment_names, start_indices, end_indices):
            segment_dfs[name][(mouse_id, brain_region)] = df.iloc[start:end]

    return segment_dfs