import os
import fnmatch
import re
import pandas as pd
from tqdm import tqdm


class DataContainer:
    def __init__(self, data_type=None):
        self.data = {}
        self.data_type = data_type

    def add_data(self, name, data):
        # Set the data type if not already set
        if self.data_type is None:
            self.data_type = type(data)
        elif not isinstance(data, self.data_type):
            raise TypeError(f"Data must be of type {self.data_type.__name__}")

        self.data[name] = data

    def get_data(self, name):
        return self.data.get(name)

    def remove_data(self, name):
        if name in self.data:
            del self.data[name]

    def fetch_all_data_names(self):
        return list(self.data.keys())

    def clear_data(self):
        self.data.clear()


class Session:
    def __init__(self, chamber_id, trial_dir, session_guide):
        self.trial_dir = trial_dir
        self.trial_id = os.path.basename(trial_dir)
        self.session_guide = session_guide
        
        self.chamber_id = chamber_id.upper()  # Ensure chamber_id is uppercase
        self.setup_id = session_guide.setup_id
        
        self.mouse_id = session_guide.mouse_id
        self.fiber_to_region = self.create_fiber_to_region_dict()
        self.brain_regions = sorted(list(self.fiber_to_region.values()))
        
        # Initialize DataContainer for DataFrame storage
        self.dfs = DataContainer()
        
        self.load_all_data()
        
    def load_data(self, file_pattern, skip_rows=None, use_cols=None, only_header=False):
        file_name = next((f for f in os.listdir(self.trial_dir) if fnmatch.fnmatch(f, file_pattern)), None)
        if file_name:
            file_path = os.path.join(self.trial_dir, file_name)
            if only_header:
                df = pd.read_csv(file_path, nrows=0)
                return df.columns.tolist()
            else:
                return pd.read_csv(file_path, skiprows=skip_rows, usecols=use_cols)
        else:
            return None

    def load_all_data(self):
        # Load photometry data
        #TODO:save memory by filtering out unnecesarry columns (usecol)
        phot_df = self.load_data('photometry_data*.csv', use_cols=self.filter_columns)
        tracking_df = self.load_data('CenterTopCam_TrackData*.csv')
        self.dfs.add_data('phot', phot_df)
        self.dfs.add_data('tracking', tracking_df)

    # create_fiber_dict creates dictionary of all fibers used and their corresponding brainregion
    def create_fiber_to_region_dict(self, fiber_pattern=re.compile(r'fiber(\d+)')):
        # Initialize an empty dictionary
        fiber_to_region_dict = {}
        
        # Iterate over the DataFrame index with enumeration for index and column label
        for idx, col in enumerate(self.session_guide.index):
            # Check if the column matches the fiber pattern and the value is not NaN
            if fiber_pattern.match(col) and pd.notna(self.session_guide[col]):
                # Ensure we don't access an index out of range
                if idx + 1 < len(self.session_guide.index):
                    # Check if the next item is NaN using iloc for safe access
                    if pd.isna(self.session_guide.iloc[idx + 1]):
                        # Extract the fiber number and add it to the dictionary
                        fiber_number = fiber_pattern.match(col).group(1)
                        region_and_side = self.session_guide[col]
                        fiber_to_region_dict[fiber_number] = tuple(region_and_side.split("_"))
        return fiber_to_region_dict

    # to save memory, all of the columns that contain data from unused fibers, are filtered out pre-loading
    def filter_columns(self, col):
        pattern = re.compile(r"[GR](\d+)")
        match = pattern.match(col)
        if match:
            # Check if the captured number is in the keys of fiber_dict
            return match.group(1) in self.fiber_to_region
        else:
            # If the column doesn't match the pattern, include it
            return True
        

# Custom sort key function that extracts the trial number from a directory name
def sort_key_func(dir_name):
    # Find all numbers following 'T' or before a period '.' and return them as a tuple of integers
    numbers = tuple(map(int, re.findall(r'T(\d+)', dir_name)))
    return numbers


# TODO: a bit limited in parameters, recommend additional ones if need be
def load_all_sessions(baseline_dir, first_n_dirs=None, remove_bad_signal_sessions=False):
    # Get a list of all subdirectories within the baseline directory
    subdirs = [d for d in os.listdir(baseline_dir) if os.path.isdir(os.path.join(baseline_dir, d))]

    # Sort the subdirectories based on the trial number
    sorted_subdirs = sorted(subdirs, key=sort_key_func)

    # Join the sorted subdirectories with the baseline path
    trial_dirs = [os.path.join(baseline_dir, sd) for sd in sorted_subdirs]
    if first_n_dirs is None:
        first_n_dirs = len(trial_dirs)

    all_sessions = []

    for trial_dir in tqdm(trial_dirs[:first_n_dirs]):
        trial_id = os.path.basename(trial_dir)
        segments = trial_id.split('_')[1].split('.')  # Assuming the format is like 'T1_23.25.29.e'
        
        for file in os.listdir(trial_dir):
            # Check if file ends with 'trial_guide.xlsx'
            if fnmatch.fnmatch(file, '*trial_guide.xlsx'):
                # Load into DataFrame
                current_trial_guide_df = pd.read_excel(os.path.join(trial_dir, file), nrows=4,
                                                    dtype={"mouse_id": str}, index_col=0, engine='openpyxl')
        for segment, chamber_id in zip(segments, "abcd"):
            if segment == 'e':
                continue  # Skip this session as it's marked as empty
            session_guide = current_trial_guide_df.loc[chamber_id]

            if session_guide.mouse_id != segment:
                raise Exception("The mouse id from the folder names and trial guide do not match")
            
            new_session = Session(chamber_id, trial_dir, session_guide)
            if len(new_session.brain_regions) > 0 or (remove_bad_signal_sessions == False):
                all_sessions.append(new_session)
    return all_sessions