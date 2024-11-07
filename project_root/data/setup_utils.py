import re

class Renamer:
    @staticmethod
    def rename_df_columns(df, pattern, replacement):
        df.rename(columns=lambda x: re.sub(pattern, replacement, x), inplace=True)

    @staticmethod
    def rename_sessions_data(sessions, patterns):
        for session in sessions:
            df_keys = session.dfs.fetch_all_data_names()
            for df_type, pattern_info in patterns:
                target_df_keys = [key for key in df_keys if df_type in key]
                for df_key in target_df_keys:            
                    df = session.dfs.get_data(df_key)
                    Renamer.rename_df_columns(df, 
                                              pattern_info.get("pattern", ""), 
                                              pattern_info.get("replacement", ""))
    
    @staticmethod
    def extract_region_number(column_name, letter):
        if letter == 'iso':
            letter = ''
        match = re.search(r'Region(\d+)' + letter, column_name)
        return match.group(1) if match else None

    @staticmethod
    def rename_sessions_fiber_to_brain_region(sessions, frequencies):
        # 
        for session in sessions:
            for letter, freq in frequencies.items():
                phot_df_key = f'phot_{freq}'
                df = session.dfs.get_data(phot_df_key)
                if df is None:
                    continue
                df.rename(columns=lambda x: 
                            session.fiber_to_region.get(Renamer.extract_region_number(x, letter), x),
                            inplace=True)
                    
    @staticmethod
    def debug_df_renames(session):
        df_names = session.dfs.fetch_all_data_names()
        print(f"Debugging session: {session.trial_id}")
        for df_name in df_names:
            df = session.dfs.get_data(df_name)
            if df is not None:
                print(f"Columns in {df_name}: {list(df.columns)}")
            else:
                print(f"{df_name} is None or does not exist.")