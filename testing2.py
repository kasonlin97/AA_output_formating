import pandas as pd
import os

# 1. Define the absolute paths (Adjust as needed for your Mac)
segment_path = '/Users/admin/S16_AA/S16A2_segment.csv'
cycle_path = '/Users/admin/S16_AA/S16A2_cycle.csv'
breakpoint_path = '/Users/admin/S16_AA/S16A2_breakpoint.csv'

# 2. Load the dataframes
segment = pd.read_csv(segment_path, index_col=False)
cycle = pd.read_csv(cycle_path, index_col=False)[['Cycle', 'Segment', 'Orientation']]
# Renamed to df_bp to avoid conflict with Python's built-in breakpoint()
df_bp = pd.read_csv(breakpoint_path, index_col=False) 

# 3. Merge cycle and segment into df1
df1 = pd.merge(cycle, segment, on='Segment', how='left')

# 4. Prepare df_bp for matching by creating clean coordinates (no +/-)
df_bp['Start_Coordinate'] = df_bp['Start_Pos'].str.rstrip('+-')
df_bp['End_Coordinate'] = df_bp['End_Pos'].str.rstrip('+-')

# 5. Initialize the output columns
df1['5_BP'] = None
df1['3_BP'] = None

# Identify occurrences
cycle_counts = df1['Cycle'].value_counts()

# --- CASE 1: Single Occurrence Cycles ---
single_occurrence_cycles = cycle_counts[cycle_counts == 1].index.tolist()
df1_single = df1[df1['Cycle'].isin(single_occurrence_cycles)].copy()

for idx, row in df1_single.iterrows():
    # Match the start coordinate
    start_match = df_bp[df_bp['Start_Coordinate'] == row['segment_coordinate_start']]
    if not start_match.empty:
        df1_single.at[idx, '5_BP'] = start_match.iloc[0]['Start_Pos']

    # Match the end coordinate
    end_match = df_bp[df_bp['End_Coordinate'] == row['segment_coordinate_end']]
    if not end_match.empty:
        df1_single.at[idx, '3_BP'] = end_match.iloc[0]['End_Pos']

# Update main df1 with Case 1 results
df1.update(df1_single[['5_BP', '3_BP']])


# --- CASE 2: Multiple Occurrence Cycles ---
multiple_occurrence_cycles = cycle_counts[cycle_counts > 1].index.tolist()
df1_multiple = df1[df1['Cycle'].isin(multiple_occurrence_cycles)].copy()

# Pre-calculate first segment start coordinate for each cycle
first_segment_starts = {}
for cycle_num in multiple_occurrence_cycles:
    cycle_rows = df1_multiple[df1_multiple['Cycle'] == cycle_num].sort_values('Segment')
    if not cycle_rows.empty:
        first_segment_starts[cycle_num] = cycle_rows.iloc[0]['segment_coordinate_start']

for idx, row in df1_multiple.iterrows():
    cycle_num = row['Cycle']
    
    if row['Orientation'] == '+':
        # Match segment_coordinate_end with df_bp's Start_Coordinate
        start_match = df_bp[df_bp['Start_Coordinate'] == row['segment_coordinate_end']]
        if not start_match.empty:
            matched_row = start_match.iloc[0]
            # Filling the columns (using the original 'Pos' to keep orientation signs)
            df1.at[idx, '5_BP'] = matched_row['Start_Pos']
            df1.at[idx, '3_BP'] = matched_row['End_Pos']

    elif row['Orientation'] == '-':
        # Orientation '-': Use cycle start and current segment start
        if cycle_num in first_segment_starts:
            df1.at[idx, '5_BP'] = first_segment_starts[cycle_num]
            df1.at[idx, '3_BP'] = row['segment_coordinate_start']

# 6. Final Result
print(df1)
df1.to_csv('/Users/admin/S16_AA/S16A2_final_output.csv', index=False)