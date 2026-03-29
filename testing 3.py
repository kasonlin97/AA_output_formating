import pandas as pd
import os

# 1. SETUP PATHS
folder_path = '/Users/admin/S16_AA/'
segment_file = os.path.join(folder_path, 'S16A2_segment.csv')
cycle_file = os.path.join(folder_path, 'S16A2_cycle.csv')
breakpoint_file = os.path.join(folder_path, 'S16A2_breakpoint.csv')

# 2. LOAD DATA
df_seg = pd.read_csv(segment_file)
df_cycle = pd.read_csv(cycle_file)[['Cycle', 'Segment', 'Orientation']]
df_bp = pd.read_csv(breakpoint_file)

# 3. PREPARE BREAKPOINT LOOKUP
df_bp['Start_Coord'] = df_bp['Start_Pos'].str.rstrip('+-')
df_bp['End_Coord'] = df_bp['End_Pos'].str.rstrip('+-')

def fetch_bp_pair(exit_coord, entry_coord, bp_df):
    """Matches the exit of current segment to the entry of the next."""
    match = bp_df[
        ((bp_df['Start_Coord'] == exit_coord) & (bp_df['End_Coord'] == entry_coord)) |
        ((bp_df['Start_Coord'] == entry_coord) & (bp_df['End_Coord'] == exit_coord))
    ]
    if not match.empty:
        row = match.iloc[0]
        # Orient the BP string so the first one matches our current segment's exit
        if row['Start_Coord'] == exit_coord:
            return row['Start_Pos'], row['End_Pos']
        else:
            return row['End_Pos'], row['Start_Pos']
    return None, None

# 4. MERGE SEGMENT DATA
df = pd.merge(df_cycle, df_seg, on='Segment', how='left')

# 5. PROCESS CYCLES WITH JUNCTION LOGIC
final_results = []

for cycle_id in df['Cycle'].unique():
    cycle_df = df[df['Cycle'] == cycle_id].copy().reset_index(drop=True)
    n = len(cycle_df)
    
    b5_vals = []
    b3_vals = []
    
    for i in range(n):
        curr = cycle_df.iloc[i]
        # Circular logic: last segment connects back to the first
        next_s = cycle_df.iloc[(i + 1) % n]
        
        # --- LOGIC FOR JUNCTION BETWEEN CURRENT AND NEXT ---
        # If current is (+), exit is 'segment_coordinate_end'
        # If current is (-), exit is 'segment_coordinate_start'
        exit_pt = curr['segment_coordinate_end'] if curr['Orientation'] == '+' else curr['segment_coordinate_start']
        
        # If next is (+), entry is 'segment_coordinate_start'
        # If next is (-), entry is 'segment_coordinate_end'
        entry_pt = next_s['segment_coordinate_start'] if next_s['Orientation'] == '+' else next_s['segment_coordinate_end']
        
        # Find the breakpoint bridge for this junction
        b5, b3 = fetch_bp_pair(exit_pt, entry_pt, df_bp)
        
        b5_vals.append(b5)
        b3_vals.append(b3)
        
    cycle_df['5_BP'] = b5_vals
    cycle_df['3_BP'] = b3_vals
    final_results.append(cycle_df)

# 6. FINAL ASSEMBLY
df_final = pd.concat(final_results).reset_index(drop=True)
df_final.fillna('None', inplace=True)

# 7. SAVE
print(df_final.to_string())
output_file = os.path.join(folder_path, 'S16A2_complete_cases.csv')
df_final.to_csv(output_file, index=False)