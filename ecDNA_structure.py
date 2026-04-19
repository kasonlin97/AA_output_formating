import pandas as pd
import os
from tqdm import tqdm

# 1. SETUP PATHS
folder_path = '/Users/admin/MSc_Project/'
# Using the Master files created from the previous aggregation steps
segment_file = os.path.join(folder_path, 'TCGA_ecDNA_Segments.csv')
cycle_file = os.path.join(folder_path, 'TCGA_ecDNA_AnnotatedCycles.csv')
breakpoint_file = os.path.join(folder_path, 'TCGA_ecDNA_Breakpoints.csv')

# 2. LOAD DATA
print("Loading master dataframes...")
df_seg = pd.read_csv(segment_file)

# Ensure columns match your master cycle file
df_cycle = pd.read_csv(cycle_file)[['Amplicon_ID', 'Cycle', 'Segment', 'Orientation']]
df_bp = pd.read_csv(breakpoint_file)

# 3. PREPARE BREAKPOINT LOOKUP (Dictionary for Speed)
# Strip the +/- to get the raw coordinate for matching
df_bp['Start_Coord'] = df_bp['Start_Pos'].str.rstrip('+-')
df_bp['End_Coord'] = df_bp['End_Pos'].str.rstrip('+-')

# Create a lookup key combining Amplicon and Coordinates to ensure sample-specific matching
# Key: (Amplicon_ID, Coord1, Coord2) -> (Start_Pos_with_sign, End_Pos_with_sign)
bp_lookup = {}
for _, row in df_bp.iterrows():
    aid = row['Amplicon_ID']
    c1, c2 = row['Start_Coord'], row['End_Coord']
    s1, s2 = row['Start_Pos'], row['End_Pos']
    
    # Store both directions so we can find the match regardless of order
    bp_lookup[(aid, c1, c2)] = (s1, s2)
    bp_lookup[(aid, c2, c1)] = (s2, s1)

# 4. MERGE DATA (Must include Amplicon_ID to keep segments unique)
df = pd.merge(df_cycle, df_seg, on=['Amplicon_ID', 'Segment'], how='left')

# 5. PROCESS CYCLES
final_results = []
# Group by both ID and Cycle to process each ecDNA structure individually
grouped = df.groupby(['Amplicon_ID', 'Cycle'])

print(f"Analyzing junctions for {len(grouped)} ecDNA cycles...")

for (amp_id, cycle_no), cycle_df in tqdm(grouped, desc="Processing Junctions", unit="cycle"):
    cycle_df = cycle_df.reset_index(drop=True)
    n = len(cycle_df)
    
    b5_vals = [] # Exit Breakpoint
    b3_vals = [] # Entry Breakpoint
    
    for i in range(n):
        curr = cycle_df.iloc[i]
        next_s = cycle_df.iloc[(i + 1) % n] # Wraps around to the start (circular)
        
        # Determine coordinates based on orientation
        # (+) Exit at End, Entry at Start | (-) Exit at Start, Entry at End
        exit_pt = str(curr['segment_coordinate_end']) if curr['Orientation'] == '+' else str(curr['segment_coordinate_start'])
        entry_pt = str(next_s['segment_coordinate_start']) if next_s['Orientation'] == '+' else str(next_s['segment_coordinate_end'])
        
        # Lookup the breakpoint pair
        # We use .get() to return (None, None) if no breakpoint exists at that junction
        b5, b3 = bp_lookup.get((amp_id, exit_pt, entry_pt), (None, None))
        
        b5_vals.append(b5)
        b3_vals.append(b3)
        
    cycle_df['5_BP_Junction'] = b5_vals
    cycle_df['3_BP_Junction'] = b3_vals
    final_results.append(cycle_df)

# 6. FINAL ASSEMBLY
if final_results:
    df_final = pd.concat(final_results).reset_index(drop=True)
    df_final.fillna('None', inplace=True)

    # 7. SAVE RESULT
    output_path = os.path.join(folder_path, 'ecDNA_BP_annotation.csv')
    df_final.to_csv(output_path, index=False)
    print(f"\nProcessing Complete.")
    print(f"Total rows in output: {len(df_final)}")
    print(f"Results saved to: {output_path}")
else:
    print("No data processed.")

                   
            
    
