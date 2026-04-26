import pandas as pd
import os
from tqdm import tqdm
import re

# 1. SETUP PATHS
folder_path = '/Users/admin/MSc_Project/ecDNA_BP_mapping/'
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

    print(f"\nProcessing Complete.")
    print(f"Total rows in output: {len(df_final)}")
else:
    print("No data processed.")


#================Merger with oncogenes, copy number & filter flag=================

df_tcga = pd.read_csv("/Users/admin/MSc_Project/Converted_Results/TCGA_ecDNA_Analysis_result_table.csv")

# 3. Clean 'Feature ID' to extract the part that matches 'Amplicon_ID'
# This regex extracts everything from the starft up to the word 'amplicon' followed by digits
# Example: '..._amplicon2_ecDNA_1' becomes '..._amplicon2'
df_tcga['Cleaned_Match_ID'] = df_tcga['Feature ID'].str.extract(r'^(.*_amplicon\d+)')

# 4. Prepare the TCGA subset
# We filter for only the columns you requested + our new match key
cols_to_keep = ['Cleaned_Match_ID', 'Oncogenes', 'Feature median copy number', 'Filter flag']
df_tcga_subset = df_tcga[cols_to_keep]

# Drop rows that didn't match the amplicon pattern and remove duplicates
# (In case multiple feature types exist for the same amplicon)
df_tcga_subset = df_tcga_subset.dropna(subset=['Cleaned_Match_ID']).drop_duplicates(subset=['Cleaned_Match_ID'])

# 5. Merge the dataframes
df_merged = pd.merge(
    df_final,
    df_tcga_subset,
    left_on='Amplicon_ID',
    right_on='Cleaned_Match_ID',
    how='left'
)

# 6. Cleanup and Save
df_merged = df_merged.drop(columns=['Cleaned_Match_ID'])
#df_merged.to_csv("/Users/admin/BP.csv", index=False)

#print(df_merged[['Amplicon_ID', 'Oncogenes', 'Feature median copy number', 'Filter flag']].head())


#===========================Liftover coordinates=================================


#==================================================CONVERT TO BED===================================== 

# Define the columns that contain coordinates
coord_columns = [
    'segment_coordinate_start', 
    'segment_coordinate_end', 
    '5_BP_Junction', 
    '3_BP_Junction'
]

def parse_to_bed(val):
    if pd.isna(val) or str(val).strip() == "": return None
    match = re.match(r'(\w+):(\d+)', str(val))
    if match:
        chrom_num = match.group(1) # e.g. "10"
        pos = int(match.group(2))
        # BED needs 'chr', but we'll save the NAME without it to match your CSV style
        return [f"chr{chrom_num}", pos - 1, pos, f"{chrom_num}:{pos}"] 
    return None

# Deduplciate to collect all unique coordinates
extracted_data = []
for col in coord_columns:
    # Apply parsing and drop None results
    results = df_merged[col].dropna().apply(parse_to_bed)
    extracted_data.extend([r for r in results if r is not None])

# Create a DataFrame and keep only unique chromosome/position pairs
bed_df = pd.DataFrame(extracted_data, columns=['chrom', 'start', 'end', 'name'])
bed_df = bed_df.drop_duplicates(subset=['chrom', 'start', 'end'])

#print(bed_df)

# Save as a tab-separated file with no header (Standard BED)
bed_df.to_csv('/Users/admin/MSc_Project/ecDNA_BP_mapping/hg19_coords.bed', sep='\t', index=False, header=False)




#============Remapping liftover bed to csv===============================

# 1. Build mapping (This stays mostly same, but names are now clean)
lifted = pd.read_csv('/Users/admin/MSc_Project/ecDNA_BP_mapping/hg38_coords.bed', sep='\t', names=['chrom', 'start', 'end', 'hg19_id'])
mapping = {}
for _, row in lifted.iterrows():
    # Store mapping as "10:103336109" -> "10:newpos"
    mapping[row['hg19_id']] = f"{row['chrom'].replace('chr', '')}:{row['end']}"

# 2. Update function
def update_coord(val, mapping_dict):
    if pd.isna(val): return val
    match = re.match(r'(\w+):(\d+)([+-]?)', str(val))
    if match:
        chrom, pos, suffix = match.group(1), match.group(2), match.group(3)
        # Create key: "10:103336109" (No 'chr', no suffix)
        lookup_key = f"{chrom}:{pos}" 
        
        if lookup_key in mapping_dict:
            # Return "10:newpos" + "-"
            return f"{mapping_dict[lookup_key]}{suffix}"
            
    return val # Returns original if not found in liftOver

# 3. Apply to original dataframe
cols_to_fix = ['segment_coordinate_start', 'segment_coordinate_end', '5_BP_Junction', '3_BP_Junction']

for col in cols_to_fix:
    df_merged[col] = df_merged[col].apply(lambda x: update_coord(x, mapping))

df_merged.to_csv('/Users/admin/MSc_Project/ecDNA_BP_mapping/ecDNA_structure.csv', index=False)
print(df_merged)
                   
            
    
