import pandas as pd
import os

# 1. Configuration
root_dir = "/Users/admin/MSc_Project/Extracted_ecDNA_Results"
column_headers = [
    'Segment_Label',           # The literal text "Segment"
    'Segment_number',          # The ID (1, 2, 3...)
    'Segment_chromosome', 
    'segment_coordinate_start', 
    'segment_coordinate_end', 
]

all_segments = []

#print("Starting batch processing of segment data...")

# 2. Recursive loop through all subfolders
for root, dirs, files in os.walk(root_dir):
    for filename in files:
        if filename.endswith("_annotated_cycles.txt"):
            file_path = os.path.join(root, filename)
            
            # Derive Amplicon_ID from filename
            amplicon_id = filename.replace("_annotated_cycles.txt", "")
            
            try:
                # 3. Read the file
                df_raw = pd.read_csv(
                    file_path, 
                    sep='\t',  
                    names=column_headers, 
                    comment='L',
                    header=None,
                    on_bad_lines='skip'
                )

                # 4. Filter for rows starting with "Segment"
                df_seg = df_raw[df_raw['Segment_Label'] == 'Segment'].copy()
                
                if df_seg.empty:
                    continue

                # 5. Fix the Float (.0) issue (Your specific logic)
                # A. Handle Chromosome (Keep X/Y, but remove .0 from numbers like 10.0)
                def clean_chr(val):
                    val = str(val).split('.')[0] # Remove decimal if it exists
                    if val.lower() == 'nan' or val == '': return '0'
                    return val
                df_seg['Segment_chromosome'] = df_seg['Segment_chromosome'].apply(clean_chr)

                # B. Handle Numbers (Segment # and coordinates)
                for col in ['Segment_number', 'segment_coordinate_start', 'segment_coordinate_end']:
                    df_seg[col] = pd.to_numeric(df_seg[col], errors='coerce').fillna(0).astype(int).astype(str)

                # 6. Format Coordinates (Chr:Start and Chr:End)
                df_seg['segment_coordinate_start'] = df_seg['Segment_chromosome'] + ":" + df_seg['segment_coordinate_start']
                df_seg['segment_coordinate_end'] = df_seg['Segment_chromosome'] + ":" + df_seg['segment_coordinate_end']

                # 7. Final Selection and Rename
                df_final = df_seg[['Segment_number', 'segment_coordinate_start', 'segment_coordinate_end']].copy()
                df_final.columns = ['Segment', 'segment_coordinate_start', 'segment_coordinate_end']
                
                # Insert the Amplicon_ID at the very beginning
                df_final.insert(0, 'Amplicon_ID', amplicon_id)

                all_segments.append(df_final)
                print(f"Processed Segments: {amplicon_id}")

            except Exception as e:
                print(f"Error processing {filename}: {e}")

# 8. Combine all amplicons into one master CSV
if all_segments:
    master_df = pd.concat(all_segments, ignore_index=True)
    output_path = "/Users/admin/MSc_Project/ecDNA_BP_mapping/TCGA_ecDNA_Segments.csv"
    master_df.to_csv(output_path, index=False)
    
    print(f"\nSuccess! Combined segments for {len(all_segments)} amplicons.")
    print(f"Final file saved to: {output_path}")
    print(master_df.head())
else:
    print("No segment data found to process.")

