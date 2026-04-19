import pandas as pd
import os

# 1. Configuration
root_dir = "/Users/admin/MSc_Project/Extracted_ecDNA_Results"
column_headers = [
    'Cycle', 
    'Copy_count', 
    'Length', 
    'IsCyclicPath', 
    'CycleClass', 
    'Segments', 
]

all_annotated_results = []

print("Starting batch processing of annotated cycles...")

# 2. Loop through the directory structure
for root, dirs, files in os.walk(root_dir):
    for filename in files:
        if filename.endswith("_annotated_cycles.txt"):
            file_path = os.path.join(root, filename)
            
            # Extract Amplicon_ID (removes '_annotated_cycles.txt')
            amplicon_id = filename.replace("_annotated_cycles.txt", "")
            
            try:
                # 3. Find where the cycle data starts
                skiprows = -1
                with open(file_path, 'r') as f:
                    for i, line in enumerate(f):
                        if line.startswith('Cycle='):
                            skiprows = i
                            break
                
                if skiprows == -1:
                    continue # Skip files with no data

                # 4. Load the data
                df_cycle = pd.read_csv(
                    file_path, 
                    sep=';', 
                    skiprows=skiprows,  
                    names=column_headers,
                    header=None
                )

                # 5. Clean data: Use .astype(str) to prevent AttributeErrors on NaNs
                for col in df_cycle.columns:
                    df_cycle[col] = df_cycle[col].str.split('=').str[-1]

                # 6. Explode Segments
                # Filter out 'nan' segments before exploding
                df_cycle = df_cycle[df_cycle['Segments'] != 'nan']
                df_cycle['Segments'] = df_cycle['Segments'].str.split(',')
                df_exploded = df_cycle.explode('Segments')

                # 7. Filter for ecDNA-like paths
                # .copy() is crucial here to avoid SettingWithCopy warnings
                filtered_df = df_exploded[df_exploded['CycleClass'] == 'ecDNA-like'].copy()

                if filtered_df.empty:
                    continue

                # 8. Extract Segment and Orientation
                filtered_df[['Segment', 'Orientation']] = filtered_df['Segments'].str.extract(r'(\d+)([+-])')
                filtered_df = filtered_df.drop(columns=['Segments'])

                # 9. Format and Insert ID
                filtered_df['Segment'] = pd.to_numeric(filtered_df['Segment'], errors='coerce')
                filtered_df['Cycle'] = pd.to_numeric(filtered_df['Cycle'], errors='coerce')
                filtered_df['Copy_count'] = pd.to_numeric(filtered_df['Copy_count'], errors='coerce')
                
                # Insert Amplicon_ID at the beginning (column 0)
                filtered_df.insert(0, 'Amplicon_ID', amplicon_id)

                all_annotated_results.append(filtered_df)
                print(f"Processed: {amplicon_id}")

            except Exception as e:
                print(f"Error processing {filename}: {e}")

# 10. Combine all files into one master dataframe
if all_annotated_results:
    master_df = pd.concat(all_annotated_results, ignore_index=True)
    
    # Save the master table
    output_path = "/Users/admin/MSc_Project/TCGA_ecDNA_AnnotatedCycles.csv"
    master_df.to_csv(output_path, index=False)
    
    print(f"\nSuccess! Combined {len(all_annotated_results)} files.")
    print(f"Final file saved to: {output_path}")
else:
    print("No matching cycle data was found.")





