import pandas as pd
import os

# 1. Configuration
root_dir = "/Users/admin/MSc_Project/Extracted_ecDNA_Results"
column_headers = [
    'BreakpointEdge', 
    'StartPosition->EndPosition', 
    'PredictedCopyCount', 
    'NumberOfReadPairs', 
    'HomologySizeIfAvailable', 
    'Homology/InsertionSequence'
]

all_results = []

print("Starting batch processing of graph files...")

# 2. Loop through the directory structure
for root, dirs, files in os.walk(root_dir):
    for filename in files:
        if filename.endswith("_graph.txt"):
            file_path = os.path.join(root, filename)
            
            # Extract Amplicon_ID from filename (removes '_graph.txt')
            amplicon_id = filename.replace("_graph.txt", "")
            
            try:
                # 3. Read file to find the header pattern (skip_rows logic)
                skip_rows = 0
                with open(file_path, 'r') as f:
                    for i, line in enumerate(f):
                        if "BreakpointEdge: StartPosition->EndPosition" in line:
                            skip_rows = i + 1
                            break
                
                # 4. Load the data
                df = pd.read_csv(
                    file_path, 
                    sep='\t', 
                    skiprows=skip_rows, 
                    names=column_headers
                )
                
                # If file is empty or header not found, skip
                if df.empty:
                    continue

                # 5. Transformation Logic
                # Split coordinates
                df[['Start_Pos', 'End_Pos']] = df['StartPosition->EndPosition'].str.split('->', expand=True)
                df = df.drop(columns=['StartPosition->EndPosition'])

                # Filter for discordant edges
                filtered_df = df[df['BreakpointEdge'] == 'discordant'].copy()

                # Reorder and Insert Amplicon_ID
                cols_to_order = [
                    'BreakpointEdge', 'Start_Pos', 'End_Pos', 
                    'PredictedCopyCount', 'NumberOfReadPairs', 
                    'HomologySizeIfAvailable', 'Homology/InsertionSequence'
                ]
                df_breakpoint = filtered_df[cols_to_order]
                df_breakpoint.insert(0, 'Amplicon_ID', amplicon_id)

                # Collect the processed dataframe
                all_results.append(df_breakpoint)
                print(f"Processed: {amplicon_id}")

            except Exception as e:
                print(f"Error processing {filename}: {e}")

# 6. Combine all amplicons into one master dataframe
if all_results:
    master_df = pd.concat(all_results, ignore_index=True)
    
    # Save the combined result
    output_path = "/Users/admin/MSc_Project/TCGA_ecDNA_Breakpoints.csv"
    master_df.to_csv(output_path, index=False)
    
    print(f"\nSuccess! Combined {len(all_results)} files into {output_path}")
    print(master_df.head())
else:
    print("No matching files were found or processed.")

