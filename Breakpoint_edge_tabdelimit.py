
import pandas as pd
import os 

file_path = os.path.join("/Users/admin/S16_AA", "sample_16_amplicon2_graph.txt")

column_headers = [
    'BreakpointEdge', 
    'StartPosition->EndPosition', 
    'PredictedCopyCount', 
    'NumberOfReadPairs', 
    'HomologySizeIfAvailable', 
    'Homology/InsertionSequence'
]

# Read the file to find where the pattern starts
with open(file_path, 'r') as f:
    lines = f.readlines()

# Find the line with your pattern
skip_rows = 0
for i, line in enumerate(lines):
    if line.startswith("BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence"):
        skip_rows = i + 1  # Skip header row too
        break

df_breakpoint = pd.read_csv(
    file_path, 
    sep='\t', 
    skiprows=skip_rows, 
    names=column_headers
)

# 3. Create the combined Start and End columns by splitting the original column
# This keeps the orientation (+/-) attached to the coordinate
df_breakpoint[['Start_Pos', 'End_Pos']] = df_breakpoint['StartPosition->EndPosition'].str.split('->', expand=True)

# 4. Drop the original combined column
df_breakpoint = df_breakpoint.drop(columns=['StartPosition->EndPosition'])

# 5. Filter for "discordant" edges
filtered_df = df_breakpoint[df_breakpoint['BreakpointEdge'] == 'discordant'].copy()

# 6. Reorder columns (using names is safer than using iloc indices)
# This puts Start and End right after the BreakpointEdge
cols_to_order = [
    'BreakpointEdge', 
    'Start_Pos', 
    'End_Pos', 
    'PredictedCopyCount', 
    'NumberOfReadPairs', 
    'HomologySizeIfAvailable', 
    'Homology/InsertionSequence'
]
df_final = filtered_df[cols_to_order]

# Print and Save
print(df_final)
df_final.to_csv('/Users/admin/S16_AA/S16A2_breakpoint.csv', index=False)