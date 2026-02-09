
import pandas as pd
import os 

file_path = os.path.join("/Users/kason/MSc_ACS_Thesis/S16_amplicon/S16A1", "sample_16_amplicon1_graph (1).txt")

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

df_breakpoint[['Start_Position', 'End_Position']] = df_breakpoint['StartPosition->EndPosition'].str.split('->', expand=True)

#Split Start & End Positions columns into number + orientation 
df_breakpoint[['Start_Coordinate', 'Start_Position_Orientation']] = df_breakpoint['Start_Position'].str.extract(r'(.*?\d+)([+-])')
df_breakpoint[['End_Coordinate', 'End_Position_Orientation']] = df_breakpoint['End_Position'].str.extract(r'(.*?\d+)([+-])')
#print(df_breakpoint_reordered[['Start_Coordinate', 'Start_Position_Orientation']])

#Drop column start_position afterwards for cleaner output 
df_breakpoint = df_breakpoint.drop(columns=['Start_Position'])
df_breakpoint= df_breakpoint.drop(columns=['End_Position'])

#FIlter "BreakpointEdge" column 
filtered_df_breakpoint = df_breakpoint[df_breakpoint['BreakpointEdge'] == 'discordant'].copy()

#print(filtered_df_breakpoint.shape)

df_breakpoint_reordered = filtered_df_breakpoint.iloc[:, [0, 6, 7, 8, 9, 1, 2, 3, 4, 5]]

print(df_breakpoint_reordered)
df_breakpoint_reordered.to_csv('/Users/kason/MSc_ACS_Thesis/S16_amplicon/S16A1/S16A1_breakpoint.csv', index=False)

