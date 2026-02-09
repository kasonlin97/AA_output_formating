
import pandas as pd 

segment = pd.read_csv('segment.csv', index_col=False)
cycle = pd.read_csv('cycle.csv', index_col=False)[['Cycle', 'Segment', 'Orientation']]
breakpoint = pd.read_csv('breakpoint.csv', index_col=False)[['BreakpointEdge', 'Start_Coordinate', 'Start_Position_Orientation', 'End_Coordinate', 'End_Position_Orientation']]

cycle['Segment'] = cycle['Segment'].astype('Int64')
cycle['Cycle'] = cycle['Cycle'].astype('Int64')

df1 = pd.merge(cycle, segment, on='Segment', how='left') #1. Merge cycle df & segment df according to cycle segment column into a new dataframe called df1

#===================================================

# **Case 1**
# #first identify cycles only appeared once 
cycle_counts = df1['Cycle'].value_counts()
single_occurrence_cycles = cycle_counts[cycle_counts == 1].index.tolist()

df1_single = df1[df1['Cycle'].isin(single_occurrence_cycles)].copy() # Filter df1 for cycles that appear only once

#then create 2 new empty columns in df1
df1_single['5_BP'] = None
df1_single['3_BP'] = None

# For each row in df1_single, find matches in df2
#First, iterate the rows with each index value 
for idx, row in df1_single.iterrows():
    start_match = breakpoint[breakpoint['Start_Coordinate'] == row['segment_coordinate_start']]
    if not start_match.empty:
        df1_single.at[idx, '5_BP'] = start_match.iloc[0]['Start_Coordinate']

    end_match = breakpoint[breakpoint['End_Coordinate'] == row['segment_coordinate_end']]
    if not end_match.empty:
        df1_single.at[idx, '3_BP'] = end_match.iloc[0]['End_Coordinate']

df1['5_BP'] = None
df1['3_BP'] = None

for idx in df1_single.index:
    df1.at[idx, '5_BP'] = df1_single.at[idx, '5_BP']
    df1.at[idx, '3_BP'] = df1_single.at[idx, '3_BP']

#===================================================

#Case 2

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
        # For Case 2: match segment_coordinate_end with df2's Start_Coordinate
        start_match = breakpoint[breakpoint['Start_Coordinate'] == row['segment_coordinate_end']]
        if not start_match.empty:
            matched_row = start_match.iloc[0]
            # Use 5_BP for the matched Start_Coordinate
            df1.at[idx, '5_BP'] = matched_row['Start_Coordinate']
            # Use 3_BP for the corresponding End_Coordinate from same row
            df1.at[idx, '3_BP'] = matched_row['End_Coordinate']

    elif row['Orientation'] == '-':
        # For '-' orientation:
        # 5_BP = segment_coordinate_start of the FIRST segment in this cycle
        # 3_BP = segment_coordinate_start of the CURRENT negative segment
        if cycle_num in first_segment_starts:
            df1.at[idx, '5_BP'] = first_segment_starts[cycle_num]
            df1.at[idx, '3_BP'] = row['segment_coordinate_start']

print(df1)

    


                   
            
    
