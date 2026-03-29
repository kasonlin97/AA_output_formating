
import pandas as pd
import os 

file_path = os.path.join("/Users/admin/S16_AA", "sample_16_amplicon2_annotated_cycles.txt")

column_headers = [
    'Segment', 
    'Segment_number', 
    'Segment_chromosome', 
    'segment_coordinate_start', 
    'segment_coordinate_end', 
]

df_segment = pd.read_csv(
    file_path, 
    sep='\t',  
    names=column_headers, 
    comment = 'L'
)

df_segment = df_segment[df_segment['Segment'] == 'Segment']

# 1. Convert numeric columns to integers first to remove the .0
df_segment['Segment_number'] = df_segment['Segment_number'].astype(int)
df_segment['segment_coordinate_start'] = df_segment['segment_coordinate_start'].astype(int)
df_segment['segment_coordinate_end'] = df_segment['segment_coordinate_end'].astype(int)

# 2. Proceed with string formatting
df_segment['Segment'] = df_segment['Segment_number'].astype(str)
# No need to drop or replace 'Segment' if you just assign the number directly

df_segment['segment_coordinate_start'] = df_segment['Segment_chromosome'].astype(str) + ":" + df_segment['segment_coordinate_start'].astype(str)
df_segment['segment_coordinate_end'] = df_segment['Segment_chromosome'].astype(str) + ":" + df_segment['segment_coordinate_end'].astype(str)

df_segment = df_segment[['Segment', 'segment_coordinate_start', 'segment_coordinate_end']]

print(df_segment)
df_segment.to_csv('/Users/admin/S16_AA/S16A2_segment.csv', index=False)