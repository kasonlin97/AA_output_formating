
import pandas as pd
import os 

file_path = os.path.join("/Users/kason/MSc_ACS_Thesis/S16_amplicon/S16A1", "sample_16_amplicon1_annotated_cycles (1).txt")

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

df_segment['Segment'] = df_segment['Segment'].astype(str) + " " + df_segment['Segment_number'].astype(str)
df_segment = df_segment.drop(columns=['Segment_number'])
df_segment['Segment'] = df_segment['Segment'].str.replace('Segment', '')

df_segment['segment_coordinate_start'] = df_segment['Segment_chromosome'].astype(str) + ":" + df_segment['segment_coordinate_start'].astype(str)
df_segment['segment_coordinate_end'] = df_segment['Segment_chromosome'].astype(str) + ":" + df_segment['segment_coordinate_end'].astype(str)
df_segment = df_segment.drop(columns=['Segment_chromosome'])

print(df_segment)
df_segment.to_csv('/Users/kason/MSc_ACS_Thesis/S16_amplicon/S16A1/S16A1_segment.csv', index=False)