

import pandas as pd
import os 

file_path = os.path.join("/Users/kason/MSc_ACS_Thesis/S16_amplicon/S16A2", "sample_16_amplicon2_annotated_cycles.txt")

# First, find where the cycle data starts
with open(file_path, 'r') as f:
    lines = f.readlines()

# Find the first line that starts with 'Cycle='
skiprows = 0
for i, line in enumerate(lines):
    if line.startswith('Cycle='):
        skiprows = i
        break

column_headers = [
    'Cycle', 
    'Copy_count', 
    'Length', 
    'IsCyclicPath', 
    'CycleClass', 
    'Segments', 
]

df_cycle = pd.read_csv(
    file_path, 
    sep=';', 
    skiprows=skiprows,  
    names=column_headers
)


for col in df_cycle.columns:
   df_cycle[col] = df_cycle[col].str.split('=').str[-1]

df_cycle['Segments'] = df_cycle['Segments'].str.split(',')
df_cycle_exploded = df_cycle.explode('Segments')

# 1. Create the filtered dataframe (make sure to use .copy())
filtered_df = df_cycle_exploded[df_cycle_exploded['CycleClass'] == 'ecDNA-like'].copy()

# 2. Use the SAME name 'filtered_df' to perform the extraction; Remove the original 'Segments' column
filtered_df[['Segment', 'Orientation']] = filtered_df['Segments'].str.extract(r'(\d+)([+-])')
filtered_df = filtered_df.drop(columns=['Segments'])

filtered_df['Segment'] = pd.to_numeric(filtered_df['Segment'], downcast='integer', errors='coerce')
filtered_df['Cycle'] = pd.to_numeric(filtered_df['Cycle'], downcast='integer', errors='coerce')
#print(filtered_df['Segment'])

# 3. Use 'filtered_df' for printing and saving
print(filtered_df)
filtered_df.to_csv('/Users/kason/MSc_ACS_Thesis/S16_amplicon/S16A1/S16A2_cycle.csv', index=False)




