
import pandas as pd
import os

# 1. Define the absolute paths
segment_path = '/Users/admin/S16_AA/S16A1_segment.csv'
cycle_path = '/Users/admin/S16_AA/S16A1_cycle.csv'

# 2. Load the dataframes using the full paths
# We use index_col=False to ensure pandas doesn't treat the first column as an index
segment = pd.read_csv(segment_path, index_col=False)
cycle = pd.read_csv(cycle_path, index_col=False)[['Cycle', 'Segment', 'Orientation']]

# 3. Merge cycle and segment dataframes on the 'Segment' column
# 'how=left' ensures all rows from 'cycle' are kept
df1 = pd.merge(cycle, segment, on='Segment', how='left')

# 4. View the result
print("Merged Dataframe (df1):")
print(df1.head())

# 5. (Optional) Save the result to the same folder
output_path = '/Users/admin/S16_AA/S16A1_merged_results.csv'
df1.to_csv(output_path, index=False)