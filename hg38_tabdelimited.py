import pandas as pd

# Load the raw text file (which uses tabs)
df = pd.read_csv('/Users/admin/MSc_Project/ecDNA_BP_mapping/hg38_reference', sep='\t', header=None)

# Assign the 12 standard BED columns
df.columns = [
    'chrom', 'txStart', 'txEnd', 'name', 'score', 'strand', 
    'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 
    'blockSizes', 'blockStarts'
]

# Save as a proper CSV with commas
# The 'quoting' parameter ensures internal data commas are handled correctly
df.to_csv('/Users/admin/MSc_Project/ecDNA_BP_mapping/hg38_reference.csv', sep=',', index=False)