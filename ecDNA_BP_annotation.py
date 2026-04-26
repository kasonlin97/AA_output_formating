import pandas as pd

def get_detailed_annotation(chrom_pos, ref_df):
    if pd.isna(chrom_pos) or ':' not in str(chrom_pos):
        return "N/A", "N/A"
    
    # Parse Coordinate
    clean_coord = str(chrom_pos).replace('+', '').replace('-', '').strip()
    chrom_num, pos_str = clean_coord.split(':')
    chrom = f"chr{chrom_num}"
    pos = int(pos_str)
    
    # 1. Find transcripts that overlap this position
    matches = ref_df[(ref_df['chrom'] == chrom) & 
                     (ref_df['txStart'] <= pos) & 
                     (ref_df['txEnd'] >= pos)]
    
    if matches.empty:
        return "Intergenic", "N/A"

    # We take the first matching transcript (usually the longest/canonical)
    row = matches.iloc[0]
    transcript_id = row['name']
    gene_label = f"{row['geneSymbol']} ({transcript_id})" if 'geneSymbol' in row else transcript_id
    
    # 2. Parse BED12 blocks to check for Exons
    # blockStarts and blockSizes are comma-separated strings
    starts = [int(x) for x in str(row['blockStarts']).strip(',').split(',')]
    sizes = [int(x) for x in str(row['blockSizes']).strip(',').split(',')]
    
    is_exon = False
    exon_num = 0
    
    # BED blocks are 0-based relative to txStart
    for i in range(len(starts)):
        exon_start = row['txStart'] + starts[i]
        exon_end = exon_start + sizes[i]
        
        if exon_start <= pos <= exon_end:
            is_exon = True
            # For forward strand (+), exons are 1, 2, 3...
            # For reverse strand (-), exons are counted backwards
            exon_num = i + 1 if row['strand'] == '+' else len(starts) - i
            break
            
    if is_exon:
        return gene_label, f"Exon {exon_num}"
    else:
        return gene_label, "Intron"

# --- Main Implementation ---
df = pd.read_csv("/Users/admin/MSc_Project/ecDNA_BP_mapping/ecDNA_structure.csv")
ref = pd.read_csv("/Users/admin/MSc_Project/ecDNA_BP_mapping/hg38_reference.csv")

# Annotate 5' Junction
results_5 = df['5_BP_Junction'].apply(lambda x: get_detailed_annotation(x, ref))
df['5_BP_Gene'], df['5_BP_Feature'] = zip(*results_5)

# Annotate 3' Junction
results_3 = df['3_BP_Junction'].apply(lambda x: get_detailed_annotation(x, ref))
df['3_BP_Gene'], df['3_BP_Feature'] = zip(*results_3)

df.to_csv("/Users/admin/MSc_Project/ecDNA_BP_mapping/ecDNA_BP_annotation.csv", index=False)
print("Done! Detailed annotation saved.")