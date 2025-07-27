# -*- coding: utf-8 -*-
"""
Created on Thu Jul 24 17:17:19 2025

@author: cariefrantz
"""

import pandas as pd

fpath = "C://Users/cariefrantz/Desktop/Rotten Ice Bioinformatics/18S_ASV-table_cc_ns_noUBABig_rarefaction_50000.csv"

depth_cutoff = 5263 # Desired sampling depth
representation_ratio = 0.9 # At the sampling depth, this is the fraction of max sequences considered acceptable
best_ratio_threshold = 0.1 # The ratio above which is 'undersampled' at the max depth possible for the sample

# Load and clean
table = pd.read_csv(fpath, header=0, index_col=0)
depth_cols = [col for col in table.columns if col.startswith('depth-')]
table_cleaned = table[depth_cols]

# Reverse for last value search
table_reversed = table_cleaned.iloc[:, ::-1]

# Get last value and depth
def get_last_value_and_depth(row):
    for col in row.index:
        val = row[col]
        if pd.notnull(val):
            depth = int(col.split('_')[0].split('-')[1])
            return pd.Series({'last_val': val, 'last_depth': depth})
    return pd.Series({'last_val': None, 'last_depth': None})

result = table_reversed.apply(get_last_value_and_depth, axis=1)
result['best_ratio'] = result['last_val'] / result['last_depth']

# Find closest depth ≤ depth_cutoff
depth_col_map = {
    int(col.split('_')[0].split('-')[1]): col for col in depth_cols
}
available_depths = sorted([d for d in depth_col_map if d <= depth_cutoff])

if not available_depths:
    raise ValueError(f"No depths ≤ {depth_cutoff} found.")
closest_depth = available_depths[-1]
closest_col = depth_col_map[closest_depth]
print(f"Using closest depth column ≤ {depth_cutoff}: {closest_col}")

# Get cutoff value, max per sample, and flags
value_at_cutoff = table_cleaned[closest_col]
max_per_sample = table_cleaned.max(axis=1)

result['value_at_cutoff'] = value_at_cutoff
result['max_val'] = max_per_sample
result['cutoff_ratio'] = result['value_at_cutoff'] / result['max_val']
result['cutoff_nan'] = result['value_at_cutoff'].isna()
result['undersampled'] = result['cutoff_ratio'] < representation_ratio

# Get all sample names in original order
all_samples = table_cleaned.index

# === List 1: best_ratio ≥ 0.1 ===
list1 = [s for s in all_samples if result.loc[s, 'best_ratio'] >= best_ratio_threshold]

# === List 2: not in list1 AND NaN at cutoff ===
list2 = [
    s for s in all_samples
    if s not in list1 and result.loc[s, 'cutoff_nan']
]

# === List 3: not in list1 or list2 AND undersampled ===
list3 = [
    s for s in all_samples
    if s not in list1 and s not in list2 and result.loc[s, 'undersampled']
]

def split_by_molecule(sample_list):
    cDNA = [s for s in sample_list if s.endswith('cDNA')]
    DNA = [s for s in sample_list if s.endswith('DNA')]
    return cDNA, DNA

# Split each list
list1_cDNA, list1_DNA = split_by_molecule(list1)
list2_cDNA, list2_DNA = split_by_molecule(list2)
list3_cDNA, list3_DNA = split_by_molecule(list3)

# Print results
print("\nList 1 (Insufficiently sampled):")
print("cDNA:")
print('\n'.join(list1_cDNA))
print("DNA:")
print('\n'.join(list1_DNA))

print("\nList 2 (Cutoff too deep — NaN at depth):")
print("cDNA:")
print('\n'.join(list2_cDNA))
print("DNA:")
print('\n'.join(list2_DNA))

print("\nList 3 (Not cut off, but undersampled at the depth cutoff):")
print("cDNA:")
print('\n'.join(list3_cDNA))
print("DNA:")
print('\n'.join(list3_DNA))