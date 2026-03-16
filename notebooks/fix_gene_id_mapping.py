#!/usr/bin/env python
"""
Fix: Apply gene ID mapping to expression data
Translate ACIAD_RSxxxxx (proteomics) to ACIADxxxx (model)
"""

import sys
import os
import pandas as pd

# Change to notebooks directory
os.chdir('/home/chenry/projects/ClaudeProjects/ADP1Research/ADP1Notebooks/notebooks')

# Import by executing util.py
exec(open('util.py').read())

print("=" * 100)
print("FIXING GENE ID MAPPING")
print("=" * 100)

# Step 1: Load gene mapping
print("\nStep 1: Loading gene mapping from ADP1Genes.csv...")
gene_mapping_df = pd.read_csv("data/ADP1Genes.csv")
print(f"✓ Loaded gene mapping: {len(gene_mapping_df)} genes")
print(f"  Columns: {gene_mapping_df.columns.tolist()}")

# Create mapping dictionary: RefSeq ID -> Old ID
# Also include Name -> Old ID mapping for cases like "DgoA"
refseq_to_old = {}
name_to_old = {}

for _, row in gene_mapping_df.iterrows():
    refseq_id = row['ID']
    old_id = row['Old ID']
    name = row['Name']

    if pd.notna(refseq_id) and pd.notna(old_id):
        refseq_to_old[refseq_id] = old_id

    if pd.notna(name) and pd.notna(old_id) and name != '':
        name_to_old[name] = old_id

print(f"  RefSeq -> Old ID mappings: {len(refseq_to_old)}")
print(f"  Name -> Old ID mappings: {len(name_to_old)}")

# Check if DgoA is in the mapping
print(f"\n  DgoA mapping:")
if "DgoA" in name_to_old:
    print(f"    DgoA -> {name_to_old['DgoA']}")
else:
    print(f"    DgoA not found in name mappings")
    # Check if it's in the CSV at all
    dgoa_rows = gene_mapping_df[gene_mapping_df['Name'] == 'DgoA']
    print(f"    DgoA rows in CSV: {len(dgoa_rows)}")
    if len(dgoa_rows) > 0:
        print(dgoa_rows[['ID', 'Old ID', 'Name']])

# Step 2: Load expression data
print("\nStep 2: Loading expression data...")
expression = MSExpression.from_spreadsheet(
    "data/UGA_Proteomics_May2025_Report.xlsx",
    sheet_name="Imputed",
    skiprows=1
)
df = expression.get_dataframe()
print(f"✓ Loaded expression data: {df.shape[0]} genes, {df.shape[1]} conditions")

# Step 3: Translate gene IDs
print("\nStep 3: Translating gene IDs...")
original_gene_ids = df.index.tolist()
translated_gene_ids = []
translation_stats = {"translated": 0, "name_mapped": 0, "unmapped": 0, "kept_as_is": 0}

unmapped_genes = []

for gene_id in original_gene_ids:
    gene_id_str = str(gene_id)

    # Try RefSeq mapping first
    if gene_id_str in refseq_to_old:
        translated_gene_ids.append(refseq_to_old[gene_id_str])
        translation_stats["translated"] += 1
    # Try name mapping (for special genes like DgoA)
    elif gene_id_str in name_to_old:
        translated_gene_ids.append(name_to_old[gene_id_str])
        translation_stats["name_mapped"] += 1
    # Keep as is (might already be in correct format, or special cases)
    else:
        translated_gene_ids.append(gene_id_str)
        translation_stats["kept_as_is"] += 1
        unmapped_genes.append(gene_id_str)

print(f"  Translation statistics:")
print(f"    Translated (RefSeq -> Old): {translation_stats['translated']}")
print(f"    Name mapped: {translation_stats['name_mapped']}")
print(f"    Kept as-is: {translation_stats['kept_as_is']}")

if len(unmapped_genes) > 0:
    print(f"\n  Unmapped genes (first 20):")
    for gene in unmapped_genes[:20]:
        print(f"    {gene}")

# Step 4: Create new DataFrame with translated IDs
print("\nStep 4: Creating translated expression DataFrame...")
translated_df = df.copy()
translated_df.index = translated_gene_ids
print(f"✓ Created translated DataFrame")

# Check for duplicates
duplicate_ids = translated_df.index[translated_df.index.duplicated()].unique()
if len(duplicate_ids) > 0:
    print(f"  ⚠️  WARNING: {len(duplicate_ids)} duplicate gene IDs after translation:")
    for dup_id in duplicate_ids[:10]:
        print(f"    {dup_id}")
    print(f"  Removing duplicates by keeping first occurrence...")
    translated_df = translated_df[~translated_df.index.duplicated(keep='first')]

# Step 5: Load model and check overlap
print("\nStep 5: Checking overlap with model genes...")
model = MSModelUtil.from_cobrapy("data/TranslatedPublishedModel.json")
model_gene_ids = set([g.id for g in model.model.genes])
translated_gene_ids_set = set(translated_df.index.tolist())

overlap = model_gene_ids.intersection(translated_gene_ids_set)

print(f"  Model genes: {len(model_gene_ids)}")
print(f"  Translated expression genes: {len(translated_gene_ids_set)}")
print(f"  Overlap: {len(overlap)} genes ({100*len(overlap)/len(model_gene_ids):.1f}% of model)")

# Step 6: Check DgoA specifically
print("\nStep 6: Checking DgoA after translation...")
if "DgoA" in translated_df.index:
    print(f"  ✓ DgoA found in translated expression data")
    dgoa_values = translated_df.loc["DgoA"]
    print(f"    Sample values: {dgoa_values.iloc[:5].values}")
else:
    print(f"  ⚠️  DgoA not found in translated data")
    # Check what the DgoA rows became
    print(f"    Checking original 'DgoA' gene...")
    if "DgoA" in df.index:
        idx_position = df.index.get_loc("DgoA")
        print(f"      Original ID: DgoA")
        print(f"      Translated to: {translated_gene_ids[idx_position]}")

# Step 7: Save translated expression for testing
print("\nStep 7: Saving translated expression data...")
translated_df.to_csv("datacache/translated_expression_data.csv")
print(f"  ✓ Saved to datacache/translated_expression_data.csv")

print("\n" + "=" * 100)
print("SUMMARY")
print("=" * 100)
print(f"Original overlap: 1 gene (0.1%)")
print(f"After translation: {len(overlap)} genes ({100*len(overlap)/len(model_gene_ids):.1f}%)")
print(f"\nThis translation should be integrated into the util.py pipeline.")
