#!/usr/bin/env python
"""
Investigation: Expression Data Integration
Check if reactions are getting nonzero expression scores in MSExpression
Specifically investigate DgoA reaction expression
"""

import sys
import os
import pandas as pd

# Change to notebooks directory
os.chdir('/home/chenry/projects/ClaudeProjects/ADP1Research/ADP1Notebooks/notebooks')

# Import by executing util.py
exec(open('util.py').read())

print("=" * 100)
print("INVESTIGATION: Expression Data Integration")
print("=" * 100)

# Step 1: Load raw proteomics data
print("\nStep 1: Loading raw proteomics data...")
expression = MSExpression.from_spreadsheet(
    "data/UGA_Proteomics_May2025_Report.xlsx",
    sheet_name="Imputed",
    skiprows=1
)
print(f"✓ Loaded expression data")
print(f"  MSExpression type: {type(expression)}")
print(f"  MSExpression attributes: {[attr for attr in dir(expression) if not attr.startswith('_')]}")

# Explore the structure
if hasattr(expression, 'genes'):
    print(f"  Number of genes: {len(expression.genes)}")
if hasattr(expression, 'conditions'):
    print(f"  Number of conditions: {len(expression.conditions)}")
    print(f"  Condition names: {[c.id for c in expression.conditions][:5]}... (showing first 5)")
if hasattr(expression, 'gene_ids'):
    print(f"  Number of gene_ids: {len(expression.gene_ids)}")
if hasattr(expression, 'condition_list'):
    print(f"  Number of conditions in condition_list: {len(expression.condition_list)}")
    print(f"  Condition names: {expression.condition_list[:5]}... (showing first 5)")

# Step 2: Examine the raw expression DataFrame
print("\nStep 2: Examining raw expression DataFrame...")
df = expression.get_dataframe()
print(f"  DataFrame shape: {df.shape}")
print(f"  Column names: {list(df.columns)[:10]}... (showing first 10)")
print(f"\nFirst few rows of expression data:")
print(df.head(10))

# Step 3: Look for DgoA specifically
print("\nStep 3: Searching for DgoA gene in expression data...")
dgoa_gene_patterns = ['dgoa', 'dgoA', 'DgoA', 'DGOA']
found_dgoa = []

# Get gene IDs from the DataFrame index
gene_ids = df.index.tolist()
print(f"  Total genes in expression data: {len(gene_ids)}")

for gene_id in gene_ids:
    gene_id_lower = str(gene_id).lower()
    if any(pattern.lower() in gene_id_lower for pattern in dgoa_gene_patterns):
        found_dgoa.append(gene_id)
        print(f"  Found gene: {gene_id}")

if found_dgoa:
    print(f"\n✓ Found {len(found_dgoa)} DgoA-related genes")
    for gene_id in found_dgoa:
        print(f"\n  Expression values for {gene_id}:")
        gene_data = df.loc[gene_id]
        print(gene_data)
else:
    print("  ⚠️  No DgoA genes found in expression data")
    print("\n  Sample of gene IDs in expression data:")
    for gene_id in gene_ids[:20]:
        print(f"    {gene_id}")

# Step 4: Load model and check DgoA reaction
print("\nStep 4: Loading model to identify DgoA reaction...")
model = MSModelUtil.from_cobrapy("data/TranslatedPublishedModel.json")
print(f"✓ Model loaded: {len(model.model.reactions)} reactions")

# Search for DgoA reaction
print("\nSearching for DgoA-related reactions...")
dgoa_reactions = []
for rxn in model.model.reactions:
    rxn_id = rxn.id.lower()
    rxn_name = rxn.name.lower() if rxn.name else ""
    if 'dgoa' in rxn_id or 'dgoa' in rxn_name:
        dgoa_reactions.append(rxn)
        print(f"  Found: {rxn.id} - {rxn.name}")
        print(f"    Genes: {[g.id for g in rxn.genes]}")

if not dgoa_reactions:
    print("  ⚠️  No reactions with 'dgoa' in ID or name")
    print("\n  Let me search for aromatic amino acid biosynthesis reactions...")
    for rxn in model.model.reactions:
        if rxn.name and ('aromatic' in rxn.name.lower() or 'dahp' in rxn.name.lower()):
            print(f"    {rxn.id} - {rxn.name}")
            print(f"      Genes: {[g.id for g in rxn.genes]}")

# Step 5: Check reaction-gene associations
print("\nStep 5: Examining gene-reaction associations...")
print(f"Total genes in model: {len(model.model.genes)}")

# Look for genes that match our proteomics data
print("\nChecking if proteomics gene IDs match model gene IDs...")
expression_gene_ids = set(df.index.tolist())
model_gene_ids = set([g.id for g in model.model.genes])

overlap = expression_gene_ids.intersection(model_gene_ids)
print(f"  Expression data genes: {len(expression_gene_ids)}")
print(f"  Model genes: {len(model_gene_ids)}")
print(f"  Overlapping genes: {len(overlap)}")
print(f"  Overlap percentage: {100*len(overlap)/len(model_gene_ids):.1f}% of model genes")

if len(overlap) > 0:
    print(f"\n  Sample overlapping genes:")
    for gene_id in list(overlap)[:10]:
        print(f"    {gene_id}")
else:
    print("\n  ⚠️  WARNING: NO OVERLAP between proteomics and model gene IDs!")
    print("\n  Sample expression gene IDs:")
    for gene_id in list(expression_gene_ids)[:10]:
        print(f"    {gene_id}")
    print("\n  Sample model gene IDs:")
    for gene_id in list(model_gene_ids)[:10]:
        print(f"    {gene_id}")

# Step 6: Check MSExpression reaction scores
print("\nStep 6: Examining MSExpression reaction scores...")
print("Attempting to get expression values for reactions...")

# Average replicates first
STRAINS = ["ACN2586", "ACN2821", "ACN3015", "ACN3468", "ACN3471", "ACN3474", "ACN3477", "ADP1"]
averaged_expression = util.average_expression_replicates(expression, STRAINS)
print(f"✓ Averaged expression data to {len(averaged_expression.conditions)} conditions")

# Try to compute reaction expression for one strain
print("\nComputing reaction expression scores for ACN2586...")
test_condition = averaged_expression.conditions.get_by_id("ACN2586")

# Check if expression has a method to calculate reaction scores
print(f"\nMSExpression object methods:")
methods = [m for m in dir(averaged_expression) if not m.startswith('_')]
print(f"  Available methods: {methods[:20]}... (showing first 20)")

# Try to get reaction expression values
if hasattr(averaged_expression, 'get_reaction_expression'):
    print("\n✓ get_reaction_expression method exists")
    # Try with a sample reaction
    sample_rxn = model.model.reactions[0]
    try:
        rxn_expr = averaged_expression.get_reaction_expression(sample_rxn.id, "ACN2586")
        print(f"  Sample reaction {sample_rxn.id} expression: {rxn_expr}")
    except Exception as e:
        print(f"  Error getting reaction expression: {e}")

# Step 7: Check the actual expression DataFrame for nonzero values
print("\nStep 7: Checking for nonzero expression values...")
averaged_df = averaged_expression.get_dataframe()
print(f"\nExpression data statistics:")
print(averaged_df.describe())

print(f"\nChecking if ANY genes have finite (non -inf) expression values...")
finite_counts = {}
for col in averaged_df.columns:
    finite_values = averaged_df[col][averaged_df[col] > -float('inf')]
    finite_counts[col] = len(finite_values)
    print(f"  {col}: {finite_counts[col]} genes with finite expression")

# Check actual values
print(f"\nSample of expression values for ACN2586:")
if "ACN2586" in averaged_df.columns:
    acn2586_data = averaged_df["ACN2586"]
    print(f"  Min: {acn2586_data.min()}")
    print(f"  Max: {acn2586_data.max()}")
    print(f"  Mean: {acn2586_data.mean()}")
    print(f"  Non -inf values: {sum(acn2586_data > -float('inf'))}")

    # Show some actual values
    print(f"\n  Top 10 highest expression values:")
    top_genes = acn2586_data.nlargest(10)
    for gene_id, value in top_genes.items():
        print(f"    {gene_id}: {value}")

print("\n" + "=" * 100)
print("Investigation complete. Summary:")
print("=" * 100)
print("1. Check gene ID overlap between proteomics and model")
print("2. Verify DgoA gene is present in proteomics data")
print("3. Confirm expression values are being loaded correctly")
print("4. Identify if issue is with gene naming or data integration")
