#!/usr/bin/env python
"""
Investigation: Check model gene ID format
Identify the gene ID naming convention in the model
"""

import sys
import os

# Change to notebooks directory
os.chdir('/home/chenry/projects/ClaudeProjects/ADP1Research/ADP1Notebooks/notebooks')

# Import by executing util.py
exec(open('util.py').read())

print("=" * 100)
print("CHECKING MODEL GENE ID FORMAT")
print("=" * 100)

# Load model
print("\nLoading model...")
model = MSModelUtil.from_cobrapy("data/TranslatedPublishedModel.json")
print(f"✓ Model loaded: {len(model.model.genes)} genes, {len(model.model.reactions)} reactions")

# Sample model gene IDs
print("\nSample of model gene IDs (first 30):")
for i, gene in enumerate(list(model.model.genes)[:30]):
    print(f"  {gene.id}")

# Check for ACIAD prefix
print("\nSearching for ACIAD-prefixed genes in model...")
aciad_genes = [g for g in model.model.genes if 'ACIAD' in g.id.upper()]
print(f"  Found {len(aciad_genes)} genes with 'ACIAD' in ID")
if aciad_genes:
    print("  Sample:")
    for gene in aciad_genes[:10]:
        print(f"    {gene.id}")

# Check DgoA reaction specifically
print("\nExamining DgoA reaction in detail...")
dgoa_rxn = model.model.reactions.get_by_id("DgoA")
print(f"  Reaction: {dgoa_rxn.id} - {dgoa_rxn.name}")
print(f"  Equation: {dgoa_rxn.reaction}")
print(f"  Genes: {[g.id for g in dgoa_rxn.genes]}")
print(f"  GPR: {dgoa_rxn.gene_reaction_rule}")

# Load proteomics and check format
print("\n" + "=" * 100)
print("COMPARING GENE ID FORMATS")
print("=" * 100)

expression = MSExpression.from_spreadsheet(
    "data/UGA_Proteomics_May2025_Report.xlsx",
    sheet_name="Imputed",
    skiprows=1
)
df = expression.get_dataframe()

print(f"\nProteomics gene IDs (first 30):")
for gene_id in df.index[:30]:
    print(f"  {gene_id}")

# Try to identify pattern
print("\n" + "=" * 100)
print("ANALYSIS")
print("=" * 100)

model_gene_ids = set([g.id for g in model.model.genes])
expression_gene_ids = set(df.index.tolist())

print(f"\nModel gene ID examples:")
for gene_id in list(model_gene_ids)[:5]:
    print(f"  {gene_id}")

print(f"\nProteomics gene ID examples:")
for gene_id in list(expression_gene_ids)[:5]:
    print(f"  {gene_id}")

print("\nHypothesis: Gene IDs need translation/mapping between formats")
print("Model uses: Simple gene names (DgoA, etc.)")
print("Proteomics uses: Locus tags (ACIAD_RSxxxxx) plus some common names")

# Check if there's a translation table in the model
print("\n" + "=" * 100)
print("CHECKING FOR GENE NAME SYNONYMS")
print("=" * 100)

# Sample a few genes and check their attributes
print("\nExamining model gene objects for synonyms/aliases...")
for gene in list(model.model.genes)[:10]:
    print(f"\nGene: {gene.id}")
    print(f"  Name: {gene.name}")
    if hasattr(gene, 'annotation'):
        print(f"  Annotation: {gene.annotation}")
    attrs = [a for a in dir(gene) if not a.startswith('_')]
    print(f"  Available attributes: {attrs}")
