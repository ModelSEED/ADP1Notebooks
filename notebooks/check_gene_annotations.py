#!/usr/bin/env python
"""Check if model genes have annotations with new locus tags"""

import sys
sys.path.insert(0, "/home/chenry/projects/KBUtilLib/src")

from cobra.io import load_json_model

print("Loading model...")
model = load_json_model("data/TranslatedPublishedModel.json")

print("\nChecking first 20 genes for notes/annotations:")
for gene in list(model.genes)[:20]:
    print(f"\nGene: {gene.id}")
    print(f"  Name: {gene.name}")
    if hasattr(gene, 'notes') and gene.notes:
        print(f"  Notes: {gene.notes}")
    if hasattr(gene, 'annotation') and gene.annotation:
        print(f"  Annotation: {gene.annotation}")

# Check if any gene has "ACIAD_RS" in its attributes
print("\n\nSearching for genes with 'ACIAD_RS' in any attribute...")
found = 0
for gene in model.genes:
    gene_str = str(gene.__dict__)
    if "ACIAD_RS" in gene_str:
        print(f"\nFound in {gene.id}:")
        print(f"  {gene.__dict__}")
        found += 1
        if found >= 5:
            break

if found == 0:
    print("No genes found with 'ACIAD_RS' in their attributes")
