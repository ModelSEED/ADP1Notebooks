#!/usr/bin/env python
"""Debug script to identify gene ID matching issue in MSExpression"""

import sys
# Use the local development version
sys.path.insert(0, "/home/chenry/projects/KBUtilLib/src/kbutillib/dependencies/ModelSEEDpy")
sys.path.insert(0, "/home/chenry/projects/KBUtilLib/src")

from cobra.io import load_json_model
from modelseedpy.core.msmodelutl import MSModelUtil
from modelseedpy.multiomics.msexpression import MSExpression

# Load the model
model = load_json_model("data/TranslatedPublishedModel.json")
pubmod = MSModelUtil(model)

# Load expression data
expression = MSExpression.from_spreadsheet(
    filename="data/UGA_Proteomics_May2025_Report.xlsx",
    sheet_name="Imputed",
    skiprows=1,
    type="Log2",
    id_column="Protein Accession"
)

print(f"Expression features loaded: {len(expression.features)}")
print(f"Model genes: {len(pubmod.model.genes)}")
print()

# Check first 10 model gene IDs
print("First 10 model gene IDs:")
for i, gene in enumerate(list(pubmod.model.genes)[:10]):
    print(f"  {gene.id}")
print()

# Check first 10 expression feature IDs
print("First 10 expression feature IDs:")
for i, feature in enumerate(list(expression.features)[:10]):
    print(f"  {feature.id}")
print()

# Check if model genes can be found in expression
print("Checking gene matching...")
matched = 0
not_matched = []
for gene in pubmod.model.genes:
    feature = expression.object.search_for_gene(gene.id)
    if feature is not None and feature.id in expression.features:
        matched += 1
    else:
        if len(not_matched) < 20:
            not_matched.append(gene.id)

print(f"Matched genes: {matched}")
print(f"Unmatched genes: {len(pubmod.model.genes) - matched}")
print()
if not_matched:
    print(f"Example unmatched gene IDs (first 20):")
    for gid in not_matched:
        print(f"  {gid}")

# Check if model genes have names or annotations that might match
print("\nChecking first 5 model genes for annotations:")
for gene in list(pubmod.model.genes)[:5]:
    print(f"\nGene ID: {gene.id}")
    print(f"  Name: {gene.name}")
    if hasattr(gene, 'annotation'):
        print(f"  Annotation: {gene.annotation}")

# Check expression genome features for aliases
print("\nChecking first 5 expression features for aliases:")
for feature in list(expression.object.features)[:5]:
    print(f"\nFeature ID: {feature.id}")
    if hasattr(feature, 'aliases'):
        print(f"  Aliases: {feature.aliases}")
