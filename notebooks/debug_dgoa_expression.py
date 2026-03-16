#!/usr/bin/env python
"""
Debug: Check DgoA expression values through the pipeline
"""

import sys
import os

# Change to notebooks directory
os.chdir('/home/chenry/projects/ClaudeProjects/ADP1Research/ADP1Notebooks/notebooks')

# Import by executing util.py
exec(open('util.py').read())

print("=" * 100)
print("DEBUGGING DgoA EXPRESSION THROUGH PIPELINE")
print("=" * 100)

# Step 1: Load raw expression
print("\nStep 1: Loading raw expression...")
raw_expression = MSExpression.from_spreadsheet(
    filename="data/UGA_Proteomics_May2025_Report.xlsx",
    sheet_name="Imputed",
    skiprows=1,
    type="Log2"
)
raw_df = raw_expression.get_dataframe()
print(f"✓ Loaded: {raw_df.shape}")

if "DgoA" in raw_df.index:
    print(f"\n  DgoA in raw expression:")
    print(f"    Sample values: {raw_df.loc['DgoA'].iloc[:5].values}")
else:
    print(f"\n  ⚠️  DgoA NOT in raw expression")

# Step 2: Translate
print("\nStep 2: Translating gene IDs...")
translated_expression = util.translate_expression_gene_ids(raw_expression)
translated_df = translated_expression.get_dataframe()
print(f"✓ Translated: {translated_df.shape}")

if "DgoA" in translated_df.index:
    print(f"\n  DgoA in translated expression:")
    print(f"    Sample values: {translated_df.loc['DgoA'].iloc[:5].values}")
else:
    print(f"\n  ⚠️  DgoA NOT in translated expression")

# Step 3: Average
print("\nStep 3: Averaging replicates...")
strains = ["ACN2586"]
averaged_expression = util.average_expression_replicates(translated_expression, strains)
averaged_df = averaged_expression.get_dataframe()
print(f"✓ Averaged: {averaged_df.shape}")
print(f"  Columns: {averaged_df.columns.tolist()}")

if "DgoA" in averaged_df.index:
    print(f"\n  DgoA in averaged expression:")
    print(f"    ACN2586 value: {averaged_df.loc['DgoA', 'ACN2586']}")
else:
    print(f"\n  ⚠️  DgoA NOT in averaged expression")

# Step 4: Check model
print("\nStep 4: Checking model...")
model = MSModelUtil.from_cobrapy("data/TranslatedPublishedModel.json")
print(f"✓ Model loaded: {len(model.model.genes)} genes, {len(model.model.reactions)} reactions")

if "DgoA" in [g.id for g in model.model.genes]:
    dgoa_gene = model.model.genes.get_by_id("DgoA")
    print(f"\n  DgoA gene in model:")
    print(f"    ID: {dgoa_gene.id}")
    print(f"    Name: {dgoa_gene.name}")
    print(f"    Reactions: {[r.id for r in dgoa_gene.reactions]}")
else:
    print(f"\n  ⚠️  DgoA gene NOT in model")

if "DgoA" in [r.id for r in model.model.reactions]:
    dgoa_rxn = model.model.reactions.get_by_id("DgoA")
    print(f"\n  DgoA reaction in model:")
    print(f"    ID: {dgoa_rxn.id}")
    print(f"    Name: {dgoa_rxn.name}")
    print(f"    Genes: {[g.id for g in dgoa_rxn.genes]}")
    print(f"    GPR: {dgoa_rxn.gene_reaction_rule}")
else:
    print(f"\n  ⚠️  DgoA reaction NOT in model")

# Step 5: Try to get expression value for DgoA using MSExpression
print("\nStep 5: Testing MSExpression.get_value()...")
try:
    dgoa_value = averaged_expression.get_value("DgoA", "ACN2586")
    print(f"  averaged_expression.get_value('DgoA', 'ACN2586') = {dgoa_value}")
except Exception as e:
    print(f"  ⚠️  Error: {e}")

# Step 6: Check what fit_model_flux_to_data does
print("\nStep 6: Testing fit_model_flux_to_data on a simple model...")
test_model = MSModelUtil.from_cobrapy("data/TranslatedPublishedModel.json")
print(f"  Model before fit: {len(test_model.model.reactions)} reactions")

# Add DGOA reaction
util.add_dgoa_reaction(test_model)
print(f"  Model after adding DGOA: {len(test_model.model.reactions)} reactions")

try:
    # Try to apply expression constraints
    print("\n  Calling fit_model_flux_to_data...")
    print(f"    Expression data shape: {averaged_expression._data.shape}")
    print(f"    Expression conditions: {[c.id for c in averaged_expression.conditions]}")
    print(f"    Target condition: ACN2586")

    on_on, off_off = averaged_expression.fit_model_flux_to_data(
        model=test_model.model,
        condition="ACN2586",
        default_coef=0.01
    )

    print(f"\n  ✓ fit_model_flux_to_data completed")
    print(f"    on_on reactions: {len(on_on) if on_on else 0}")
    print(f"    off_off reactions: {len(off_off) if off_off else 0}")

    if on_on:
        print(f"    Sample on_on (first 10): {list(on_on.keys())[:10]}")
    if off_off:
        print(f"    Sample off_off (first 10): {list(off_off.keys())[:10]}")

    # Check DgoA reaction expression after fit
    dgoa_rxn = test_model.model.reactions.get_by_id("DgoA")
    if hasattr(dgoa_rxn, 'gene_reaction_rule') and dgoa_rxn.gene_reaction_rule:
        print(f"\n  DgoA reaction after fit:")
        print(f"    GPR: {dgoa_rxn.gene_reaction_rule}")
        print(f"    Genes: {[g.id for g in dgoa_rxn.genes]}")

        # Try to manually calculate expression
        for gene in dgoa_rxn.genes:
            try:
                gene_expr = averaged_expression.get_value(gene.id, "ACN2586")
                print(f"    Gene {gene.id} expression: {gene_expr}")
            except:
                print(f"    Gene {gene.id} expression: NOT FOUND")

except Exception as e:
    print(f"  ❌ Error in fit_model_flux_to_data: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 100)
print("DEBUG COMPLETE")
print("=" * 100)
