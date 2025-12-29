---
prd_number: 0001
feature_name: expression-flux-analysis
created_date: 2025-10-22
status: draft
---

# PRD: Expression-Based Flux Profile Analysis for ADP1 Strains

## Introduction/Overview

Complete the `ADP1ExpressionAnalysis.ipynb` notebook to predict metabolic flux profiles for all strains in the UGA proteomics dataset using expression-constrained flux balance analysis. The notebook will integrate proteomics data with metabolic models to predict how protein abundance levels influence metabolic flux distributions in wild-type and evolved ADP1 strains.

**Problem Statement**: Understanding how proteome-level changes in evolved ADP1 strains (with DgoA enzyme pathway) affect metabolic flux is critical for identifying metabolic adaptations. Current analysis lacks predicted flux profiles for each condition that account for measured protein abundances.

**Goal**: Generate comprehensive, expression-constrained flux predictions for 8 ADP1 strains (with and without DGOA enzyme) using ModelSEEDpy's MSExpression framework, producing both machine-readable outputs and interactive visualizations.

## Goals

1. **Complete Notebook Implementation**: Fully implement ADP1ExpressionAnalysis.ipynb with tested, working code
2. **Expression Integration**: Use MSExpression.fit_model_flux_to_data() to constrain metabolic models based on proteomics
3. **Multi-Condition Analysis**: Generate flux profiles for 8 strains averaging all replicates per strain
4. **DGOA Comparison**: Analyze each strain with and without DGOA enzyme to understand pathway dependencies
5. **Data Persistence**: Save results in multiple formats for downstream analysis
6. **Visual Validation**: Create Escher metabolic maps showing flux distributions for each condition
7. **Quality Assurance**: Validate all predictions produce biologically feasible results

## User Stories

1. **As a metabolic engineer**, I want to see predicted flux profiles for each ADP1 strain so that I can understand metabolic differences between wild-type and evolved strains.

2. **As a systems biologist**, I want to compare flux predictions with and without the DGOA enzyme so that I can quantify the enzyme's impact on central metabolism.

3. **As a data analyst**, I want flux data in both JSON and Excel formats so that I can perform custom statistical analyses and integrate with other datasets.

4. **As a researcher**, I want interactive Escher maps for each condition so that I can visually explore metabolic pathway usage and identify bottlenecks.

5. **As a developer**, I want the notebook to handle errors gracefully so that analysis failures are informative and debuggable.

## Functional Requirements

### FR1: Data Loading and Preprocessing
1.1. Load proteomics data from `data/UGA_Proteomics_May2025_Report.xlsx` (Imputed sheet, skiprows=1)
1.2. Use MSExpression.from_spreadsheet() with parameters:
   - `sheet_name="Imputed"`
   - `skiprows=1`
   - `type="Log2"`
   - `id_column="Protein Accession"`
1.3. Average all replicates for each of 8 strains: ACN2586, ACN2821, ACN3015, ACN3468, ACN3471, ACN3474, ACN3477, ADP1
1.4. Load TranslatedPublishedModel.json using MSModelUtil.from_cobrapy()
1.5. Load pyruvate media: `KBaseMedia/Carbon-Pyruvic-Acid`

### FR2: Model Configuration
2.1. Create two model variants for each strain:
   - Variant A: With DGOA enzyme added (util.add_dgoa_reaction())
   - Variant B: Without DGOA enzyme (native model)
2.2. For all models, knock out native DAHP synthase (rxn01332_c0):
   - Set lower_bound = 0
   - Set upper_bound = 0
2.3. Set media to pyruvate for all models
2.4. Constrain biomass objective (GROWTH_DASH_RXN) to minimum 25% of optimal using util.constrain_objective_to_fraction_of_optimum()

### FR3: Expression-Constrained Flux Prediction
3.1. For each of 8 strains:
   - Create condition-specific model copies (16 total: 8 strains × 2 DGOA variants)
   - Apply MSExpression.fit_model_flux_to_data() with parameters:
     - `model`: condition-specific model
     - `condition`: strain name (e.g., "ACN2586")
     - `default_coef=0.01`
     - `activation_threshold=None` (disable activation constraints)
     - `deactivation_threshold=0.000001` (reactions with very low expression)
3.2. Extract from fit_model_flux_to_data() output:
   - Flux solution (pFBA result)
   - Reactions turned off (off_off list)
   - Objective value (biomass production)
   - Gene-reaction activity mapping
3.3. Apply deactivation constraints:
   - For reactions in off_off list, set bounds to 0 (both lower and upper)
3.4. Run final pFBA to get expression-constrained flux distribution

### FR4: Data Output - JSON Format
4.1. Save individual JSON files for each condition to `datacache/`:
   - `{strain}_with_dgoa_fluxes.json` - Flux dictionary (16 files for with DGOA)
   - `{strain}_no_dgoa_fluxes.json` - Flux dictionary (16 files for without DGOA)
4.2. JSON structure for flux files:
   ```json
   {
     "rxn_id": flux_value,
     "rxn01234_c0": 1.234,
     ...
   }
   ```
4.3. Save comprehensive analysis summary:
   - `expression_flux_analysis_summary.json` with structure:
   ```json
   {
     "strain_name": {
       "with_dgoa": {
         "biomass": 0.XX,
         "active_reactions": 123,
         "off_reactions": ["rxn_id", ...],
         "dgoa_flux": 0.XX
       },
       "without_dgoa": {
         "biomass": 0.XX,
         "active_reactions": 123,
         "off_reactions": ["rxn_id", ...]
       }
     }
   }
   ```

### FR5: Data Output - Excel Format
5.1. Create `nboutput/expression_flux_analysis.xlsx` with sheets:
   - Sheet "Summary": Overview table with columns:
     - Strain
     - DGOA_Status (with/without)
     - Biomass_Flux
     - Active_Reactions_Count
     - Off_Reactions_Count
     - DGOA_Flux (if applicable)
   - Sheet per strain (16 sheets): "{Strain}_with_DGOA" / "{Strain}_no_DGOA"
     - Columns: Reaction_ID, Reaction_Name, Flux, Gene_Association
     - Include only non-zero fluxes (|flux| > 1e-9)
5.2. Use pandas ExcelWriter to create multi-sheet workbook

### FR6: Visualization - Escher Maps
6.1. Use util.create_map_html() from EscherUtils (KBUtilLib) for each condition
6.2. Parameters for create_map_html():
   - `model`: condition-specific model
   - `flux_solution`: flux dictionary from pFBA
   - `map_name="Core"` (central metabolism map)
   - `output_file="nboutput/escher_{strain}_{dgoa_status}.html"`
   - `title="{Strain} - {DGOA Status} - Pyruvate Media"`
6.3. Generate 16 HTML files total (8 strains × 2 DGOA variants)
6.4. Maps should visualize:
   - Reaction fluxes (color intensity and arrow width)
   - Active vs inactive reactions
   - DGOA flux if present

### FR7: Validation and Quality Checks
7.1. For each flux solution, validate:
   - Biomass production > 0 (must be viable)
   - Solution status is "optimal"
   - Number of active reactions is reasonable (>50, <500)
   - DGOA flux > 0 when DGOA present and rxn01332 knocked out
7.2. Generate validation report printed to notebook output:
   ```
   Validation Summary:
   -------------------
   ACN2586_with_DGOA: ✓ Biomass=0.XX, Active=123, DGOA=0.XX
   ACN2586_no_DGOA: ✓ Biomass=0.XX, Active=120
   ...
   ```
7.3. If any validation fails, print ERROR message with details
7.4. Create correlation analysis between protein abundance and flux for major pathways

### FR8: Error Handling
8.1. Wrap all major operations in try-except blocks
8.2. Log errors with context (which strain, which step)
8.3. Continue processing other conditions if one fails
8.4. At end, report summary of successes and failures
8.5. Save error log to `nboutput/analysis_errors.log`

### FR9: Progress Tracking
9.1. Use tqdm or print statements to show progress:
   - "Processing strain 1/8: ACN2586..."
   - "Running with DGOA variant..."
   - "Running without DGOA variant..."
   - "Generating Escher map..."
9.2. Display execution time for each major step

### FR10: Code Organization
10.1. Structure notebook with clear markdown sections:
   - Introduction and setup
   - Data loading and preprocessing
   - Model configuration
   - Expression integration (main loop)
   - Validation and statistics
   - Output generation
   - Summary and conclusions
10.2. Use helper functions for repeated operations:
   - `process_strain_condition(strain, expression, model, media, with_dgoa)`
   - `validate_flux_solution(solution, strain, dgoa_status)`
   - `save_flux_outputs(strain, dgoa_status, fluxes, analysis_data)`
10.3. Include docstrings for all helper functions
10.4. Add inline comments explaining complex operations

## Non-Goals (Out of Scope)

1. **Statistical Comparison**: Not performing statistical tests between conditions (future work)
2. **Parameter Optimization**: Not optimizing MSExpression thresholds (using standard values)
3. **Alternative Media**: Only analyzing pyruvate media (other media conditions are future work)
4. **Mutant Growth Data**: Not integrating MutantGrowthRatesData.xls (different analysis)
5. **Custom Map Design**: Not creating custom Escher maps (using existing Core map)
6. **Web Interface**: Not creating interactive dashboard (static HTML outputs only)
7. **Pathway Enrichment**: Not performing pathway enrichment analysis
8. **Time-Series Analysis**: Not analyzing temporal dynamics (single time point)

## Design Considerations

### Notebook Structure
```python
# Section 1: Setup and Imports
%run util.py
import pandas as pd
from tqdm import tqdm

# Section 2: Load Data
proteomics_data = MSExpression.from_spreadsheet(...)
model = MSModelUtil.from_cobrapy("data/TranslatedPublishedModel.json")
media = util.get_media("KBaseMedia/Carbon-Pyruvic-Acid")

# Section 3: Define Strains
strains = ["ACN2586", "ACN2821", "ACN3015", "ACN3468",
           "ACN3471", "ACN3474", "ACN3477", "ADP1"]

# Section 4: Main Analysis Loop
results = {}
for strain in tqdm(strains):
    for dgoa_status in ["with_dgoa", "without_dgoa"]:
        # Process condition
        results[f"{strain}_{dgoa_status}"] = process_strain_condition(...)

# Section 5: Validation
validate_all_results(results)

# Section 6: Save Outputs
save_json_outputs(results)
save_excel_output(results)
generate_escher_maps(results)

# Section 7: Summary Statistics
display_summary_table(results)
```

### Helper Function Signatures
```python
def process_strain_condition(
    strain: str,
    expression: MSExpression,
    base_model: MSModelUtil,
    media: MSMedia,
    with_dgoa: bool
) -> dict:
    """Process single strain condition with expression data.

    Returns:
        dict with keys: fluxes, biomass, active_rxns, off_rxns, dgoa_flux
    """
    pass

def validate_flux_solution(
    solution: cobra.Solution,
    strain: str,
    dgoa_status: str,
    dgoa_flux: float = None
) -> bool:
    """Validate that flux solution is biologically feasible."""
    pass

def save_flux_outputs(
    strain: str,
    dgoa_status: str,
    fluxes: dict,
    analysis_data: dict
) -> None:
    """Save flux data to JSON and add to summary."""
    pass
```

### EscherUtils Integration
Based on KBUTILLIB_CODEBASE_RESEARCH.md (lines 440-462), use:
```python
util.create_map_html(
    model=condition_model.model,  # COBRApy model object
    flux_solution=flux_dict,       # {rxn_id: flux_value}
    map_name="Core",               # Central metabolism map
    output_file=f"nboutput/escher_{strain}_{dgoa_status}.html",
    title=f"{strain} - {dgoa_status} - Pyruvate"
)
```

Note: `create_map_html` is from EscherUtils class (escher_utils.py:29-303)

## Technical Considerations

### Dependencies
- **Required**: modelseedpy (MSExpression, MSModelUtil), cobrakbase, cobra, pandas, openpyxl
- **KBUtilLib**: Must have EscherUtils available (via util.py inheritance)
- **Data Files**: TranslatedPublishedModel.json, UGA_Proteomics_May2025_Report.xlsx must exist

### Model Constraints
- **Biomass Constraint**: Minimum 25% of optimal ensures viable but allows metabolic flexibility
- **Expression Thresholds**:
  - activation_threshold=None (disable to avoid over-constraining)
  - deactivation_threshold=0.000001 (very low expression → off)
- **DGOA Bounds**: When added, allow 0-1000 mmol/gDW/hr (util.add_dgoa_reaction default)

### Performance
- **Expected Runtime**: ~5-10 minutes for 16 conditions (depends on model size ~874 reactions)
- **Memory**: Expect ~2-4 GB RAM (16 model copies + flux solutions)
- **Parallel Processing**: Not implemented (sequential for clarity, can parallelize in future)

### Data Integrity
- **Proteomics Validation**: Check that all 8 strains have data in spreadsheet
- **Gene-Reaction Mapping**: MSExpression handles automatic mapping via model gene IDs
- **Reaction Bounds**: Verify deactivated reactions don't break model feasibility

## Success Metrics

### Functional Success
1. ✅ Notebook runs end-to-end without errors
2. ✅ All 16 conditions produce optimal flux solutions
3. ✅ All biomass fluxes > 0.01 mmol/gDW/hr
4. ✅ DGOA flux > 0 for all "with_dgoa" conditions
5. ✅ All output files created (32 JSON + 1 Excel + 16 HTML = 49 files)

### Validation Success
6. ✅ Correlation between protein abundance and flux > 0.3 for well-covered pathways
7. ✅ Number of active reactions similar across conditions (within 20% variance)
8. ✅ Escher maps visually show expected central metabolism activity
9. ✅ DGOA-dependent strains show different flux through shikimate pathway vs wild-type

### Code Quality
10. ✅ All helper functions have docstrings
11. ✅ No hardcoded file paths (use relative paths)
12. ✅ Clear progress indicators during execution
13. ✅ Validation report clearly indicates pass/fail for each condition

## Testing Requirements

### Unit Testing (within notebook)
1. **Test Data Loading**: Verify expression object has 8 conditions
2. **Test Model Setup**: Verify rxn01332 is knocked out, DGOA is added correctly
3. **Test Expression Integration**: Verify fit_model_flux_to_data returns expected keys
4. **Test Flux Extraction**: Verify flux dictionaries have correct structure
5. **Test File Output**: Verify JSON files are valid and readable

### Integration Testing
6. **End-to-End Test**: Run full notebook on all 8 strains
7. **Comparison Test**: Verify with_dgoa vs without_dgoa shows expected differences
8. **Visualization Test**: Open generated Escher maps in browser, verify they load
9. **Excel Test**: Open Excel file, verify all sheets present and formatted correctly

### Statistical Validation
10. **Biomass Validation**: All biomass values in reasonable range (0.01-0.5)
11. **Flux Balance**: For each solution, verify mass balance is maintained
12. **Expression Correlation**: Calculate Pearson correlation between abundance and flux
13. **Outlier Detection**: Identify any conditions with anomalous results

### Edge Case Testing
14. **Missing Genes**: Verify graceful handling if gene not in model
15. **Zero Expression**: Verify reactions with zero expression are handled correctly
16. **Infeasible Constraints**: Verify error handling if expression makes model infeasible

## Open Questions

1. **Replication Averaging**: Should we average replicates before creating MSExpression object, or does MSExpression have built-in averaging? Need to verify MSExpression API.

2. **DGOA Model Feasibility**: For wild-type ADP1, is the model feasible with DGOA added but rxn01332 knocked out? May need to only knock out rxn01332 for non-wild-type strains.

3. **Expression Mapping**: Do all genes in the proteomics data map to genes in the model? Should we report unmapped genes?

4. **Threshold Optimization**: Is deactivation_threshold=0.000001 appropriate for log2-transformed proteomics data? May need empirical testing.

5. **Map Availability**: Is the "Core" Escher map available in the KBUtilLib data? May need to verify map exists or use different map name.

6. **Strain Metadata**: Should we include strain descriptions (from strain_data in research doc) in the outputs?

7. **Comparison Baseline**: Should we include a "no expression constraint" baseline for comparison?

8. **FVA Integration**: Should we run FVA on expression-constrained models to get flux ranges, or just pFBA point solutions?

---

## Acceptance Criteria

### Must Have
- [x] Notebook loads proteomics data for 8 strains
- [x] Creates 16 model variants (8 strains × 2 DGOA conditions)
- [x] Applies MSExpression.fit_model_flux_to_data() to all conditions
- [x] Saves 32 JSON flux files
- [x] Saves 1 Excel file with summary + per-condition sheets
- [x] Generates 16 Escher HTML maps
- [x] Validates all solutions have positive biomass
- [x] Prints clear validation summary
- [x] Code is well-commented and organized

### Should Have
- [x] Helper functions for repeated operations
- [x] Progress indicators (tqdm or print statements)
- [x] Error handling with informative messages
- [x] Correlation analysis between expression and flux
- [x] Execution time reporting

### Nice to Have
- [ ] Comparison plots (scatter: expression vs flux)
- [ ] Heatmap of flux differences across conditions
- [ ] Statistical summary table with mean/std of fluxes
- [ ] Export reusable functions to util.py for future use

---

**Document Status**: Ready for Implementation
**Priority**: High
**Estimated Effort**: 4-6 hours development + 2-3 hours testing
**Dependencies**: ADP1_Project_Research.md, KBUTILLIB_CODEBASE_RESEARCH.md, existing util.py module
