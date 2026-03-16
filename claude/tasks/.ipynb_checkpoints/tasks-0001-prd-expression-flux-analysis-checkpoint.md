# Tasks: Expression-Based Flux Profile Analysis for ADP1 Strains

Based on PRD: `0001-prd-expression-flux-analysis.md`

## Relevant Files

- `ADP1Notebooks/notebooks/util.py` - Core utility module extending KBUtilLib (ALL imports and helper functions go here)
- `ADP1Notebooks/notebooks/ADP1ExpressionAnalysis.ipynb` - Main notebook (should only contain `%run util.py` and calls to util methods)
- `ADP1Notebooks/notebooks/data/TranslatedPublishedModel.json` - Published ADP1 metabolic model with ModelSEED compound IDs
- `ADP1Notebooks/notebooks/data/UGA_Proteomics_May2025_Report.xlsx` - Proteomics data with Imputed sheet containing strain replicates
- `ADP1Notebooks/notebooks/datacache/*.json` - Cached intermediate results (auto-saved via util.save())
- `ADP1Notebooks/notebooks/nboutput/*.html` - Escher map visualizations
- `ADP1Notebooks/notebooks/nboutput/expression_flux_analysis.xlsx` - Excel output with summary and flux tables

### Notes

- **Design Pattern**: util.py contains ALL imports and helper methods; notebooks call `%run util.py` then use `util.method()`
- util.py already imports: MSExpression, MSModelUtil, pfba, pandas, cobra, etc. (lines 20-27)
- util.py extends: MSFBAUtils, AICurationUtils, NotebookUtils, KBGenomeUtils, EscherUtils (line 30)
- **Singleton pattern**: `util = NotebookUtil()` instance created at bottom of util.py (line 163)
- **Caching pattern**: Save intermediate results with `util.save("name", data)` - enables resuming from any step
- **Inherited methods available**: get_model(), get_media(), run_fba(), run_fva(), constrain_objective_to_fraction_of_optimum(), save(), load(), display_dataframe(), create_map_html()
- **New methods needed in util.py**: Expression analysis workflow, strain processing, validation, Excel export

## Tasks

- [ ] 1.0 Add Helper Methods to util.py
  - [x] 1.1 Add method `process_strain_with_expression(self, strain, expression_data, base_model, media, with_dgoa)` - Processes single strain condition with expression constraints, returns dict with fluxes, biomass, active_reactions, off_reactions, dgoa_flux, solution_status
  - [x] 1.2 Add method `average_expression_replicates(self, expression, strain_list)` - Averages replicates per strain (e.g., ACN2586_1 through ACN2586_5 → ACN2586), returns new MSExpression object
  - [x] 1.3 Add method `validate_expression_flux_solution(self, solution_dict, strain, dgoa_status)` - Validates biomass >0, status=="optimal", active_reactions >50, dgoa_flux >0 if with_dgoa, returns (pass/fail, message)
  - [x] 1.4 Add method `create_expression_flux_summary(self, results_dict)` - Creates summary dict with structure from PRD FR4.3, returns summary dict
  - [x] 1.5 Add method `export_expression_flux_to_excel(self, results_dict, output_file)` - Creates multi-sheet Excel workbook with summary + per-condition sheets, saves to nboutput/
  - [ ] 1.6 Add method `generate_all_escher_maps(self, results_dict, models_dict)` - Loops through all conditions and generates Escher maps using create_map_html(), handles errors gracefully
  - [ ] 1.7 Add method `run_expression_flux_analysis(self, strains, proteomics_file, model_file, media_id)` - Main orchestrator method that runs entire pipeline, caches at each step, returns results_dict

- [ ] 2.0 Implement util.py Helper Methods
  - [ ] 2.1 Implement `process_strain_with_expression()`: Create model copy, apply DGOA if needed, knock out rxn01332, set media, constrain biomass 25%, call fit_model_flux_to_data(), apply deactivation constraints, run pFBA, extract fluxes/stats
  - [ ] 2.2 Implement `average_expression_replicates()`: Group expression columns by strain prefix, average replicates, return new MSExpression with averaged data
  - [ ] 2.3 Implement `validate_expression_flux_solution()`: Check biomass >0, status=="optimal", active_reactions >50, dgoa_flux >0 if with_dgoa, return (True/False, validation_message)
  - [ ] 2.4 Implement `create_expression_flux_summary()`: Loop through results_dict, extract key metrics (biomass, active_reactions, off_reactions, dgoa_flux), structure per PRD FR4.3
  - [ ] 2.5 Implement `export_expression_flux_to_excel()`: Create summary DataFrame, create per-condition DataFrames with reaction details, use pd.ExcelWriter to create multi-sheet workbook
  - [ ] 2.6 Implement `generate_all_escher_maps()`: Loop through results_dict, call self.create_map_html() for each condition, wrap in try-except, save to nboutput/, return list of created files
  - [ ] 2.7 Implement `run_expression_flux_analysis()`: Orchestrator that loads data, averages replicates, processes all strains × DGOA variants, caches intermediate results, returns complete results_dict

- [ ] 3.0 Create Simple Notebook Implementation
  - [ ] 3.1 Add markdown cell: "# ADP1 Expression-Constrained Flux Analysis" with description and objectives
  - [ ] 3.2 Add cell: `%run util.py` (this is the ONLY import needed)
  - [ ] 3.3 Add markdown cell: "## Define Analysis Parameters"
  - [ ] 3.4 Add cell defining: `STRAINS = ["ACN2586", "ACN2821", "ACN3015", "ACN3468", "ACN3471", "ACN3474", "ACN3477", "ADP1"]`
  - [ ] 3.5 Add markdown cell: "## Run Expression Flux Analysis Pipeline"
  - [ ] 3.6 Add cell calling: `results = util.run_expression_flux_analysis(STRAINS, "data/UGA_Proteomics_May2025_Report.xlsx", "data/TranslatedPublishedModel.json", "KBaseMedia/Carbon-Pyruvic-Acid")`
  - [ ] 3.7 Add markdown cell: "## Validate Results"
  - [ ] 3.8 Add cell: Loop through results, call `util.validate_expression_flux_solution()`, print validation summary table
  - [ ] 3.9 Add markdown cell: "## Generate Outputs"
  - [ ] 3.10 Add cell: Save individual flux files with `util.save(f"{strain}_{dgoa}_fluxes", results[key]["fluxes"])` for each condition
  - [ ] 3.11 Add cell: `summary = util.create_expression_flux_summary(results); util.save("expression_flux_analysis_summary", summary)`
  - [ ] 3.12 Add cell: `util.export_expression_flux_to_excel(results, "nboutput/expression_flux_analysis.xlsx")`
  - [ ] 3.13 Add cell: `escher_files = util.generate_all_escher_maps(results, models_dict); print(f"Generated {len(escher_files)} Escher maps")`
  - [ ] 3.14 Add markdown cell: "## Summary Statistics" with analysis of results

- [ ] 4.0 Test and Debug util.py Methods
  - [ ] 4.1 Test `average_expression_replicates()`: Load raw proteomics, average, verify 8 conditions result
  - [ ] 4.2 Test `process_strain_with_expression()`: Run on one strain (ACN2586), both DGOA variants, verify output structure
  - [ ] 4.3 Test validation method: Pass valid and invalid solutions, verify correct pass/fail results
  - [ ] 4.4 Test summary creation: Pass mock results_dict, verify summary structure matches PRD
  - [ ] 4.5 Test Excel export: Create small results dict, export, verify file creation and sheet structure
  - [ ] 4.6 Test Escher generation: Generate 2 maps, verify HTML files created and valid
  - [ ] 4.7 Test full orchestrator with 2 strains: Verify caching, intermediate outputs, final results

- [ ] 5.0 Run Full Pipeline and Validate
  - [ ] 5.1 Clear all outputs, restart kernel
  - [ ] 5.2 Run all notebook cells sequentially
  - [ ] 5.3 Monitor for errors, check intermediate caching is working
  - [ ] 5.4 Verify datacache/ has: 32 flux JSON files, 1 summary JSON, expression intermediate caches
  - [ ] 5.5 Verify nboutput/ has: 1 Excel file (expression_flux_analysis.xlsx), 16 Escher HTML files
  - [ ] 5.6 Open Excel: Verify Summary sheet with 16 rows, verify 16 condition sheets with flux data
  - [ ] 5.7 Validation checks: All biomass >0.01, all DGOA flux >0 for with_dgoa, all status=="optimal"
  - [ ] 5.8 Open 3 Escher maps in browser: Verify maps load, display flux colors/widths correctly
  - [ ] 5.9 Check validation table in notebook: Document any failures, investigate causes
  - [ ] 5.10 Calculate expression-flux correlation for central metabolism pathways
  - [ ] 5.11 Add final summary cell: Mean biomass by strain, DGOA flux statistics, condition comparison
  - [ ] 5.12 Document edge cases: Any strains with warnings, missing genes, low coverage
  - [ ] 5.13 Save notebook with all outputs visible
  - [ ] 5.14 Create backup: `cp ADP1ExpressionAnalysis.ipynb ADP1ExpressionAnalysis_$(date +%Y%m%d).ipynb`
