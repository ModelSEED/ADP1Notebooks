# ModelReconciliation Notebook - Development Summary

## Purpose

Compare the published ADP1 metabolic model against the ModelSEED-reconstructed model, identify unique published reactions, merge them into the ModelSEED model with ATP safety filtering, and visualize the merged model on an Escher map with reaction class badges.

Both models have been translated to the same namespace (ModelSEED compound/reaction IDs), enabling direct reaction-by-reaction comparison of flux, genes, stoichiometry, and directionality.

## Files

### Notebook
- `notebooks/ModelReconciliation.ipynb` - 9 code cells (plus markdown headers)

### Model Files
- **Published model**: `notebooks/models/FullyTranslatedPublishedModel.json`
  - Local COBRApy JSON, loaded via `MSModelUtil.from_cobrapy()`
  - Biomass reaction: `GROWTH_DASH_RXN`
  - Exchange reactions use `EXF_DASH_*_LBRACKET_Extraorganism_RBRACKET_` naming
  - ~990 reactions, gene IDs in `ACIAD_RS*` format
  - Compartments: `Cytosol`, `Extraorganism`, `Periplasm`, `e`
  - Created by translating the original published XML model (see `ADP1ModelReview.ipynb` Cell 1)

- **ModelSEED model**: KBase workspace ref `179225/Abaylyi_ADP1_RASTMS2_OMEGGA_Abaylyi_ADP1_RAST.mdlMS2_OMEGGA_iAbaylyi_Carbon_Succinic.gf`
  - Loaded via `util.get_model()` from KBase
  - Also saved locally as `notebooks/models/ModelSEED_ADP1.json` (Cell 5)
  - Biomass reaction: `bio1`
  - Exchange reactions use `EX_cpd*_e0` naming
  - ~1417 reactions, gene IDs in `ACIAD_RS*` format
  - Has `mRNA_` prefixed genes that must be removed after loading

- **Merged model**: `notebooks/models/MergedADP1Model.json`
  - ModelSEED model + retained published reactions (after ATP safe expansion filtering)
  - Created in Cell 8

### Datacache Files (in `notebooks/datacache/ModelReconciliation/`)
- `model_info.json` - Basic model statistics
- `published_pyr_results.json` - pFBA fluxes and FVA results for published model in pyruvate media
- `modelseed_pyr_results.json` - pFBA fluxes and FVA results for ModelSEED model in pyruvate media
- `comparison_dataframe.json` - Full comparison dataframe as list of records
- `published_only_unique_gene_rxns.json` - Published-only reactions with genes absent from ModelSEED
- `merge_results.json` - Summary of merge + ATP filtering (retained/filtered reaction IDs, counts)
- `merged_pyr_results.json` - pFBA fluxes and FVA results for merged model in pyruvate media
- `reaction_badges.json` - Reaction class badge assignments for Escher map

### Output Files
- `notebooks/nboutput/model_reconciliation.tsv` - TSV export of comparison dataframe
- `notebooks/nboutput/ModelReconciliation/merged_model_escher.html` - Escher map with badges

## Notebook Cell Structure

### Cell 1: Load Both Models
- Loads published model from local JSON via `MSModelUtil.from_cobrapy("models/FullyTranslatedPublishedModel.json")`
- Loads ModelSEED model from KBase via `util.get_model("179225/...")`
- Removes `mRNA_` prefixed genes from ModelSEED model
- Saves `ModelReconciliation/model_info` to datacache

### Cell 2: Run pFBA and FVA
- Runs `util.run_fba()` (pFBA) on both models with `media="KBaseMedia/Carbon-Pyruvic-Acid"`
- Runs `util.run_fva()` at `fraction_of_optimum=0.5` on both models
- Published model objective: `MAX{GROWTH_DASH_RXN}`
- ModelSEED model objective: `MAX{bio1}`
- Saves results to `ModelReconciliation/published_pyr_results` and `ModelReconciliation/modelseed_pyr_results`

### Cell 3: Build Comparison Dataframe
- Loads cached pFBA/FVA results and reloads models for structure inspection
- Standardizes exchange reaction IDs across models to `EX_cpd<ID>` format
- Builds comprehensive dataframe with columns: ID, published_pyr_flux, modelseed_pyr_flux, membership, shared_genes, extra_published_genes, extra_modelseed_genes, directionality, equation, stoichiometry_differences
- Saves to datacache and TSV

### Cell 4: Identify Published-Only Reactions with Non-Overlapping Genes
- Filters published model reactions to find those that:
  1. Are NOT in the ModelSEED model
  2. Have genes that do NOT appear anywhere in the ModelSEED model
  3. Are NOT exchange reactions (single-metabolite) or spontaneous (no genes)
- Saves list to `ModelReconciliation/published_only_unique_gene_rxns`

### Cell 5: Save ModelSEED Model Locally as JSON
- Pulls ModelSEED model from KBase, removes `mRNA_` genes
- Saves to `models/ModelSEED_ADP1.json` for local loading in subsequent cells

### Cell 6: Merge Published Reactions and Run ATP Safe Expansion
- Loads local ModelSEED model and published model
- Excludes diffusion/transport reactions using `util.is_diffusion_reaction()` (reactions that move metabolites between compartments with no chemical change)
- Adds remaining candidate reactions from Cell 4 to the ModelSEED model
- Also adds exchange reactions for any new extracellular metabolites
- Uses `MSATPCorrection` to:
  1. Evaluate ATP production on default media (no gapfilling)
  2. Build ATP threshold tests
  3. Run `reaction_expansion_test` on the added reactions
- Removes reactions that break ATP production thresholds
- Saves merge statistics to `ModelReconciliation/merge_results`

### Cell 7: Run pFBA and FVA on Merged Model
- Rebuilds merged model from local JSON + retained reactions
- Runs pFBA and FVA (50% optimum) in pyruvate media
- Compares growth rates: published vs ModelSEED vs merged
- Classifies all reactions by FVA class
- Saves to `ModelReconciliation/merged_pyr_results`

### Cell 8: Save the Final Merged Model
- Rebuilds merged model and saves as `models/MergedADP1Model.json`

### Cell 9: Escher Map with Reaction Class Badges
- Loads merged model and pFBA/FVA results
- Classifies each reaction into one of 6 badge categories:
  - **Essential (ModelSEED)**: Must carry flux, from original ModelSEED model
  - **Variable (ModelSEED)**: Flexible flux, from original ModelSEED model
  - **Blocked (ModelSEED)**: No flux possible, from original ModelSEED model
  - **Essential (Published)**: Must carry flux, added from published model
  - **Variable (Published)**: Flexible flux, added from published model
  - **Blocked (Published)**: No flux possible, added from published model
- Uses `util.create_map_html2()` from KBUtilLib's `EscherUtils` to generate interactive HTML
- Flux data painted on map, badges injected as colored SVG rectangles
- Auto-selects best map by model reaction coverage

## Utility Functions in util.py

All added as methods of `NotebookUtil` class in `notebooks/util.py`:

### `classify_fva_flux(fva_result, flux_value, zero_tol=1e-5)`
Classifies a reaction based on FVA min/max and pFBA flux into one of:
`blocked`, `essential_fwd`, `essential_rev`, `optionally_active_fwd`, `optionally_active_rev`, `variable_fwd`, `variable_rev`, `variable_zero`, `optionally_active_zero_fwd`, `optionally_active_zero_rev`

### `standardize_exchange_id(rxn)`
Returns `EX_<cpd_id>` for single-metabolite reactions, original ID otherwise.

### `get_exchange_map(model)`
Returns `{standardized_exchange_id: original_rxn_id}` for all single-metabolite reactions.

### `compare_reaction_stoichiometry(rxn1, rxn2)`
Returns list of difference strings comparing metabolite coefficients between two reactions.

### `get_reaction_directionality(rxn)`
Returns `"forward"`, `"reverse"`, or `"reversible"` based on reaction bounds.

### `is_diffusion_reaction(rxn)`
Returns `True` if a reaction is a pure single-compound diffusion reaction — exactly one unique compound (by name) moving between compartments with no co-transported species. H+-coupled symporters (e.g., `h+ + hexanoate <=> h+ + hexanoate`) are NOT flagged.

### `build_gene_reaction_map(model)`
Returns `{gene_id: [rxn_id, ...]}` mapping, filtering out `mRNA_` prefixed genes.

## Key Design Decisions

### ATP Safe Expansion (Cell 6)
Uses `MSATPCorrection` from ModelSEEDpy to ensure added reactions don't create spurious ATP production. The workflow:
1. Temporarily disables all candidate reactions (bounds set to 0)
2. Evaluates ATP production on default media without gapfilling
3. Determines which media produce ATP reliably
4. Builds threshold tests (max ATP allowed = multiplier × baseline ATP)
5. Re-enables reactions one-by-one, filtering any that exceed ATP thresholds
6. Uses `reaction_expansion_test` with binary search for efficiency

Key method: `MSModelUtil.reaction_expansion_test(reaction_list, condition_list)` — tests each reaction direction independently and filters those that cause ATP to exceed thresholds.

### Diffusion Reaction Filtering
Pure single-compound diffusion reactions are excluded before merging. These reactions move a single metabolite between compartments with no co-transported species (e.g., `formate <=> formate`, `no2- <=> no2-`). Detected programmatically by `util.is_diffusion_reaction()`, which checks that only one unique compound name appears and has net zero stoichiometry. H+-coupled symporters (e.g., `h+ + hexanoate <=> h+ + hexanoate`) are NOT excluded since the H+ co-transport is meaningful.

### Manual Reaction Exclusions
- `rxn08856_c0`: Excluded because it creates a flux loop with existing reaction `rxn08854_c0`. Both reactions carried -47.354 flux when merged, an emergent coupling not present in either parent model.

### Exchange Reaction Standardization
Both models use different exchange naming conventions, standardized to `EX_cpd<ID>` for comparison.

### Model Rebuilding Pattern
Cells 7-9 rebuild the merged model from components (base JSON + retained reaction list from datacache) rather than serializing the modified model from Cell 6. This maintains cell independence.

### Escher Map Badges
Uses `EscherUtils.create_map_html2()` with `reaction_classes` parameter. Badges are injected as colored SVG rectangles via JavaScript post-processing of the Escher HTML. A legend is added in the bottom-right corner.

## Observed Results (from last run)
- Published model: 990 reactions, pFBA growth rate 0.446864
- ModelSEED model: 1417 reactions, pFBA growth rate 0.209254
- Total unique canonical reactions: 1804
  - Published only: 387
  - ModelSEED only: 814
  - In both: 603

## Dependencies

### Libraries
- `modelseedpy` (MSModelUtil, MSPackageManager, MSMedia, MSATPCorrection, etc.)
- `cobra` (COBRApy for model I/O and FBA)
- `cobrakbase` (KBase model format support)
- `kbutillib` (MSFBAUtils, NotebookUtils, EscherUtils, etc.)
- `escher` (Escher map visualization)
- `pandas` (DataFrame construction and export)

### External Services
- KBase workspace API (for loading ModelSEED model and media)
- Media: `KBaseMedia/Carbon-Pyruvic-Acid` (pyruvate minimal media from KBase)
