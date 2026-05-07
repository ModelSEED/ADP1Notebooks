# Phase 4c-iv Migration Log: Model Review and Debugging Quartet

Date: 2026-05-06

## ADP1ModelReview.ipynb

**Original cells:** 10 (6 code, 4 markdown/empty)
**Migrated cells:** 8 (5 code, 3 markdown)

Consolidated from 10 cells to 8. Removed two empty code cells. Three code cells
with legacy `util.*` calls migrated to `_legacy` shim pattern.

**Key changes:**
- All `util.save/load` -> `session.cache.save/load` (plain keys: `gene_translation`, `match_stats`, `model_match_stats`, `mdl_rxn_hash`)
- KBase API calls (`load_kbase_gene_container`, `alias_to_ftrs`, `get_model`, `model_standardization`, `_parse_rxn_stoichiometry`, `ai_analysis_of_model_reactions`, `add_dgoa_reaction`, `get_media`, `set_media`, `display_dataframe`) -> `_legacy` shim
- Cell 5 (directionality flux impact): replaced `model.set_media(pyruvate_media)` with `_legacy.set_media(model, pyruvate_media)` (AP2 fix)
- Cell 5: already used `"bio1"` for biomass objective, no AP1 fix needed
- Cell 5: replaced `wt_solution.fluxes[rxn.id]` access with `.fluxes.get(rxn.id, 0.0)` (AP5 safety)

**Antipattern scorecard:** AP1 clean (already used `bio1`), AP2 fixed (media set via `_legacy.set_media`), AP3 clean (no reimplementations), AP4 clean (no slash prefixes), AP5 fixed (`.get()` instead of direct index)

**Deferrals:** None.

## ADP1AIModelCompare.ipynb

**Original cells:** 5 (3 code, 1 markdown, 1 empty)
**Migrated cells:** 4 (3 code, 1 markdown)

Removed one empty cell. All three code cells use AI curation methods via `_legacy`.

**Key changes:**
- `util.load/save` -> `session.cache.load/save` (keys: `reactiondata`, `rxn_mapping`)
- KBase/AI API calls (`get_reaction_by_id`, `evaluate_reaction_equivalence`, `reaction_id_to_msid`, `_load_cached_curation`) -> `_legacy` shim
- No biomass/FBA/FVA calls in this notebook -> AP1, AP2 not applicable
- No `.fluxes` usage -> AP5 not applicable

**Antipattern scorecard:** All 5 clean (AP1-5: not applicable or explicitly avoided)

**Deferrals:** None.

## ModelDebugging.ipynb

**Original cells:** 23 (15 code, 8 markdown/empty)
**Migrated cells:** 23 (15 code, 8 markdown)

Preserved 1:1 cell structure. This notebook is primarily pure-Python model comparison
(pandas DataFrames, cobra model attributes). Minimal KBase interaction.

**Key changes:**
- Cells 6-8 (FBA tests): replaced `objective="GROWTH_DASH_RXN"` with `objective="MAX{bio1}"` (AP1 fix)
- Cells 6-8: `util.get_media(...)` -> `_legacy.get_media(...)`, `util.run_fba(...)` -> `_legacy.run_fba(...)` with media parameter (AP2 addressed)
- Cells 10-22 (model comparison): pure Python/cobra/pandas, no migration needed beyond preamble
- Cell 1: preamble now includes `_legacy` initialization for cells that need it

**Antipattern scorecard:** AP1 fixed (3 instances of `GROWTH_DASH_RXN` -> `bio1`), AP2 addressed (media via `_legacy.run_fba(media=...)`), AP3 clean, AP4 clean, AP5 clean (no flux iteration)

**Deferrals:** None.

## ModelReconciliation.ipynb

**Original cells:** 20 (10 code, 10 markdown)
**Migrated cells:** 19 (10 code, 9 markdown)

Removed one empty trailing cell. The most complex notebook in this quartet with
extensive model merging, ATP expansion testing, FVA classification, and Escher map generation.

**Key changes:**
- All `util.save/load` -> `session.cache.save/load`
- Cache keys use `ModelReconciliation/` prefix (AP4 exception: intra-notebook save+load pattern)
- Cell 2 (original cell 1 header): replaced `"biomass_rxn": "GROWTH_DASH_RXN"` with `"bio1"` (AP1 fix)
- Cell 4 (pFBA/FVA): replaced `objective="MAX{GROWTH_DASH_RXN}"` with `objective="MAX{bio1}"` (AP1 fix)
- All FBA/FVA calls go through `_legacy.run_fba/_legacy.run_fva` with media parameter (AP2)
- `util.classify_fva_flux`, `util.get_exchange_map`, etc. -> `_legacy.classify_fva_flux`, `_legacy.get_exchange_map` for legacy dict-style FVA returns
- Cell 6 (comparison dataframe): uses util.py's `build_gene_reaction_map`, `get_reaction_directionality`, `reaction_equation_with_names`, `compare_reaction_stoichiometry` (already ported)
- Cell 8 (published-only reactions): uses `reaction_equation_with_names` from util.py, `_legacy.standardize_exchange_id` for legacy model-method version
- Cell 12 (merge): uses `is_diffusion_reaction` from util.py, `_legacy.normalize_compartment` for legacy compartment mapping
- Cell 17 (Escher): delegates entirely to `_legacy.create_map_html2()`

**Antipattern scorecard:** AP1 fixed (2 `GROWTH_DASH_RXN` -> `bio1`), AP2 addressed (media via `_legacy`), AP3 clean (no reimplementations; ATP expansion uses modelseedpy's MSATPCorrection directly), AP4 exception (intra-notebook `ModelReconciliation/` prefix), AP5 clean (`.fluxes.items()` used correctly)

**Deferrals:** None.

## util.py additions

None. All four notebooks' logic is either:
- Pure inline data wrangling (cobra model attribute comparison, pandas DataFrames) that stays in the cell
- KBase-service-touching methods delegated to `_legacy`
- ModelSEEDpy methods called directly (MSATPCorrection, MSModelUtil)
- Already-ported util.py helpers (build_gene_reaction_map, get_reaction_directionality, reaction_equation_with_names, compare_reaction_stoichiometry, is_diffusion_reaction)

No functions crossed the "pure Python helper worth porting" threshold.

## Slash-prefixed loads in OTHER notebooks

Searched for patterns `"ADP1ModelReview/"`, `"ADP1AIModelCompare/"`, `"ModelDebugging/"`, `"ModelReconciliation/"` in all notebooks. `ModelReconciliation/` keys are used only within ModelReconciliation.ipynb itself (intra-notebook, AP4 exception). No cross-notebook consumers found.

## Lessons for playbook section 13

- The `_legacy.classify_fva_flux()` method in util_legacy.py accepts a dict `{"MIN": ..., "MAX": ...}` as its first argument, while `util.py`'s `classify_fva_flux()` takes separate `minimum` and `maximum` float arguments. When migrating cells that use FVA results from `_legacy.run_fva()` (which returns dict-style FVA), use `_legacy.classify_fva_flux()` for consistency. Only switch to `util.py`'s version when the FVA source also changes to cobra's native DataFrame format.
- ModelReconciliation uses a significant amount of model merging logic (adding reactions, metabolites, exchanges) inline in cells. This code is not a candidate for util.py extraction because it relies heavily on `_legacy.normalize_compartment()` (which uses a different mapping than `util.py`'s `normalize_compartment()`) and modelseedpy's `MSATPCorrection`. The inline pattern is correct here.
- Notebooks that build and compare two models (ModelDebugging, ModelReconciliation) tend to be cell-heavy but mostly pure Python + cobra. Migration is largely a preamble swap with `GROWTH_DASH_RXN` -> `bio1` fixes. Most cells require no structural changes.
