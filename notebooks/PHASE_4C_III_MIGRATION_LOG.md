# Phase 4c-iii Migration Log: Phenotype/Essentiality Trio

Date: 2026-05-06

## ADP1DGOAAnalysis.ipynb

**Original cells:** 22 (17 code)
**Migrated cells:** 18 (9 code)

Consolidated heavily. The original had many exploratory/debugging cells for Escher map generation that were collapsed into a single clean cell.

**Key changes:**
- All `util.load/save` -> `session.cache.load/save` (plain keys, no slash prefix)
- All KBase API calls (get_model, get_media, run_fba, run_fva, add_dgoa_reaction, translate_model_to_ms_namespace) -> `_legacy` shim
- Biomass objective consistently set to `"bio1"` (AP1 fixed)
- Media explicitly applied via `_legacy.set_media(model, media)` or `_legacy.run_fba(model, media=...)` (AP2 addressed)
- No algorithm reimplementations (AP3 clear)
- No slash-prefixed cache keys (AP4 clear)
- No `.fluxes.values()` calls (AP5 clear)
- Escher maps delegated entirely to `_legacy.create_map_html2()` (no reimplementation)

**Deferrals:** None.

**Deprecations:** Cells 15-21 (Escher map debugging/workaround cells) removed -- their functionality is now a single clean cell using `_legacy.create_map_html2`.

## ADP1EssentialityAnalysis.ipynb

**Original cells:** 21 (10 code)
**Migrated cells:** 21 (10 code)

Near 1:1 migration. This notebook is primarily pure-Python statistics (pandas, scipy, matplotlib) with minimal KBase interaction.

**Key changes:**
- All `util.load/save` -> `session.cache.load/save` (plain keys)
- `gene_translation` loaded from cache key `"adp1_gene_translation"` (no slash prefix)
- `essentiality_gene_lists`, `tnseq_parameters`, `essentiality_mutant_overlap`, `essentiality_stats_comparison` all plain keys
- No FBA/FVA in this notebook -- purely statistical analysis of pre-existing growth rate data
- MSExpression used for data loading only (from_spreadsheet), not for model constraint fitting
- matplotlib/scipy usage left as-is (pure Python, not a KBUtilLib method)

**Antipattern scorecard:** All 5 clean (AP1-5: not applicable or explicitly avoided).

**Deferrals:** None.

## ADP1MutantPhenotypeAnalysis.ipynb

**Original cells:** 18 (8 code)
**Migrated cells:** 18 (8 code + 1 empty)

The most complex of the three. Heavy use of `_legacy` for expression-constrained FBA, media handling, and model manipulation.

**Key changes:**
- All `util.load/save` -> `session.cache.load/save`
- Cache keys use flat names: `"ADP1-AA-MGR"`, `"ADP1-NR-MGR"`, `"ADP1-NR-RxnMGR"`, `"ADP1-AA-MGR-conditions"`, `"ADP1-AA-MGR-stats"`, `"ADP1-AA-MGR-media"`, `"ADP1-AA-MGR-growth"`, `"FittingMutantPhenotypeData"` (AP4 clear)
- Biomass objective: `"MAX{bio1}"` for `constrain_objective_to_fraction_of_optimum`, `"bio1"` for `biomass_id` in coupling analysis (AP1 fixed from legacy `GROWTH_DASH_RXN`)
- Media applied via `_legacy.constrain_objective_to_fraction_of_optimum(model, media=media, ...)` which internally sets media (AP2 addressed)
- `MSExpression.fit_flux_to_mutant_growth_rate_data()` called directly (ModelSEEDpy method, not reimplemented) (AP3 clear)
- `_legacy.analyzed_reaction_objective_coupling()` used for KO analysis (not reimplemented) (AP3 clear)
- `_legacy.get_msgenome_from_dict()` for genome object creation (AP3 clear)
- No `.fluxes.values()` -- uses `.fluxes.to_dict()`, `.fluxes.get()` (AP5 clear)
- Escher maps via `_legacy.create_map_html2()` (AP3 clear)
- `compute_expression_flux_statistics()` defined inline in notebook cell (pure Python dict manipulation, no model deps -- not ported to util.py since it's single-use and cell-local)

**Deferrals:** None.

**Cross-notebook dependencies:**
- Loads `"ADP1Genome"` from cache (produced by ADP1AnnotationAnalysis in 4c-ii)
- Loads `"essentiality_gene_lists"` (produced by ADP1EssentialityAnalysis earlier in this trio)
- Produces `"FittingMutantPhenotypeData"` consumed by no other notebooks currently

## util.py additions

None. All three notebooks' logic is either:
- Pure inline data wrangling (statistics, plotting) that stays in the cell
- KBase-service-touching methods delegated to `_legacy`
- ModelSEEDpy methods called directly (MSExpression, MSMedia, MSModelUtil)

No functions crossed the "pure Python helper worth porting" threshold.

## Slash-prefixed loads in OTHER notebooks

Searched for patterns `"ADP1DGOAAnalysis/"`, `"ADP1EssentialityAnalysis/"`, `"ADP1MutantPhenotypeAnalysis/"` in all notebooks -- zero hits found.

## Lessons for playbook section 13

- The `compute_expression_flux_statistics` function in MutantPhenotypeAnalysis is a pure-Python helper but is cell-local and single-use. Decision: leave inline rather than port to util.py. Rule of thumb: only port to util.py if the function is (a) used by multiple notebooks OR (b) >20 lines and benefits from unit testing.
- The `objective="MAX{bio1}"` syntax for `constrain_objective_to_fraction_of_optimum` is a KBUtilLib convention where the `MAX{}` wrapper signals maximization direction. This is NOT the same as passing `objective="bio1"` (which some methods also accept). When migrating from `GROWTH_DASH_RXN`, check whether the legacy call used `MAX{GROWTH_DASH_RXN}` wrapper or bare ID.
- MutantPhenotypeAnalysis uses a compound cache key naming scheme (`"ADP1-AA-MGR"`, `"ADP1-NR-MGR"`) that encodes data transformations. This is a valid flat-key convention (no slash), but future consumers need to know the schema: `<project>-<normalization>-<datatype>`.
