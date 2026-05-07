# Phase 4c-ii Migration Log: Expression / Annotation Trio

**Date**: 2026-05-06
**Scope**: ADP1AnnotationAnalysis, ADP1ExpressionAnalysis, ADP1AIGeneAssociation

## ADP1AnnotationAnalysis.ipynb (31 cells total, 16 code cells migrated)

### Cells migrated
- **Cell 1**: KBase annotation API call. `util.anno_client()` + `AnnotationOntology.from_kbase_data()` → `_legacy` shim. `util.save` → `session.cache.save`.
- **Cell 3**: `util.get_model_and_simulate()` → `_legacy.get_model_and_simulate()`. Cache load/save migrated.
- **Cell 5**: Same pattern as cell 3 for published model.
- **Cell 7**: Pure pandas (Tukey proteomics). `util.load/save` → `session.cache`.
- **Cell 9**: Dictionary merge (EC/KO/PF/SSO). Cache only.
- **Cell 11**: Pure Python term-name replacement. Cache only.
- **Cell 13**: Pandas → Excel export. Cache load only.
- **Cell 15**: **Canonical ADP1Genome producer.** `util.get_object(id, ws)` → `_legacy.get_object(id, ws)` + `session.cache.save("ADP1Genome", data)`. Plain key, no slash prefix.
- **Cell 17**: Genome + annotation → CSV. `util.load` → `session.cache.load`. Contains inline `convert_role_to_searchrole` helper (pure regex, left inline).
- **Cell 18**: GFF editing. No util calls — pure pandas + file I/O. Fixed leading whitespace on import.
- **Cell 19**: Proteomics comparison → Excel. `util.load` → `session.cache.load`.
- **Cell 20**: Strain metadata dict. `util.save` → `session.cache.save`.
- **Cell 21**: Cluster sizes from text file. `util.save` → `session.cache.save`.
- **Cell 22**: Cluster + protein assignments merge. Cache load/save.
- **Cell 23**: Reaction/flux data merge into gene hash. Cache load/save.
- **Cell 24**: Raw proteome values. Cache load/save.
- **Cell 26**: Comprehensive gene dataframe build. `util.load` → `session.cache.load`.
- **Cell 28**: Pairwise comparisons — pure pandas, no util calls. Unchanged.
- **Cell 30**: PLM API queries. `util.query_plm_api()` → `_legacy.query_plm_api()`. Cache load/save.

### Deferrals
- `convert_role_to_searchrole` left inline in cell 17 (pure regex helper, only used once; not worth porting to util.py).

### Deprecations
- None. All original functionality preserved.

## ADP1ExpressionAnalysis.ipynb (16 cells total, 8 code cells migrated)

### Cells migrated
- **Cell 1**: Setup cell. Added `_legacy = NotebookUtil()` alongside `%run util.py`.
- **Cell 3**: Strain list definition. Unchanged (no util calls).
- **Cell 5**: `util.run_expression_flux_analysis(...)` → `_legacy.run_expression_flux_analysis(...)`. This method internally calls `get_media`, `process_strain_with_expression`, etc. — all via the `_legacy` instance. Media application (Antipattern 2) is handled inside the legacy method.
- **Cell 7**: `util.validate_expression_flux_solution(...)` → `_legacy.validate_expression_flux_solution(...)`.
- **Cell 9**: `util.save(f"{key}_fluxes", ...)` → `session.cache.save(...)`. Plain key, no slash prefix.
- **Cell 10**: `util.create_expression_flux_summary(...)` → `_legacy.create_expression_flux_summary(...)` + `session.cache.save`.
- **Cell 11**: `util.export_expression_flux_to_excel(...)` → `_legacy.export_expression_flux_to_excel(...)`.
- **Cell 12**: `util.generate_all_escher_maps(...)` → `_legacy.generate_all_escher_maps(...)`.
- **Cell 14**: Summary stats — pure pandas. Unchanged.

### Deferrals
- The `run_expression_flux_analysis` method in util_legacy.py still references `GROWTH_DASH_RXN` internally (line 219 of `process_strain_with_expression`). This is Antipattern 1 inside the legacy code. NOT fixed here — it's inside the legacy shim which will be deprecated in Phase 4d. The caller passes the correct `media_id` and the legacy method handles the rest.

### Deprecations
- None.

## ADP1AIGeneAssociation.ipynb (3 cells total, 1 code cell migrated)

### Cells migrated
- **Cell 1**: Model load + gene association evaluation. `MSModelUtil.from_cobrapy` (stays), `util.load` → `session.cache.load`, `util.evaluate_reaction_gene_association(...)` → `_legacy.evaluate_reaction_gene_association(...)` (AICurationUtils method), `util.save` → `session.cache.save`.
- **Cell 2**: Empty. Unchanged.

### Deferrals
- None.

### Deprecations
- None.

## util.py additions

**None.** No new functions were ported. All KBase-touching methods are accessed via `_legacy`. The existing helpers (from Phase 4a/4b) already cover the pure-Python needs of these notebooks.

## Antipattern audit results

| # | Antipattern | Status |
|---|---|---|
| 1 | Hardcoded `GROWTH_DASH_RXN` | Not present in migrated cells. Present inside `_legacy.process_strain_with_expression` (deferred to Phase 4d). |
| 2 | Missing media application | N/A — all FBA-calling cells delegate to `_legacy` which handles media internally. |
| 3 | Algorithm reimplementations | No new util.py functions added. All complex methods accessed via `_legacy`. |
| 4 | Slash-prefixed cache keys | None found. All cache keys are plain (e.g., `"ADP1Genome"`, `"gene_term_hash"`, `"rxn_gene_mapping"`). |
| 5 | `sol.fluxes.values()` | Not present in any migrated cell. |

## Cross-notebook slash-prefixed loads (consumers to fix later)

**None found.** Searched all `.ipynb` and `.py` files for `ADP1AnnotationAnalysis/`, `ADP1ExpressionAnalysis/`, `ADP1AIGeneAssociation/` — zero hits.

## Notes for playbook section 13

- The AnnotationAnalysis notebook is unusually cell-heavy (16 code cells) but each cell is thin — predominantly cache load/save with inline data wrangling.
- `ExpressionAnalysis` is a "thick _legacy" notebook: almost all logic lives in `NotebookUtil.run_expression_flux_analysis`. Migration is trivial (swap `util.` → `_legacy.`) but the GROWTH_DASH_RXN bug (Antipattern 1) remains inside the legacy code.
- `AIGeneAssociation` is minimal (1 real cell). The `evaluate_reaction_gene_association` method comes from `AICurationUtils` parent class — no reimplementation needed.
