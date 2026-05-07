# Phase 4d Cleanup Log — ADP1Notebooks

Date: 2026-05-06
Executor: AgentForge worker emailmac-2

---

## Section 1 — Dropped Notebooks

| Notebook | Rationale |
|----------|-----------|
| `ADP1ExpressionAnalysis_20251024.ipynb` | Dated copy; superseded by `ADP1ExpressionAnalysis.ipynb` |
| `ADP1FoldChangeAnalysis.ipynb` | Non-BERDL lactate analysis; superseded by `ADP1BERDLFoldChangeAnalysis.ipynb` |
| `ADP1OldExpressionAnalysis.ipynb` | Explicitly marked "Old" |
| `BERDLMockup.ipynb` | Mockup notebook, never used for real analysis |

**Dependency check:** No surviving notebook calls `session.cache.load("ADP1FoldChangeAnalysis_*")` or
otherwise programmatically depends on the deleted notebooks' outputs. The only references are
documentary comments in `ADP1BERDLFoldChangeAnalysis.ipynb` (lines 41, 801) noting historical context.

---

## Section 2 — Legacy Scaffolding

### util_legacy.py — Reduced to Thin Shim (NOT deleted)

**Reason:** 12 of 14 remaining notebooks still contain `from util_legacy import NotebookUtil`
(the `_legacy` shim cells from migration phases 4b/4c). Deleting the file would break them.

Notebooks with `_legacy` imports:
1. ADP1AIGeneAssociation.ipynb
2. ADP1AIModelCompare.ipynb
3. ADP1AnnotationAnalysis.ipynb
4. ADP1BERDLCrossSampleAnalysis.ipynb
5. ADP1BERDLFitnessFluxFitting.ipynb
6. ADP1BERDLFoldChangeAnalysis.ipynb
7. ADP1DGOAAnalysis.ipynb
8. ADP1ExpressionAnalysis.ipynb
9. ADP1ModelReview.ipynb
10. ADP1MutantPhenotypeAnalysis.ipynb
11. ModelDebugging.ipynb
12. ModelReconciliation.ipynb

**Action taken:** Replaced the 1133-line god-class with a ~320-line thin shim that:
- Inherits from the same KBUtilLib mixin classes (`MSFBAUtils, AICurationUtils, NotebookUtils, KBPLMUtils, EscherUtils`)
- Preserves all custom methods not found in KBUtilLib (DGOA helpers, expression processing, FVA classification, model comparison utilities)
- Exposes the module-level `util = NotebookUtil()` for backward compat
- Removes dead code, incomplete pipeline methods, and duplicated logic already in base classes

**Test:** `test_util_legacy.py` written to verify import and instantiation. Cannot run on email-mac (project venv not present).

### datacache/ — Deleted

All 21 JSON files in `notebooks/datacache/` were `git rm`-ed. These were legacy cache outputs.

**Surviving references:** `ADP1BERDLFoldChangeAnalysis.ipynb` (lines 306-307) contains a
try/fallback path checking `datacache/ADP1BERDLFitnessFluxFitting/berdl_constrained_pyruvate_result.json`
and `datacache/berdl_constrained_pyruvate_result.json`. This is a graceful fallback in already-migrated
code (it tries `.kbcache/` first via `session.cache.load`, then falls back to datacache). Since the
notebook is migrated and this is a non-modifiable file, the fallback will simply fail silently and the
notebook will proceed with the cache miss path. No action required.

---

## Section 3 — Papermill Smoke Results

**Status: SKIPPED — environment unavailable**

The ADP1Notebooks project venv (`kbu.nb-adp1notebooks-py3.10`) is not installed on
email-mac. Only `AgentForge-py3.12`, `ClaudeCommands-py3.12`, and `EmailAssistant-py3.12`
are present. The `kbutillib` package is not available in any installed environment.

Papermill smoke tests must be run on the primary laptop where the project venv exists.

**Notebooks that would be tested (14):**
1. ADP1AIGeneAssociation.ipynb — likely requires KBase auth
2. ADP1AIModelCompare.ipynb — likely requires KBase auth
3. ADP1AnnotationAnalysis.ipynb
4. ADP1BERDLAnalysis.ipynb
5. ADP1BERDLCrossSampleAnalysis.ipynb — likely requires KBase auth
6. ADP1BERDLFitnessFluxFitting.ipynb — likely requires KBase auth
7. ADP1BERDLFoldChangeAnalysis.ipynb — likely requires KBase auth
8. ADP1DGOAAnalysis.ipynb
9. ADP1EssentialityAnalysis.ipynb
10. ADP1ExpressionAnalysis.ipynb — likely requires KBase auth
11. ADP1ModelReview.ipynb — likely requires KBase auth
12. ADP1MutantPhenotypeAnalysis.ipynb — likely requires KBase auth
13. ModelDebugging.ipynb — likely requires KBase auth
14. ModelReconciliation.ipynb — likely requires KBase auth

---

## Section 4 — Open Follow-ups

1. **ModelDebugging.ipynb cells 2-3 use bare `util.get_media()` / `util.run_fba()` via `%run util.py`.**
   The new `util.py` does not expose a module-level `util` object with those methods — it exposes
   `session` (a `NotebookSession`). Cells 2 and beyond that use `util.<method>` will fail at
   runtime. This was already flagged by the 4c-iv reviewer. Fix in Phase 5 (either add
   compatibility shim to util.py or update the notebook cells).

2. **12 notebooks retain `_legacy` shim cells.** Phase 5 should remove these cells and delete
   `util_legacy.py` once all notebooks are fully migrated to the `session`-based API.

3. **Papermill smoke tests not executed.** Must be run on primary-laptop with the correct venv.
   Recommend adding a CI job or follow-up task.

4. **`ADP1BERDLFoldChangeAnalysis.ipynb` datacache fallback path** (lines 306-307) references
   deleted `datacache/` files. The code is structured to handle the miss gracefully, but the
   dead fallback paths should be cleaned up in Phase 5 if the notebook is touched.

5. **`test_util_legacy.py` not validated.** The shim test was written but cannot execute on
   email-mac. Must be run on primary-laptop.
