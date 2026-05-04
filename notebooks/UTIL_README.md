# ADP1 Notebooks — util.py

## Overview

`util.py` provides thin, testable helper functions for ADP1 analysis notebooks.
It replaces the legacy 1100-line multi-inheritance god-class with composable
functions built on `kbutillib.notebook.NotebookSession`.

## Usage

```python
from util import session_for, get_reaction_directionality, classify_fva_flux

session = session_for()  # auto-detects notebook context
```

## Design Principles

1. **Functions, not classes** — each helper is a standalone function.
2. **NotebookSession for state** — caching, provenance, and experiment
   tracking go through the session object.
3. **No pickle** — all serialization uses JSON/parquet via the session cache.
4. **Testable** — every function works on plain cobra objects or DataFrames.

## Running Tests

```bash
cd notebooks/
pytest -v test_util.py
```

## Migration

- `util_legacy.py` is preserved for reference during the Phase 4 refactor.
- Notebooks will be updated in Phase 4b to import from the new `util.py`.
- `util_legacy.py` will be removed in Phase 4d.

## Functions

| Function | Purpose |
|----------|---------|
| `session_for` | Create a NotebookSession for ADP1 notebooks |
| `normalize_compartment` | Map compartment IDs to readable names |
| `get_reaction_directionality` | Classify reaction bounds |
| `standardize_exchange_id` | Normalize exchange reaction IDs |
| `get_exchange_map` | Build exchange-id-to-reaction dict |
| `build_gene_reaction_map` | Map gene IDs to their reactions |
| `reaction_equation_with_names` | Human-readable equation string |
| `is_diffusion_reaction` | Detect simple transport reactions |
| `compare_reaction_stoichiometry` | Diff two reactions' metabolites |
| `find_significant_differences` | Filter FVA results by flux span |
| `classify_fva_flux` | Categorize FVA variability |
