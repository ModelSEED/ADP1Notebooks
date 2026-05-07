"""adp1notebooks util -- runs as %run util.py at the top of every cell.

Provides: common imports, project NotebookSession (named `session`), and
project-specific helper functions (added below the marker by migration tasks).
"""
from __future__ import annotations

# === System path bootstrap =================================================
# Prepend machine-specific Python paths to sys.path BEFORE any heavy imports.
# Reads ~/.kbu-sys-paths (plain text, one path per line, # comments OK).
# Silent no-op if the file isn't present or unreadable.
import sys as _sys
from pathlib import Path as _Path

def _bootstrap_sys_paths() -> None:
    user_file = _Path.home() / ".kbu-sys-paths"
    if not user_file.exists():
        return
    try:
        for raw in user_file.read_text().splitlines():
            s = raw.split("#", 1)[0].strip()
            if not s:
                continue
            expanded = str(_Path(s).expanduser())
            if expanded and expanded not in _sys.path:
                _sys.path.insert(0, expanded)
    except Exception:
        pass

_bootstrap_sys_paths()
# ===========================================================================

import json, math, os, re
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd
try:
    import cobra
    from cobra import Reaction, Metabolite
except ImportError:
    cobra = None  # type: ignore

from kbutillib.notebook import NotebookSession
from kbutillib.notebook.schema import (
    Sample, Computation, ExternalDataset, Experiment,
    Strain, Mutation, Media,
    Vector, VectorType, EntityKind, EntityRef,
)

from modelseedpy import MSModelUtil

# Auto-detect notebook name from IPython runtime.
session: NotebookSession = NotebookSession.for_notebook(
    None, project_name="adp1notebooks"
)


def session_for(notebook_file: Optional[str] = None) -> NotebookSession:
    """Back-compat shim. Prefer the module-level `session`."""
    return NotebookSession.for_notebook(notebook_file, project_name="adp1notebooks")


# === project-specific helpers below ===

# ---------------------------------------------------------------------------
# Compartment utilities
# ---------------------------------------------------------------------------

COMPARTMENT_MAP: dict[str, str] = {
    "c": "cytoplasm",
    "c0": "cytoplasm",
    "e": "extracellular",
    "e0": "extracellular",
    "p": "periplasm",
    "p0": "periplasm",
    "m": "mitochondria",
    "m0": "mitochondria",
}


def normalize_compartment(compartment_id: str) -> str:
    """Map a compartment suffix/id to a canonical human-readable name.

    Returns the original id unchanged if no mapping is found.
    """
    return COMPARTMENT_MAP.get(compartment_id, compartment_id)


# ---------------------------------------------------------------------------
# Reaction directionality
# ---------------------------------------------------------------------------

def get_reaction_directionality(reaction: cobra.Reaction) -> str:
    """Classify reaction directionality based on its bounds.

    Returns one of: 'reversible', 'forward', 'reverse', 'blocked'.
    """
    if reaction.lower_bound < 0 and reaction.upper_bound > 0:
        return "reversible"
    elif reaction.lower_bound >= 0 and reaction.upper_bound > 0:
        return "forward"
    elif reaction.lower_bound < 0 and reaction.upper_bound <= 0:
        return "reverse"
    else:
        return "blocked"


# ---------------------------------------------------------------------------
# Exchange reactions
# ---------------------------------------------------------------------------

_EXCHANGE_RE = re.compile(r"^EX_(.+?)(?:_e0?)?$")


def standardize_exchange_id(reaction_id: str) -> str:
    """Normalize exchange reaction IDs to the form 'EX_<met_id>_e0'.

    If *reaction_id* does not look like an exchange reaction, return it as-is.
    """
    m = _EXCHANGE_RE.match(reaction_id)
    if m:
        met_base = m.group(1)
        # Strip trailing compartment tag if already present
        met_base = re.sub(r"_e0?$", "", met_base)
        return f"EX_{met_base}_e0"
    return reaction_id


def get_exchange_map(model: cobra.Model) -> dict[str, cobra.Reaction]:
    """Return a mapping of standardized exchange id -> Reaction for all exchanges."""
    result: dict[str, cobra.Reaction] = {}
    for rxn in model.exchanges:
        std_id = standardize_exchange_id(rxn.id)
        result[std_id] = rxn
    return result


# ---------------------------------------------------------------------------
# Gene-reaction mapping
# ---------------------------------------------------------------------------

def build_gene_reaction_map(model: cobra.Model) -> dict[str, list[str]]:
    """Build a dict mapping gene id -> list of reaction ids that use that gene."""
    gene_rxn_map: dict[str, list[str]] = {}
    for gene in model.genes:
        gene_rxn_map[gene.id] = [rxn.id for rxn in gene.reactions]
    return gene_rxn_map


# ---------------------------------------------------------------------------
# Reaction display
# ---------------------------------------------------------------------------

def reaction_equation_with_names(reaction: cobra.Reaction) -> str:
    """Build a human-readable reaction equation using metabolite names."""
    def _format_met(met: cobra.Metabolite, coeff: float) -> str:
        abs_coeff = abs(coeff)
        name = met.name if met.name else met.id
        if abs_coeff == 1.0:
            return name
        # Use integer display when possible
        if abs_coeff == int(abs_coeff):
            return f"{int(abs_coeff)} {name}"
        return f"{abs_coeff} {name}"

    reactants = []
    products = []
    for met, coeff in reaction.metabolites.items():
        if coeff < 0:
            reactants.append(_format_met(met, coeff))
        else:
            products.append(_format_met(met, coeff))

    arrow = " <=> " if reaction.reversibility else " => "
    return " + ".join(reactants) + arrow + " + ".join(products)


# ---------------------------------------------------------------------------
# Diffusion detection
# ---------------------------------------------------------------------------

def is_diffusion_reaction(reaction: cobra.Reaction) -> bool:
    """Return True if the reaction is a simple diffusion (same metabolite in two compartments)."""
    mets = list(reaction.metabolites.keys())
    if len(mets) != 2:
        return False
    # Check same base metabolite in different compartments
    id0 = re.sub(r"_[a-z]\d?$", "", mets[0].id)
    id1 = re.sub(r"_[a-z]\d?$", "", mets[1].id)
    if id0 != id1:
        return False
    # Coefficients should be +1 and -1
    coeffs = sorted(reaction.metabolites.values())
    return coeffs == [-1.0, 1.0] or coeffs == [-1, 1]


# ---------------------------------------------------------------------------
# Stoichiometry comparison
# ---------------------------------------------------------------------------

def compare_reaction_stoichiometry(
    rxn_a: cobra.Reaction,
    rxn_b: cobra.Reaction,
) -> dict[str, Any]:
    """Compare stoichiometry of two reactions.

    Returns a dict with keys:
        - identical: bool
        - only_in_a: list of met ids
        - only_in_b: list of met ids
        - coefficient_diffs: dict of met_id -> (coeff_a, coeff_b)
    """
    mets_a = {met.id: coeff for met, coeff in rxn_a.metabolites.items()}
    mets_b = {met.id: coeff for met, coeff in rxn_b.metabolites.items()}

    only_in_a = sorted(set(mets_a) - set(mets_b))
    only_in_b = sorted(set(mets_b) - set(mets_a))

    coefficient_diffs: dict[str, tuple[float, float]] = {}
    for met_id in set(mets_a) & set(mets_b):
        if not np.isclose(mets_a[met_id], mets_b[met_id]):
            coefficient_diffs[met_id] = (mets_a[met_id], mets_b[met_id])

    identical = not only_in_a and not only_in_b and not coefficient_diffs
    return {
        "identical": identical,
        "only_in_a": only_in_a,
        "only_in_b": only_in_b,
        "coefficient_diffs": coefficient_diffs,
    }


# ---------------------------------------------------------------------------
# FVA analysis helpers
# ---------------------------------------------------------------------------

def find_significant_differences(
    fva_df: pd.DataFrame,
    threshold: float = 1e-6,
) -> pd.DataFrame:
    """Filter an FVA result DataFrame to rows where min and max flux differ significantly.

    Expects columns 'minimum' and 'maximum' (standard cobra FVA output).
    """
    span = (fva_df["maximum"] - fva_df["minimum"]).abs()
    return fva_df[span > threshold].copy()


def classify_fva_flux(
    minimum: float,
    maximum: float,
    tol: float = 1e-6,
) -> str:
    """Classify flux variability for a single reaction.

    Returns one of:
        'fixed_zero'     — both bounds effectively zero
        'fixed_nonzero'  — flux is fixed at a nonzero value
        'variable'       — flux can vary (min != max)
        'blocked'        — both bounds are exactly zero (same as fixed_zero but used
                           when there is genuinely no flux possible)
    """
    min_zero = abs(minimum) <= tol
    max_zero = abs(maximum) <= tol
    span = abs(maximum - minimum)

    if min_zero and max_zero:
        return "blocked"
    elif span <= tol:
        return "fixed_nonzero"
    else:
        return "variable"


# ---------------------------------------------------------------------------
# Gene ID translation (ported from util_legacy.py)
# ---------------------------------------------------------------------------

def translate_expression_gene_ids(
    df: pd.DataFrame,
    gene_mapping_file: str = "data/ADP1Genes.csv",
) -> pd.DataFrame:
    """Translate gene IDs from RefSeq format (ACIAD_RSxxxxx) to old format (ACIADxxxx).

    Pure DataFrame transform. The input DataFrame should have gene IDs as index.
    Returns a new DataFrame with translated index; duplicates are dropped (first kept).

    Parameters
    ----------
    df : pd.DataFrame
        Expression data with RefSeq gene IDs as index.
    gene_mapping_file : str
        Path to CSV with columns 'ID' (RefSeq), 'Old ID', 'Name'.
    """
    gene_mapping_df = pd.read_csv(gene_mapping_file)

    refseq_to_old: dict[str, str] = {}
    name_to_old: dict[str, str] = {}

    for _, row in gene_mapping_df.iterrows():
        refseq_id = row.get("ID")
        old_id = row.get("Old ID")
        name = row.get("Name")

        if pd.notna(refseq_id) and pd.notna(old_id):
            refseq_to_old[str(refseq_id)] = str(old_id)
        if pd.notna(name) and pd.notna(old_id) and str(name).strip():
            name_to_old[str(name)] = str(old_id)

    translated_ids = []
    for gene_id in df.index:
        gid = str(gene_id)
        if gid in refseq_to_old:
            translated_ids.append(refseq_to_old[gid])
        elif gid in name_to_old:
            translated_ids.append(name_to_old[gid])
        else:
            translated_ids.append(gid)

    out = df.copy()
    out.index = translated_ids
    out = out[~out.index.duplicated(keep="first")]
    return out
