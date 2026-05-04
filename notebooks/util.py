"""ADP1 Notebook utilities — thin helpers over kbutillib.notebook.NotebookSession."""
from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import cobra
import numpy as np
import pandas as pd

from kbutillib.notebook import NotebookSession


# ---------------------------------------------------------------------------
# Session helper
# ---------------------------------------------------------------------------

def session_for(notebook_file: str | None = None) -> NotebookSession:
    """Return a configured NotebookSession for the calling notebook.

    Parameters
    ----------
    notebook_file : str, optional
        Explicit path to the .ipynb file. If *None*, auto-detection is used.
    """
    return NotebookSession.for_notebook(notebook_file, project_name="ADP1Notebooks")


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


# ---------------------------------------------------------------------------
# Fold-change flux simulation (ported from util_legacy.process_strain_with_expression)
# ---------------------------------------------------------------------------

def run_fold_change_simulation(
    model: cobra.Model,
    fold_change_vector: dict[str, float],
    *,
    reference_flux: dict[str, float],
    max_fold_change: float = 3.0,
    max_flux: float = 20.0,
    biomass_reaction: str = "GROWTH_DASH_RXN",
    biomass_fraction: float = 0.25,
    zero_flux: float = 0.001,
) -> dict[str, Any]:
    """Run fold-change-constrained FBA on a cobra model.

    Takes explicit arguments (no hidden state). Returns a dict with:
        - fluxes: dict[str, float]  (only non-zero)
        - biomass: float
        - objective_value: float
        - status: str
        - target_flux: dict[str, float]  (the computed target fluxes)
        - n_active_reactions: int

    Parameters
    ----------
    model : cobra.Model
        A model copy (will be mutated). Caller should provide a deep copy.
    fold_change_vector : dict[str, float]
        Gene-level linear fold-change values (gene_id -> fold_change).
    reference_flux : dict[str, float]
        Reference flux distribution (reaction_id -> flux).
    max_fold_change : float
        Cap on fold change magnitude.
    max_flux : float
        Cap on reference flux magnitude before applying fold changes.
    biomass_reaction : str
        ID of the biomass reaction for constraining.
    biomass_fraction : float
        Minimum fraction of unconstrained optimum to enforce.
    zero_flux : float
        Small flux value for reactions with zero reference flux.
    """
    from cobra.flux_analysis import pfba

    log2_fc_cap = math.log2(max_fold_change)

    # Build gene -> reaction mapping
    gene_rxn_map: dict[str, list[str]] = {}
    for gene in model.genes:
        gene_rxn_map[gene.id] = [rxn.id for rxn in gene.reactions]

    # Compute target fluxes from fold changes
    target_flux: dict[str, float] = {}
    for gene_id, fc in fold_change_vector.items():
        if gene_id not in gene_rxn_map:
            continue
        # Cap the fold change
        log2_fc = math.log2(fc) if fc > 0 else 0.0
        log2_fc = max(-log2_fc_cap, min(log2_fc_cap, log2_fc))
        capped_fc = 2 ** log2_fc

        for rxn_id in gene_rxn_map[gene_id]:
            ref = reference_flux.get(rxn_id, 0.0)
            # Cap reference flux
            ref = max(-max_flux, min(max_flux, ref))
            if abs(ref) < 1e-9:
                ref = zero_flux
            target_flux[rxn_id] = ref * capped_fc

    # Constrain biomass to fraction of optimum
    if biomass_reaction in [r.id for r in model.reactions]:
        # Quick unconstrained solve for optimum
        model.objective = biomass_reaction
        sol0 = model.optimize()
        if sol0.status == "optimal":
            min_biomass = sol0.objective_value * biomass_fraction
            bio_rxn = model.reactions.get_by_id(biomass_reaction)
            bio_rxn.lower_bound = max(bio_rxn.lower_bound, min_biomass)

    # Apply target flux as constraints (soft: adjust bounds toward target)
    for rxn_id, target in target_flux.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
        except KeyError:
            continue
        # Nudge bounds toward target
        if target > 0:
            rxn.lower_bound = max(rxn.lower_bound, min(target * 0.1, rxn.upper_bound))
        elif target < 0:
            rxn.upper_bound = min(rxn.upper_bound, max(target * 0.1, rxn.lower_bound))

    # Run pFBA
    try:
        solution = pfba(model)
        fluxes = {
            rxn_id: flux
            for rxn_id, flux in solution.fluxes.items()
            if abs(flux) > 1e-9
        }
        biomass = solution.fluxes.get(biomass_reaction, 0.0)
        return {
            "fluxes": fluxes,
            "biomass": biomass,
            "objective_value": solution.objective_value,
            "status": solution.status,
            "target_flux": target_flux,
            "n_active_reactions": len(fluxes),
        }
    except Exception as e:
        return {
            "fluxes": {},
            "biomass": 0.0,
            "objective_value": 0.0,
            "status": "error",
            "target_flux": target_flux,
            "n_active_reactions": 0,
            "error": str(e),
        }


# ---------------------------------------------------------------------------
# Escher map generation (ported from util_legacy.generate_all_escher_maps)
# ---------------------------------------------------------------------------

def generate_escher_map(
    model: cobra.Model,
    flux_dict: dict[str, float],
    map_name: str,
    output_path: str | Path,
    *,
    title: str | None = None,
) -> Path:
    """Generate an Escher map HTML file for a given flux distribution.

    This is a thin wrapper; actual Escher rendering requires the ``escher``
    package and a map JSON file. If escher is not importable the function
    writes a minimal placeholder HTML.

    Parameters
    ----------
    model : cobra.Model
        COBRA model for reaction metadata.
    flux_dict : dict[str, float]
        Reaction-id -> flux mapping.
    map_name : str
        Map identifier (e.g. 'full', 'Core'). Used to locate a JSON file.
    output_path : str | Path
        Where to write the HTML file.
    title : str, optional
        Title shown in the HTML page header.

    Returns
    -------
    Path
        Resolved path to the written HTML file.
    """
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        import escher

        builder = escher.Builder(
            map_json=map_name if map_name.endswith(".json") else None,
            map_name=map_name if not map_name.endswith(".json") else None,
            model=model,
            reaction_data=flux_dict,
        )
        html = builder._repr_html_()
        output_path.write_text(html, encoding="utf-8")
    except ImportError:
        # Escher not available — write a summary placeholder
        _write_placeholder_map(flux_dict, output_path, title=title)

    return output_path


def _write_placeholder_map(
    flux_dict: dict[str, float],
    output_path: Path,
    *,
    title: str | None = None,
) -> None:
    """Write a minimal HTML summary when Escher is not installed."""
    title = title or "Escher Map Placeholder"
    nonzero = {k: v for k, v in flux_dict.items() if abs(v) > 1e-9}
    top10 = sorted(nonzero.items(), key=lambda kv: abs(kv[1]), reverse=True)[:10]
    rows = "".join(
        f"<tr><td>{rid}</td><td>{flux:.6f}</td></tr>" for rid, flux in top10
    )
    html = (
        f"<html><head><title>{title}</title></head><body>"
        f"<h2>{title}</h2>"
        f"<p>Total reactions with flux data: {len(flux_dict)}, "
        f"non-zero: {len(nonzero)}</p>"
        f"<table border='1'><tr><th>Reaction</th><th>Flux</th></tr>{rows}</table>"
        f"<p><em>Escher not installed — this is a placeholder.</em></p>"
        f"</body></html>"
    )
    output_path.write_text(html, encoding="utf-8")
