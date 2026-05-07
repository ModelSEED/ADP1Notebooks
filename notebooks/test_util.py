"""Tests for the new ADP1 notebooks/util.py module."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from util import (
    COMPARTMENT_MAP,
    build_gene_reaction_map,
    classify_fva_flux,
    compare_reaction_stoichiometry,
    find_significant_differences,
    get_exchange_map,
    get_reaction_directionality,
    is_diffusion_reaction,
    normalize_compartment,
    reaction_equation_with_names,
    session_for,
    standardize_exchange_id,
    translate_expression_gene_ids,
)


# ---------------------------------------------------------------------------
# session_for
# ---------------------------------------------------------------------------

class TestSessionFor:
    def test_returns_session(self):
        # Just verify it returns without error and has expected attributes
        session = session_for()
        assert hasattr(session, "cache")
        assert session._project_name == "ADP1Notebooks"

    def test_session_for_fitness_flux_fitting(self):
        session = session_for("ADP1BERDLFitnessFluxFitting.ipynb")
        assert hasattr(session, "cache")
        assert session._project_name == "ADP1Notebooks"

    def test_session_for_berdl_analysis(self):
        session = session_for("ADP1BERDLAnalysis.ipynb")
        assert hasattr(session, "cache")
        assert session._project_name == "ADP1Notebooks"

    def test_session_for_cross_sample(self):
        session = session_for("ADP1BERDLCrossSampleAnalysis.ipynb")
        assert hasattr(session, "cache")
        assert session._project_name == "ADP1Notebooks"


# ---------------------------------------------------------------------------
# Compartment utilities
# ---------------------------------------------------------------------------

class TestCompartment:
    def test_normalize_known(self):
        assert normalize_compartment("c") == "cytoplasm"
        assert normalize_compartment("e0") == "extracellular"
        assert normalize_compartment("p") == "periplasm"
        assert normalize_compartment("m0") == "mitochondria"

    def test_normalize_unknown(self):
        assert normalize_compartment("x") == "x"

    def test_map_contents(self):
        assert "c" in COMPARTMENT_MAP
        assert "e" in COMPARTMENT_MAP


# ---------------------------------------------------------------------------
# Reaction directionality
# ---------------------------------------------------------------------------

class TestDirectionality:
    def test_reversible(self, mini_model):
        rxn = mini_model.reactions.get_by_id("R1")
        assert get_reaction_directionality(rxn) == "reversible"

    def test_forward(self, mini_model):
        rxn = mini_model.reactions.get_by_id("R2")
        assert get_reaction_directionality(rxn) == "forward"

    def test_reverse(self, mini_model):
        rxn = mini_model.reactions.get_by_id("R_reverse")
        assert get_reaction_directionality(rxn) == "reverse"

    def test_blocked(self, mini_model):
        rxn = mini_model.reactions.get_by_id("R_blocked")
        assert get_reaction_directionality(rxn) == "blocked"


# ---------------------------------------------------------------------------
# Exchange ID standardization
# ---------------------------------------------------------------------------

class TestExchangeId:
    def test_already_standard(self):
        assert standardize_exchange_id("EX_glc__D_e0") == "EX_glc__D_e0"

    def test_without_compartment_suffix(self):
        assert standardize_exchange_id("EX_glc__D") == "EX_glc__D_e0"

    def test_with_e_suffix(self):
        assert standardize_exchange_id("EX_glc__D_e") == "EX_glc__D_e0"

    def test_non_exchange(self):
        assert standardize_exchange_id("R_PFK") == "R_PFK"


# ---------------------------------------------------------------------------
# Exchange map
# ---------------------------------------------------------------------------

class TestExchangeMap:
    def test_builds_map(self, mini_model):
        emap = get_exchange_map(mini_model)
        assert "EX_a_e0" in emap
        assert emap["EX_a_e0"].id == "EX_a_e0"


# ---------------------------------------------------------------------------
# Gene-reaction map
# ---------------------------------------------------------------------------

class TestGeneReactionMap:
    def test_map_contents(self, mini_model):
        grm = build_gene_reaction_map(mini_model)
        assert "gene1" in grm
        assert "gene2" in grm
        assert "R1" in grm["gene1"]
        assert "R2" in grm["gene1"]
        assert "R1" in grm["gene2"]
        assert "R2" not in grm["gene2"]


# ---------------------------------------------------------------------------
# Reaction equation with names
# ---------------------------------------------------------------------------

class TestReactionEquation:
    def test_reversible(self, mini_model):
        rxn = mini_model.reactions.get_by_id("R1")
        eq = reaction_equation_with_names(rxn)
        assert "<=>" in eq
        assert "A" in eq
        assert "B" in eq
        assert "2 C" in eq

    def test_forward_only(self, mini_model):
        rxn = mini_model.reactions.get_by_id("R2")
        eq = reaction_equation_with_names(rxn)
        assert "=>" in eq


# ---------------------------------------------------------------------------
# Diffusion detection
# ---------------------------------------------------------------------------

class TestDiffusion:
    def test_is_diffusion(self, mini_model):
        rxn = mini_model.reactions.get_by_id("diffusion_a")
        assert is_diffusion_reaction(rxn) is True

    def test_not_diffusion_multiple_mets(self, mini_model):
        rxn = mini_model.reactions.get_by_id("R1")
        assert is_diffusion_reaction(rxn) is False

    def test_not_diffusion_exchange(self, mini_model):
        rxn = mini_model.reactions.get_by_id("EX_a_e0")
        assert is_diffusion_reaction(rxn) is False


# ---------------------------------------------------------------------------
# Stoichiometry comparison
# ---------------------------------------------------------------------------

class TestStoichiometryCompare:
    def test_identical(self, mini_model):
        rxn = mini_model.reactions.get_by_id("R2")
        result = compare_reaction_stoichiometry(rxn, rxn)
        assert result["identical"] is True
        assert result["only_in_a"] == []
        assert result["only_in_b"] == []
        assert result["coefficient_diffs"] == {}

    def test_different(self, mini_model):
        rxn_a = mini_model.reactions.get_by_id("R1")
        rxn_b = mini_model.reactions.get_by_id("R2")
        result = compare_reaction_stoichiometry(rxn_a, rxn_b)
        assert result["identical"] is False


# ---------------------------------------------------------------------------
# FVA helpers
# ---------------------------------------------------------------------------

class TestFVA:
    def test_find_significant_differences(self):
        df = pd.DataFrame(
            {"minimum": [0.0, -5.0, 0.0], "maximum": [0.0, 10.0, 0.0]},
            index=["rxn1", "rxn2", "rxn3"],
        )
        result = find_significant_differences(df)
        assert len(result) == 1
        assert "rxn2" in result.index

    def test_threshold(self):
        df = pd.DataFrame(
            {"minimum": [0.0, -1e-8], "maximum": [1e-8, 1e-8]},
            index=["rxn1", "rxn2"],
        )
        result = find_significant_differences(df, threshold=1e-6)
        assert len(result) == 0


class TestClassifyFvaFlux:
    def test_blocked(self):
        assert classify_fva_flux(0.0, 0.0) == "blocked"
        assert classify_fva_flux(1e-10, -1e-10) == "blocked"

    def test_fixed_nonzero(self):
        assert classify_fva_flux(5.0, 5.0) == "fixed_nonzero"
        assert classify_fva_flux(-3.0, -3.0) == "fixed_nonzero"

    def test_variable(self):
        assert classify_fva_flux(-5.0, 10.0) == "variable"
        assert classify_fva_flux(0.0, 5.0) == "variable"
        assert classify_fva_flux(-10.0, 0.0) == "variable"

    def test_edge_near_zero(self):
        # One bound is zero-ish, other is nonzero
        assert classify_fva_flux(0.0, 1.0) == "variable"


# ---------------------------------------------------------------------------
# translate_expression_gene_ids
# ---------------------------------------------------------------------------

class TestTranslateExpressionGeneIds:
    def test_basic_translation(self, tmp_path):
        # Create a synthetic mapping CSV
        mapping = pd.DataFrame({
            "ID": ["ACIAD_RS00010", "ACIAD_RS00020", "ACIAD_RS00030"],
            "Old ID": ["ACIAD0001", "ACIAD0002", "ACIAD0003"],
            "Name": ["geneA", "geneB", "geneC"],
        })
        mapping_file = tmp_path / "mapping.csv"
        mapping.to_csv(mapping_file, index=False)

        # Create expression DataFrame with RefSeq IDs as index
        df = pd.DataFrame(
            {"cond1": [1.0, 2.0, 3.0], "cond2": [4.0, 5.0, 6.0]},
            index=["ACIAD_RS00010", "ACIAD_RS00020", "ACIAD_RS00030"],
        )

        result = translate_expression_gene_ids(df, gene_mapping_file=str(mapping_file))
        assert list(result.index) == ["ACIAD0001", "ACIAD0002", "ACIAD0003"]
        assert result.loc["ACIAD0001", "cond1"] == 1.0

    def test_unmapped_ids_kept(self, tmp_path):
        mapping = pd.DataFrame({
            "ID": ["ACIAD_RS00010"],
            "Old ID": ["ACIAD0001"],
            "Name": ["geneA"],
        })
        mapping_file = tmp_path / "mapping.csv"
        mapping.to_csv(mapping_file, index=False)

        df = pd.DataFrame(
            {"cond1": [1.0, 2.0]},
            index=["ACIAD_RS00010", "DgoA"],
        )

        result = translate_expression_gene_ids(df, gene_mapping_file=str(mapping_file))
        assert "ACIAD0001" in result.index
        assert "DgoA" in result.index

    def test_name_mapping_fallback(self, tmp_path):
        mapping = pd.DataFrame({
            "ID": ["ACIAD_RS00010"],
            "Old ID": ["ACIAD0001"],
            "Name": ["SpecialGene"],
        })
        mapping_file = tmp_path / "mapping.csv"
        mapping.to_csv(mapping_file, index=False)

        df = pd.DataFrame({"cond1": [1.0]}, index=["SpecialGene"])
        result = translate_expression_gene_ids(df, gene_mapping_file=str(mapping_file))
        assert "ACIAD0001" in result.index

    def test_duplicates_dropped(self, tmp_path):
        mapping = pd.DataFrame({
            "ID": ["ACIAD_RS00010", "ACIAD_RS00020"],
            "Old ID": ["SAME_ID", "SAME_ID"],
            "Name": ["", ""],
        })
        mapping_file = tmp_path / "mapping.csv"
        mapping.to_csv(mapping_file, index=False)

        df = pd.DataFrame(
            {"cond1": [10.0, 20.0]},
            index=["ACIAD_RS00010", "ACIAD_RS00020"],
        )

        result = translate_expression_gene_ids(df, gene_mapping_file=str(mapping_file))
        assert len(result) == 1
        assert result.loc["SAME_ID", "cond1"] == 10.0  # first kept
