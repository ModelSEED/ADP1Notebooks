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
    generate_escher_map,
    get_exchange_map,
    get_reaction_directionality,
    is_diffusion_reaction,
    normalize_compartment,
    reaction_equation_with_names,
    run_fold_change_simulation,
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


# ---------------------------------------------------------------------------
# run_fold_change_simulation
# ---------------------------------------------------------------------------

class TestRunFoldChangeSimulation:
    @pytest.fixture
    def fba_model(self):
        """A tiny model that can actually run FBA."""
        import cobra

        model = cobra.Model("fba_test")

        # Metabolites
        a_e = cobra.Metabolite("a_e", compartment="e")
        a_c = cobra.Metabolite("a_c", compartment="c")
        b_c = cobra.Metabolite("b_c", compartment="c")
        biomass_met = cobra.Metabolite("biomass_c", compartment="c")

        # Exchange: EX_a_e0 -> a_e
        ex = cobra.Reaction("EX_a_e0")
        ex.lower_bound = -10
        ex.upper_bound = 1000
        ex.add_metabolites({a_e: -1})

        # Transport: a_e -> a_c
        transport = cobra.Reaction("transport_a")
        transport.lower_bound = 0
        transport.upper_bound = 1000
        transport.add_metabolites({a_e: -1, a_c: 1})
        transport.gene_reaction_rule = "gene1"

        # Reaction: a_c -> b_c
        r1 = cobra.Reaction("R1")
        r1.lower_bound = 0
        r1.upper_bound = 1000
        r1.add_metabolites({a_c: -1, b_c: 1})
        r1.gene_reaction_rule = "gene2"

        # Biomass: b_c -> biomass
        bio = cobra.Reaction("GROWTH_DASH_RXN")
        bio.lower_bound = 0
        bio.upper_bound = 1000
        bio.add_metabolites({b_c: -1, biomass_met: 1})

        # Demand for biomass
        demand = cobra.Reaction("DM_biomass")
        demand.lower_bound = 0
        demand.upper_bound = 1000
        demand.add_metabolites({biomass_met: -1})

        model.add_reactions([ex, transport, r1, bio, demand])
        model.objective = "GROWTH_DASH_RXN"
        return model

    def test_produces_nonzero_biomass(self, fba_model):
        reference_flux = {
            "EX_a_e0": -5.0,
            "transport_a": 5.0,
            "R1": 5.0,
            "GROWTH_DASH_RXN": 5.0,
            "DM_biomass": 5.0,
        }
        fold_changes = {"gene1": 1.5, "gene2": 0.8}

        result = run_fold_change_simulation(
            fba_model,
            fold_changes,
            reference_flux=reference_flux,
        )

        assert result["status"] == "optimal"
        assert result["biomass"] > 0
        assert result["n_active_reactions"] > 0

    def test_flux_vector_shape(self, fba_model):
        reference_flux = {
            "EX_a_e0": -5.0,
            "transport_a": 5.0,
            "R1": 5.0,
            "GROWTH_DASH_RXN": 5.0,
            "DM_biomass": 5.0,
        }
        fold_changes = {"gene1": 1.0, "gene2": 1.0}

        result = run_fold_change_simulation(
            fba_model,
            fold_changes,
            reference_flux=reference_flux,
        )

        assert isinstance(result["fluxes"], dict)
        assert isinstance(result["target_flux"], dict)
        assert result["status"] == "optimal"

    def test_empty_fold_changes(self, fba_model):
        reference_flux = {"R1": 5.0}
        result = run_fold_change_simulation(
            fba_model,
            {},
            reference_flux=reference_flux,
        )
        # Should still solve successfully
        assert result["status"] == "optimal"


# ---------------------------------------------------------------------------
# generate_escher_map
# ---------------------------------------------------------------------------

class TestGenerateEscherMap:
    def test_writes_html_file(self, mini_model, tmp_path):
        flux = {"R1": 1.5, "R2": -0.3}
        out = tmp_path / "subdir" / "map.html"
        result = generate_escher_map(
            mini_model, flux, "test_map", out, title="Test Map"
        )
        assert result == out
        assert out.exists()
        content = out.read_text()
        assert "R1" in content
        assert "<html>" in content

    def test_creates_parent_directories(self, mini_model, tmp_path):
        out = tmp_path / "a" / "b" / "c" / "map.html"
        generate_escher_map(mini_model, {"R1": 1.0}, "test", out)
        assert out.exists()

    def test_empty_flux(self, mini_model, tmp_path):
        out = tmp_path / "empty.html"
        generate_escher_map(mini_model, {}, "test", out, title="Empty")
        assert out.exists()
