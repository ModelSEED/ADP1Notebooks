"""Thin compatibility shim for notebooks that still import NotebookUtil from util_legacy.

This replaces the original 1132-line god-class with a minimal wrapper that
inherits from KBUtilLib's mixin classes (the same hierarchy the old class used)
and preserves the custom helper methods notebooks depend on.

Phase 4d will eventually remove all `_legacy` shim cells from notebooks, at
which point this file can be deleted entirely.
"""
import sys
import os

# === Path bootstrap (same as the old file) ===================================
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
base_dir = os.path.dirname(os.path.dirname(script_dir))
folder_name = os.path.basename(script_dir)

sys.path = [
    base_dir + "/KBUtilLib/src",
    base_dir + "/cobrakbase",
    base_dir + "/ModelSEEDpy/",
    "/home/chenry/MyEnvs/modelseed_cplex",
] + sys.path
# =============================================================================

from kbutillib import (
    MSFBAUtils,
    AICurationUtils,
    NotebookUtils,
    EscherUtils,
    KBPLMUtils,
    KBModelUtils,
    MSBiochemUtils,
    KBAnnotationUtils,
)

import pandas as pd
import numpy as np
from cobra import Reaction, Metabolite
from cobra.flux_analysis import pfba


class NotebookUtil(MSFBAUtils, AICurationUtils, NotebookUtils, KBPLMUtils, EscherUtils):
    """Legacy-compatible multi-inheritance util class.

    Provides the same public surface as the original util_legacy.NotebookUtil so
    that notebooks with ``from util_legacy import NotebookUtil`` continue to work
    without modification.
    """

    # Map from published model compartment names to ModelSEED compartment IDs
    COMPARTMENT_MAP = {
        "Cytosol": "c0",
        "Extraorganism": "e0",
        "Periplasm": "p0",
        "e": "e0",
    }

    def __init__(self, **kwargs):
        super().__init__(
            notebook_folder=script_dir,
            name="ADP1NotebookUtils",
            user="chenry",
            retries=5,
            proxy_port=None,
            **kwargs,
        )

    # ------------------------------------------------------------------
    # Custom helpers that are NOT in KBUtilLib — kept here verbatim
    # ------------------------------------------------------------------

    def add_dgoa_reaction(self, model):
        """Add DGOA reaction to the model."""
        rxn = Reaction("DgoA")
        rxn.name = "DgoA"
        rxn.lower_bound = 0
        rxn.upper_bound = 1000
        rxn.gene_reaction_rule = "DgoA"
        cpd00236 = model.find_met("cpd00236")[0]
        cpd00020 = model.find_met("cpd00020")[0]
        cpd02857 = model.find_met("cpd02857")[0]
        rxn.add_metabolites({cpd00236: -1, cpd00020: -1, cpd02857: 1})
        model.model.add_reactions([rxn])

    def reaction_equation_with_names(self, rxn):
        """Format reaction equation with metabolite names."""
        lhs, rhs = [], []
        for met, coeff in rxn.metabolites.items():
            name = met.name
            if coeff < 0:
                lhs.append(f"{-coeff:g} {name}" if abs(coeff) != 1 else name)
            elif coeff > 0:
                rhs.append(f"{coeff:g} {name}" if abs(coeff) != 1 else name)
        return " + ".join(lhs) + " <=> " + " + ".join(rhs)

    def get_model_and_simulate(self, model_id, media_id):
        """Load a model and simulate it with the specified media."""
        model = self.msrecon.get_model(model_id)
        if "cpd00236_c0" in model.model.metabolites:
            self.add_dgoa_reaction(model)
            model.model.reactions.get_by_id("rxn01332_c0").lower_bound = 0.0
            model.model.reactions.get_by_id("rxn01332_c0").upper_bound = 0
        media = self.msrecon.get_media(media_id)
        self._set_model_media(model, media)
        solution = pfba(model.model)
        rxn_flux_data = {}
        for rxn in model.model.reactions:
            if abs(solution.fluxes[rxn.id]) > 1e-6:
                oldupper = rxn.upper_bound
                oldlower = rxn.lower_bound
                rxn.upper_bound = 0
                rxn.lower_bound = 0
                kosol = pfba(model.model)
                rxn.upper_bound = oldupper
                rxn.lower_bound = oldlower
                rxn_flux_data[rxn.id] = {
                    "flux": solution.fluxes[rxn.id],
                    "ko_ratio": kosol.fluxes["bio1"] / solution.fluxes["bio1"],
                }
        return {"model": model, "rxn_flux_data": rxn_flux_data, "solution": solution}

    def add_model_data_to_annotations(self, model, annotations, rxn_flux_data, label):
        """Add model flux data to gene annotations dict."""
        for rxn in model.model.reactions:
            rxnstr = rxn.name + ":" + self.reaction_equation_with_names(rxn)
            for gene in rxn.genes:
                gene = str(gene)
                if gene.startswith("mRNA_"):
                    continue
                if gene not in annotations:
                    annotations[gene] = {}
                if label not in annotations[gene]:
                    annotations[gene][label] = {}
                annotations[gene][label][rxn.id] = rxnstr
                if label + "_flux" not in annotations[gene]:
                    annotations[gene][label + "_flux"] = {}
                if rxn.id in rxn_flux_data:
                    annotations[gene][label + "_flux"][rxn.id] = (
                        rxn.id + ":" + str(rxn_flux_data[rxn.id]["flux"]) + ";" + str(rxn_flux_data[rxn.id]["ko_ratio"])
                    )
                else:
                    annotations[gene][label + "_flux"][rxn.id] = rxn.id + ":0;1"
        return annotations

    def find_significant_differences(self, flux1, flux2, threshold=0.01):
        """Find reactions with significantly different fluxes."""
        sig_reactions = []
        for rxn_id in flux1.index:
            if abs(flux1[rxn_id] - flux2[rxn_id]) > threshold:
                sig_reactions.append({
                    "Reaction": rxn_id,
                    "Flux1": flux1[rxn_id],
                    "Flux2": flux2[rxn_id],
                    "Difference": flux2[rxn_id] - flux1[rxn_id],
                })
        return sig_reactions

    def reactions_to_genes(self, reaction_list, model):
        """Convert reaction list to gene list."""
        genes = set()
        reaction_gene_map = {}
        for rxn_data in reaction_list:
            rxn_id = rxn_data["Reaction"]
            try:
                rxn = model.model.reactions.get_by_id(rxn_id)
                rxn_genes = [str(gene) for gene in rxn.genes if not str(gene).startswith("mRNA_")]
                reaction_gene_map[rxn_id] = rxn_genes
                genes.update(rxn_genes)
            except Exception:
                reaction_gene_map[rxn_id] = []
        return list(genes), reaction_gene_map

    def translate_expression_gene_ids(self, expression, gene_mapping_file="data/ADP1Genes.csv"):
        """Translate expression gene IDs from RefSeq format to old format."""
        gene_mapping_df = pd.read_csv(gene_mapping_file)
        refseq_to_old = {}
        name_to_old = {}
        for _, row in gene_mapping_df.iterrows():
            refseq_id = row["ID"]
            old_id = row["Old ID"]
            name = row["Name"]
            if pd.notna(refseq_id) and pd.notna(old_id):
                refseq_to_old[refseq_id] = old_id
            if pd.notna(name) and pd.notna(old_id) and name != "":
                name_to_old[name] = old_id

        df = expression.get_dataframe()
        original_gene_ids = df.index.tolist()
        translated_gene_ids = []
        for gene_id in original_gene_ids:
            gene_id_str = str(gene_id)
            if gene_id_str in refseq_to_old:
                translated_gene_ids.append(refseq_to_old[gene_id_str])
            elif gene_id_str in name_to_old:
                translated_gene_ids.append(name_to_old[gene_id_str])
            else:
                translated_gene_ids.append(gene_id_str)

        translated_df = df.copy()
        translated_df.index = translated_gene_ids
        translated_df = translated_df[~translated_df.index.duplicated(keep="first")]

        import copy
        from cobra.core import DictList

        translated_expression = copy.deepcopy(expression)
        translated_expression._data = translated_df
        translated_expression.features = DictList()

        class ExpressionFeature:
            def __init__(self, feature_id, expression_obj):
                self.id = feature_id
                self._id = feature_id
                self.feature_id = feature_id
                self.expression = expression_obj

        for gene_id in translated_df.index:
            feature = ExpressionFeature(gene_id, translated_expression)
            translated_expression.features.append(feature)

        class SimpleGenome:
            def __init__(self, gene_ids):
                self.gene_ids = set(gene_ids)

            def search_for_gene(self, gene_id):
                if gene_id in self.gene_ids:
                    class GeneFeature:
                        def __init__(self, gid):
                            self.id = gid
                            self._id = gid
                    return GeneFeature(gene_id)
                return None

        translated_expression.object = SimpleGenome(translated_df.index.tolist())
        self.logger.info(f"  Translated {len(original_gene_ids)} -> {len(translated_df)} unique gene IDs")
        self.logger.info(f"  Updated {len(translated_expression.features)} features in MSExpression")
        return translated_expression

    def average_expression_replicates(self, expression, strain_list):
        """Average expression replicates for each strain."""
        expression_df = expression._data.copy()
        averaged_data = {"index": expression_df.index}
        for strain in strain_list:
            replicate_cols = [col for col in expression_df.columns if col.startswith(f"{strain}_")]
            if replicate_cols:
                averaged_data[strain] = expression_df[replicate_cols].mean(axis=1)
                self.logger.info(f"Averaged {len(replicate_cols)} replicates for strain {strain}")
            elif strain in expression_df.columns:
                averaged_data[strain] = expression_df[strain]
                self.logger.info(f"No replicates found for {strain}, using existing column")
            else:
                self.logger.warning(f"No data found for strain {strain}")

        averaged_df = pd.DataFrame(averaged_data)
        averaged_df.set_index("index", inplace=True)

        import copy
        from cobra.core import DictList

        averaged_expression = copy.deepcopy(expression)
        averaged_expression._data = averaged_df

        class ExpressionCondition:
            def __init__(self, condition_id):
                self.id = condition_id
                self._id = condition_id

        averaged_expression.conditions = DictList()
        for strain in strain_list:
            if strain in averaged_df.columns:
                condition = ExpressionCondition(strain)
                averaged_expression.conditions.append(condition)

        self.logger.info(f"Created averaged expression data with {len(averaged_expression.conditions)} conditions")
        return averaged_expression

    def process_strain_with_expression(self, strain, expression_data, base_model, media, with_dgoa,
                                       knockout_dahp=True, biomass_fraction=0.25,
                                       default_coef=0.01, activation_threshold=None,
                                       deactivation_threshold=0.000001, minimal_flux=0.001):
        """Process a single strain condition with expression constraints."""
        try:
            import cobra.io
            from modelseedpy import MSModelUtil

            model_copy = MSModelUtil.from_cobrapy(cobra.io.json.to_json(base_model.model))
            model_copy.util = self

            if with_dgoa:
                self.add_dgoa_reaction(model_copy)

            if knockout_dahp and "rxn01332_c0" in [rxn.id for rxn in model_copy.model.reactions]:
                model_copy.model.reactions.get_by_id("rxn01332_c0").lower_bound = 0
                model_copy.model.reactions.get_by_id("rxn01332_c0").upper_bound = 0

            self._set_model_media(model_copy, media)
            self.constrain_objective_to_fraction_of_optimum(
                model_copy, objective="GROWTH_DASH_RXN", fraction=biomass_fraction
            )

            analysis_result = expression_data.fit_model_flux_to_data(
                model=model_copy, condition=strain, default_coef=default_coef,
                activation_threshold=activation_threshold,
                deactivation_threshold=deactivation_threshold,
            )

            off_off_list = analysis_result.get("off_off", [])
            on_on_list = analysis_result.get("on_on", [])

            for rxn_id in off_off_list:
                try:
                    rxn = model_copy.model.reactions.get_by_id(rxn_id)
                    rxn.lower_bound = 0
                    rxn.upper_bound = 0
                except Exception:
                    pass

            for rxn_id in on_on_list:
                try:
                    rxn = model_copy.model.reactions.get_by_id(rxn_id)
                    temp_sol = pfba(model_copy.model)
                    flux_direction = temp_sol.fluxes.get(rxn_id, 0)
                    if flux_direction > 0:
                        rxn.lower_bound = max(rxn.lower_bound, minimal_flux)
                    elif flux_direction < 0:
                        rxn.upper_bound = min(rxn.upper_bound, -minimal_flux)
                    else:
                        if rxn.upper_bound >= minimal_flux:
                            rxn.lower_bound = max(rxn.lower_bound, minimal_flux)
                except Exception:
                    pass

            solution = pfba(model_copy.model)
            fluxes = {rxn.id: solution.fluxes[rxn.id] for rxn in model_copy.model.reactions if abs(solution.fluxes[rxn.id]) > 1e-9}
            dgoa_flux = solution.fluxes.get("DgoA") if with_dgoa else None
            biomass = solution.fluxes.get("GROWTH_DASH_RXN", 0)

            return {
                "fluxes": fluxes,
                "biomass": biomass,
                "active_reactions": len(fluxes),
                "off_reactions": off_off_list,
                "on_reactions": on_on_list,
                "dgoa_flux": dgoa_flux,
                "solution_status": solution.status,
            }
        except Exception as e:
            self.logger.error(f"Error processing strain {strain} (with_dgoa={with_dgoa}): {str(e)}")
            return {
                "fluxes": {},
                "biomass": 0,
                "active_reactions": 0,
                "off_reactions": [],
                "on_reactions": [],
                "dgoa_flux": None,
                "solution_status": "error",
                "error": str(e),
            }

    def validate_expression_flux_solution(self, solution_dict, strain, dgoa_status):
        """Validate expression-constrained flux solution."""
        validation_messages = []
        passed = True

        if "error" in solution_dict:
            return False, f"FAIL {strain}_{dgoa_status}: Solution error - {solution_dict['error']}"

        biomass = solution_dict.get("biomass", 0)
        if biomass <= 0:
            validation_messages.append(f"Biomass={biomass:.6f} (FAIL)")
            passed = False
        else:
            validation_messages.append(f"Biomass={biomass:.4f} OK")

        status = solution_dict.get("solution_status", "unknown")
        if status != "optimal":
            validation_messages.append(f"Status={status} (FAIL)")
            passed = False
        else:
            validation_messages.append(f"Status={status} OK")

        active_count = solution_dict.get("active_reactions", 0)
        if active_count < 50:
            validation_messages.append(f"Active_Rxns={active_count} (FAIL)")
            passed = False
        else:
            validation_messages.append(f"Active_Rxns={active_count} OK")

        if dgoa_status == "with_dgoa":
            dgoa_flux = solution_dict.get("dgoa_flux")
            if dgoa_flux is None or dgoa_flux <= 0:
                validation_messages.append(f"DGOA_Flux={dgoa_flux} (FAIL)")
                passed = False
            else:
                validation_messages.append(f"DGOA_Flux={dgoa_flux:.4f} OK")
        else:
            validation_messages.append("DGOA_Flux=N/A")

        symbol = "PASS" if passed else "FAIL"
        message = f"{symbol} {strain}_{dgoa_status}: " + ", ".join(validation_messages)
        return passed, message

    def create_expression_flux_summary(self, results_dict):
        """Create summary dictionary from expression flux analysis results."""
        summary = {}
        for condition_key, result_data in results_dict.items():
            if "_with_dgoa" in condition_key:
                strain = condition_key.replace("_with_dgoa", "")
                dgoa_status = "with_dgoa"
            elif "_without_dgoa" in condition_key:
                strain = condition_key.replace("_without_dgoa", "")
                dgoa_status = "without_dgoa"
            else:
                continue
            if strain not in summary:
                summary[strain] = {}
            condition_summary = {
                "biomass": result_data.get("biomass", 0),
                "active_reactions": result_data.get("active_reactions", 0),
                "off_reactions": result_data.get("off_reactions", []),
                "on_reactions": result_data.get("on_reactions", []),
                "solution_status": result_data.get("solution_status", "unknown"),
                "dgoa_flux": result_data.get("dgoa_flux", 0),
            }
            if "error" in result_data:
                condition_summary["error"] = result_data["error"]
            summary[strain][dgoa_status] = condition_summary
        self.logger.info(f"Created summary for {len(summary)} strains")
        return summary

    def classify_fva_flux(self, fva_result, flux_value, zero_tol=1e-5):
        """Classify a reaction based on FVA min/max and pFBA flux."""
        mn = fva_result.get("MIN", 0) or 0
        mx = fva_result.get("MAX", 0) or 0
        if abs(mn) < zero_tol and abs(mx) < zero_tol:
            return "blocked"
        if mn > zero_tol:
            return "essential_fwd"
        if mx < -zero_tol:
            return "essential_rev"
        if abs(mn) < zero_tol and mx > zero_tol:
            return "optionally_active_fwd" if abs(flux_value) > zero_tol else "optionally_active_zero_fwd"
        if mn < -zero_tol and abs(mx) < zero_tol:
            return "optionally_active_rev" if abs(flux_value) > zero_tol else "optionally_active_zero_rev"
        if flux_value > zero_tol:
            return "variable_fwd"
        if flux_value < -zero_tol:
            return "variable_rev"
        return "variable_zero"

    def standardize_exchange_id(self, rxn):
        """Return standardized exchange ID based on the metabolite being exchanged."""
        mets = list(rxn.metabolites.keys())
        if len(mets) == 1:
            met_id = mets[0].id
            base_id = met_id.rsplit("_", 1)[0] if "_" in met_id else met_id
            return f"EX_{base_id}"
        return rxn.id

    def get_exchange_map(self, model):
        """Build map from standardized exchange ID to original reaction ID."""
        ex_map = {}
        for rxn in model.reactions:
            mets = list(rxn.metabolites.keys())
            if len(mets) == 1:
                std_id = self.standardize_exchange_id(rxn)
                ex_map[std_id] = rxn.id
        return ex_map

    def compare_reaction_stoichiometry(self, rxn1, rxn2):
        """Compare stoichiometry between two reactions."""
        diffs = []
        mets1 = {m.id: c for m, c in rxn1.metabolites.items()}
        mets2 = {m.id: c for m, c in rxn2.metabolites.items()}
        for met_id in sorted(set(mets1.keys()) | set(mets2.keys())):
            c1 = mets1.get(met_id)
            c2 = mets2.get(met_id)
            if c1 is None:
                diffs.append(f"{met_id}: missing in published, {c2:g} in modelseed")
            elif c2 is None:
                diffs.append(f"{met_id}: {c1:g} in published, missing in modelseed")
            elif abs(c1 - c2) > 1e-9:
                diffs.append(f"{met_id}: {c1:g} in published vs {c2:g} in modelseed")
        return diffs

    def get_reaction_directionality(self, rxn):
        """Return directionality string for a reaction based on bounds."""
        if rxn.lower_bound < 0 and rxn.upper_bound > 0:
            return "reversible"
        if rxn.lower_bound >= 0:
            return "forward"
        if rxn.upper_bound <= 0:
            return "reverse"
        return "unknown"

    def normalize_compartment(self, compartment):
        """Normalize compartment name to ModelSEED convention."""
        return self.COMPARTMENT_MAP.get(compartment, compartment)

    def is_diffusion_reaction(self, rxn):
        """Check if a reaction is a pure single-compound diffusion reaction."""
        name_stoich = {}
        for met, coeff in rxn.metabolites.items():
            name_stoich[met.name] = name_stoich.get(met.name, 0) + coeff
        if not all(abs(v) < 1e-9 for v in name_stoich.values()):
            return False
        return len(name_stoich) == 1

    def build_gene_reaction_map(self, model):
        """Build map from gene ID to list of reaction IDs."""
        gene_rxn_map = {}
        for rxn in model.reactions:
            for gene in rxn.genes:
                gid = str(gene)
                if gid.startswith("mRNA_"):
                    continue
                if gid not in gene_rxn_map:
                    gene_rxn_map[gid] = []
                gene_rxn_map[gid].append(rxn.id)
        return gene_rxn_map


# Module-level instance for backwards compatibility (some cells do `util = NotebookUtil()`)
util = NotebookUtil()
