import sys
import os
import json
from os import path
from zipfile import ZipFile

# Add the parent directory to the sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
script_path = os.path.abspath(__file__)
script_dir = os.path.dirname(script_path)
base_dir = os.path.dirname(os.path.dirname(script_dir))
folder_name = os.path.basename(script_dir)

print(base_dir+"/KBUtilLib/src")
sys.path = [base_dir+"/KBUtilLib/src",base_dir+"/cobrakbase/cobrakbase",base_dir+"/ModelSEEDpy/","/home/chenry/MyEnvs/modelseed_cplex"] + sys.path

# Import utilities with error handling
from kbutillib import MSFBAUtils, AICurationUtils, NotebookUtils, EscherUtils, KBPLMUtils,KBModelUtils,MSBiochemUtils,KBAnnotationUtils

import hashlib
import pandas as pd
from pandas import DataFrame, read_csv, concat, set_option
import cobra
from cobra import Reaction, Metabolite
from cobra.flux_analysis import pfba
from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem, MSExpression
import numpy as np
from scipy import stats as scipy_stats
import matplotlib.pyplot as plt
from IPython.display import HTML, display

# Define the base classes based on what's available
# Note: KBPLMUtils inherits from KBGenomeUtils, so we use KBPLMUtils instead of KBGenomeUtils
class NotebookUtil(MSFBAUtils, AICurationUtils, NotebookUtils, KBPLMUtils, EscherUtils):
    def __init__(self,**kwargs):
        super().__init__(
            notebook_folder=script_dir,
            name="ADP1NotebookUtils",
            user="chenry",
            retries=5,
            proxy_port=None,
            **kwargs
        )

    def get_model_and_simulate(self, model_id, media_id):
        """
        Load a model and simulate it with the specified media.
        """
        model = self.msrecon.get_model(model_id)
        if "cpd00236_c0" in model.model.metabolites:
            util.add_dgoa_reaction(model)
            model.model.reactions.get_by_id("rxn01332_c0").lower_bound = 0.0
            model.model.reactions.get_by_id("rxn01332_c0").upper_bound = 0
        media = self.msrecon.get_media(media_id)
        self._set_model_media(model,media)
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
    
    def add_model_data_to_annotations(self, model, annotations,rxn_flux_data,label):
        """
        Add model data to annotations.
        """
        for rxn in model.model.reactions:
            rxnstr = rxn.name+":"+self.reaction_equation_with_names(rxn)
            for gene in rxn.genes:
                gene = str(gene)
                if gene.startswith("mRNA_"):
                    continue
                if gene not in annotations:
                    annotations[gene] = {}
                if label not in annotations[gene]:
                    annotations[gene][label] = {}
                annotations[gene][label][rxn.id] = rxnstr
                if label+"_flux" not in annotations[gene]:
                    annotations[gene][label+"_flux"] = {}
                if rxn.id in rxn_flux_data:
                    annotations[gene][label+"_flux"][rxn.id] = rxn.id+":"+str(rxn_flux_data[rxn.id]["flux"])+";"+str(rxn_flux_data[rxn.id]["ko_ratio"])
                else:
                    annotations[gene][label+"_flux"][rxn.id] = rxn.id+":0;1"
        return annotations

    # Function to add DGOA reaction to the model
    def add_dgoa_reaction(self, model):
        """
        Add DGOA reaction to the model.
        """
        # Create the reaction
        rxn = Reaction('DgoA')
        rxn.name = 'DgoA'
        rxn.lower_bound = 0  # irreversible
        rxn.upper_bound = 1000
        rxn.gene_reaction_rule = 'DgoA'
        # Get metabolites from the model by their IDs
        cpd00236 = model.find_met('cpd00236')[0]
        cpd00020 = model.find_met('cpd00020')[0]
        cpd02857 = model.find_met('cpd02857')[0]

        # Add metabolites to the reaction (negative for reactants, positive for products)
        rxn.add_metabolites({
            cpd00236: -1,
            cpd00020: -1,
            cpd02857: 1
        })

        # Add the reaction to the model
        model.model.add_reactions([rxn])

    # Function to format reaction equations with metabolite names
    def reaction_equation_with_names(self,rxn):
        lhs = []
        rhs = []
        for met, coeff in rxn.metabolites.items():
            name = met.name
            if coeff < 0:
                lhs.append(f"{-coeff:g} {name}" if abs(coeff) != 1 else name)
            elif coeff > 0:
                rhs.append(f"{coeff:g} {name}" if abs(coeff) != 1 else name)
        return " + ".join(lhs) + " <=> " + " + ".join(rhs)
    
    # Function to identify significantly different fluxes
    def find_significant_differences(self,flux1, flux2, threshold=0.01):
        """Find reactions with significantly different fluxes"""
        sig_reactions = []
        for rxn_id in flux1.index:
            if abs(flux1[rxn_id] - flux2[rxn_id]) > threshold:
                sig_reactions.append({
                    'Reaction': rxn_id,
                    'Flux1': flux1[rxn_id],
                    'Flux2': flux2[rxn_id],
                    'Difference': flux2[rxn_id] - flux1[rxn_id]
                })
        return sig_reactions
    
    # Function to convert reactions to genes
    def reactions_to_genes(self,reaction_list, model):
        """Convert reaction list to gene list"""
        genes = set()
        reaction_gene_map = {}

        for rxn_data in reaction_list:
            rxn_id = rxn_data['Reaction']
            try:
                rxn = model.model.reactions.get_by_id(rxn_id)
                rxn_genes = [str(gene) for gene in rxn.genes if not str(gene).startswith("mRNA_")]
                reaction_gene_map[rxn_id] = rxn_genes
                genes.update(rxn_genes)
            except:
                reaction_gene_map[rxn_id] = []

        return list(genes), reaction_gene_map

    def process_strain_with_expression(self, strain, expression_data, base_model, media, with_dgoa,
                                      knockout_dahp=True, biomass_fraction=0.25,
                                      default_coef=0.01, activation_threshold=None,
                                      deactivation_threshold=0.000001, minimal_flux=0.001):
        """
        Process a single strain condition with expression constraints.

        Args:
            strain (str): Strain name (e.g., "ACN2586")
            expression_data (MSExpression): Expression data object
            base_model (MSModelUtil): Base metabolic model
            media (MSMedia): Media conditions
            with_dgoa (bool): Whether to add DGOA reaction
            knockout_dahp (bool): Whether to knock out native DAHP synthase (default: True)
            biomass_fraction (float): Minimum fraction of optimal biomass (default: 0.25)
            default_coef (float): Default coefficient for expression fitting (default: 0.01)
            activation_threshold (float): Threshold for activation constraints (default: None)
            deactivation_threshold (float): Threshold for deactivation constraints (default: 0.000001)
            minimal_flux (float): Minimal flux to enforce for on_on reactions (default: 0.001)

        Returns:
            dict: Dictionary containing:
                - fluxes: Dict of reaction_id -> flux_value for non-zero fluxes
                - biomass: Biomass flux value
                - active_reactions: Count of reactions with non-zero flux
                - off_reactions: List of reaction IDs turned off by expression
                - on_reactions: List of reaction IDs forced on by expression
                - dgoa_flux: DGOA reaction flux (if with_dgoa=True, else None)
                - solution_status: Optimization status string
        """
        try:
            # Create a deep copy of the model for this condition
            import cobra.io
            model_copy = MSModelUtil.from_cobrapy(cobra.io.json.to_json(base_model.model))
            model_copy.util = self
            
            # Apply DGOA modification if requested
            if with_dgoa:
                self.add_dgoa_reaction(model_copy)

            # Knock out native DAHP synthase if requested
            if knockout_dahp and "rxn01332_c0" in [rxn.id for rxn in model_copy.model.reactions]:
                model_copy.model.reactions.get_by_id("rxn01332_c0").lower_bound = 0
                model_copy.model.reactions.get_by_id("rxn01332_c0").upper_bound = 0

            # Set media conditions
            self._set_model_media(model,media)

            # Constrain biomass to specified fraction of optimal
            self.constrain_objective_to_fraction_of_optimum(
                model_copy,
                objective="GROWTH_DASH_RXN",
                fraction=biomass_fraction
            )

            # Fit model flux to expression data
            analysis_result = expression_data.fit_model_flux_to_data(
                model=model_copy,
                condition=strain,
                default_coef=default_coef,
                activation_threshold=activation_threshold,
                deactivation_threshold=deactivation_threshold
            )

            # Extract off_off list (reactions to be deactivated)
            off_off_list = analysis_result.get("off_off", [])

            # Extract on_on list (reactions to be activated)
            on_on_list = analysis_result.get("on_on", [])

            # Apply deactivation constraints
            for rxn_id in off_off_list:
                try:
                    rxn = model_copy.model.reactions.get_by_id(rxn_id)
                    rxn.lower_bound = 0
                    rxn.upper_bound = 0
                except:
                    # Reaction might not exist in model
                    pass

            # Apply activation constraints based on flux direction
            for rxn_id in on_on_list:
                try:
                    rxn = model_copy.model.reactions.get_by_id(rxn_id)
                    # Run a quick FBA to determine flux direction
                    temp_sol = pfba(model_copy.model)
                    flux_direction = temp_sol.fluxes.get(rxn_id, 0)

                    # Set minimal flux based on direction
                    if flux_direction > 0:
                        # Positive flux - set lower bound
                        rxn.lower_bound = max(rxn.lower_bound, minimal_flux)
                    elif flux_direction < 0:
                        # Negative flux - set upper bound
                        rxn.upper_bound = min(rxn.upper_bound, -minimal_flux)
                    else:
                        # Zero flux - try setting positive minimal flux
                        if rxn.upper_bound >= minimal_flux:
                            rxn.lower_bound = max(rxn.lower_bound, minimal_flux)
                except:
                    # Reaction might not exist in model
                    pass

            # Run final pFBA with expression constraints applied
            solution = pfba(model_copy.model)

            # Extract fluxes (only non-zero)
            fluxes = {}
            for rxn in model_copy.model.reactions:
                if abs(solution.fluxes[rxn.id]) > 1e-9:
                    fluxes[rxn.id] = solution.fluxes[rxn.id]

            # Get DGOA flux if applicable
            dgoa_flux = None
            if with_dgoa and "DgoA" in solution.fluxes:
                dgoa_flux = solution.fluxes["DgoA"]

            # Get biomass flux
            biomass = solution.fluxes.get("GROWTH_DASH_RXN", 0)

            return {
                "fluxes": fluxes,
                "biomass": biomass,
                "active_reactions": len(fluxes),
                "off_reactions": off_off_list,
                "on_reactions": on_on_list,
                "dgoa_flux": dgoa_flux,
                "solution_status": solution.status
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
                "error": str(e)
            }

    def translate_expression_gene_ids(self, expression, gene_mapping_file="data/ADP1Genes.csv"):
        """
        Translate expression gene IDs from RefSeq format (ACIAD_RSxxxxx) to old format (ACIADxxxx).

        Args:
            expression (MSExpression): Expression data with RefSeq gene IDs
            gene_mapping_file (str): Path to CSV file with ID mapping (default: "data/ADP1Genes.csv")

        Returns:
            MSExpression: New expression object with translated gene IDs
        """
        import pandas as pd

        # Load gene mapping
        gene_mapping_df = pd.read_csv(gene_mapping_file)

        # Create mapping dictionaries
        refseq_to_old = {}
        name_to_old = {}

        for _, row in gene_mapping_df.iterrows():
            refseq_id = row['ID']
            old_id = row['Old ID']
            name = row['Name']

            if pd.notna(refseq_id) and pd.notna(old_id):
                refseq_to_old[refseq_id] = old_id

            if pd.notna(name) and pd.notna(old_id) and name != '':
                name_to_old[name] = old_id

        # Get expression DataFrame
        df = expression.get_dataframe()

        # Translate gene IDs
        original_gene_ids = df.index.tolist()
        translated_gene_ids = []

        for gene_id in original_gene_ids:
            gene_id_str = str(gene_id)

            # Try RefSeq mapping first
            if gene_id_str in refseq_to_old:
                translated_gene_ids.append(refseq_to_old[gene_id_str])
            # Try name mapping (for special genes like DgoA)
            elif gene_id_str in name_to_old:
                translated_gene_ids.append(name_to_old[gene_id_str])
            # Keep as is (might already be in correct format, or introduced genes like DgoA, DgoA_Ec, Kan)
            else:
                translated_gene_ids.append(gene_id_str)

        # Create new DataFrame with translated IDs
        translated_df = df.copy()
        translated_df.index = translated_gene_ids

        # Remove duplicates by keeping first occurrence
        translated_df = translated_df[~translated_df.index.duplicated(keep='first')]

        # Create a deep copy of the expression object and replace its data
        import copy
        from cobra.core import DictList
        translated_expression = copy.deepcopy(expression)
        translated_expression._data = translated_df

        # CRITICAL: Update the features DictList to match translated gene IDs
        # The features DictList is used by build_reaction_expression to look up genes
        translated_expression.features = DictList()

        # Recreate feature objects with translated IDs
        for gene_id in translated_df.index:
            # Create a simple feature object
            class ExpressionFeature:
                def __init__(self, feature_id, expression_obj):
                    self.id = feature_id
                    self._id = feature_id
                    self.feature_id = feature_id
                    self.expression = expression_obj

            feature = ExpressionFeature(gene_id, translated_expression)
            translated_expression.features.append(feature)

        # Create a simple object that implements search_for_gene for build_reaction_expression
        class SimpleGenome:
            def __init__(self, gene_ids):
                self.gene_ids = set(gene_ids)

            def search_for_gene(self, gene_id):
                # Return a simple object with an id attribute if gene exists
                if gene_id in self.gene_ids:
                    class GeneFeature:
                        def __init__(self, gid):
                            self.id = gid
                            self._id = gid
                    return GeneFeature(gene_id)
                return None

        # Set the object attribute to our simple genome
        translated_expression.object = SimpleGenome(translated_df.index.tolist())

        self.logger.info(f"  Translated {len(original_gene_ids)} -> {len(translated_df)} unique gene IDs")
        self.logger.info(f"  Updated {len(translated_expression.features)} features in MSExpression")

        return translated_expression

    def average_expression_replicates(self, expression, strain_list):
        """
        Average expression replicates for each strain.

        Takes an MSExpression object with replicate columns (e.g., ACN2586_1, ACN2586_2, ...)
        and averages them to create single columns per strain (e.g., ACN2586).

        Args:
            expression (MSExpression): Expression data object with replicates
            strain_list (list): List of strain names (e.g., ["ACN2586", "ACN2821", ...])

        Returns:
            MSExpression: New MSExpression object with averaged data per strain
        """
        try:
            # Access the underlying DataFrame
            expression_df = expression._data.copy()

            # Create new DataFrame for averaged data
            averaged_data = {}

            # Keep the index (gene/protein IDs)
            averaged_data['index'] = expression_df.index

            # For each strain, find and average its replicates
            for strain in strain_list:
                # Find columns that match this strain pattern (e.g., ACN2586_1, ACN2586_2, ...)
                replicate_cols = [col for col in expression_df.columns if col.startswith(f"{strain}_")]

                if replicate_cols:
                    # Average the replicates
                    averaged_data[strain] = expression_df[replicate_cols].mean(axis=1)
                    self.logger.info(f"Averaged {len(replicate_cols)} replicates for strain {strain}")
                else:
                    # No replicates found - check if strain column exists as-is
                    if strain in expression_df.columns:
                        averaged_data[strain] = expression_df[strain]
                        self.logger.info(f"No replicates found for {strain}, using existing column")
                    else:
                        self.logger.warning(f"No data found for strain {strain}")

            # Create new DataFrame from averaged data
            averaged_df = pd.DataFrame(averaged_data)
            averaged_df.set_index('index', inplace=True)

            # Create a new MSExpression object with averaged data
            # We need to save to a temporary file and reload, or create from scratch
            # The simplest approach is to modify the _data attribute of a copy

            # Create a deep copy of the expression object
            import copy
            averaged_expression = copy.deepcopy(expression)

            # Replace the data with averaged data
            averaged_expression._data = averaged_df

            # Update conditions list to match new columns
            # DictList is from cobra.core, we need to create condition objects properly
            from cobra.core import DictList

            # Create a simple class to represent expression conditions
            class ExpressionCondition:
                def __init__(self, condition_id):
                    self.id = condition_id
                    self._id = condition_id

            # Clear and rebuild conditions
            averaged_expression.conditions = DictList()
            for strain in strain_list:
                if strain in averaged_df.columns:
                    condition = ExpressionCondition(strain)
                    averaged_expression.conditions.append(condition)

            self.logger.info(f"Created averaged expression data with {len(averaged_expression.conditions)} conditions")

            return averaged_expression

        except Exception as e:
            self.logger.error(f"Error averaging expression replicates: {str(e)}")
            raise

    def validate_expression_flux_solution(self, solution_dict, strain, dgoa_status):
        """
        Validate expression-constrained flux solution.

        Checks that the solution meets biological feasibility criteria:
        - Biomass production is positive
        - Optimization status is optimal
        - Sufficient reactions are active
        - DGOA flux is positive when DGOA is present

        Args:
            solution_dict (dict): Solution dictionary from process_strain_with_expression()
            strain (str): Strain name for reporting
            dgoa_status (str): Either "with_dgoa" or "without_dgoa"

        Returns:
            tuple: (pass/fail bool, validation message string)
        """
        validation_messages = []
        passed = True

        # Check if solution has an error
        if "error" in solution_dict:
            return False, f"❌ {strain}_{dgoa_status}: Solution error - {solution_dict['error']}"

        # Validate biomass production
        biomass = solution_dict.get("biomass", 0)
        if biomass <= 0:
            validation_messages.append(f"Biomass={biomass:.6f} (FAIL: ≤0)")
            passed = False
        else:
            validation_messages.append(f"Biomass={biomass:.4f} ✓")

        # Validate optimization status
        status = solution_dict.get("solution_status", "unknown")
        if status != "optimal":
            validation_messages.append(f"Status={status} (FAIL: not optimal)")
            passed = False
        else:
            validation_messages.append(f"Status={status} ✓")

        # Validate active reactions count
        active_count = solution_dict.get("active_reactions", 0)
        if active_count < 50:
            validation_messages.append(f"Active_Rxns={active_count} (FAIL: <50)")
            passed = False
        else:
            validation_messages.append(f"Active_Rxns={active_count} ✓")

        # Validate DGOA flux if applicable
        if dgoa_status == "with_dgoa":
            dgoa_flux = solution_dict.get("dgoa_flux")
            if dgoa_flux is None or dgoa_flux <= 0:
                validation_messages.append(f"DGOA_Flux={dgoa_flux} (FAIL: ≤0)")
                passed = False
            else:
                validation_messages.append(f"DGOA_Flux={dgoa_flux:.4f} ✓")
        else:
            validation_messages.append("DGOA_Flux=N/A")

        # Construct final message
        status_symbol = "✓" if passed else "❌"
        message = f"{status_symbol} {strain}_{dgoa_status}: " + ", ".join(validation_messages)

        return passed, message

    def create_expression_flux_summary(self, results_dict):
        """
        Create summary dictionary from expression flux analysis results.

        Reorganizes results from flat structure (strain_dgoa_status keys) to nested
        structure organized by strain, then DGOA status.

        Args:
            results_dict (dict): Results dictionary with keys like "ACN2586_with_dgoa", etc.
                Each value is a dict from process_strain_with_expression()

        Returns:
            dict: Nested summary with structure:
                {
                    "strain_name": {
                        "with_dgoa": {
                            "biomass": float,
                            "active_reactions": int,
                            "off_reactions": list,
                            "on_reactions": list,
                            "dgoa_flux": float
                        },
                        "without_dgoa": {
                            "biomass": float,
                            "active_reactions": int,
                            "off_reactions": list,
                            "on_reactions": list
                        }
                    }
                }
        """
        summary = {}

        for condition_key, result_data in results_dict.items():
            # Parse condition key (e.g., "ACN2586_with_dgoa" -> strain="ACN2586", dgoa="with_dgoa")
            if "_with_dgoa" in condition_key:
                strain = condition_key.replace("_with_dgoa", "")
                dgoa_status = "with_dgoa"
            elif "_without_dgoa" in condition_key:
                strain = condition_key.replace("_without_dgoa", "")
                dgoa_status = "without_dgoa"
            else:
                self.logger.warning(f"Unexpected condition key format: {condition_key}")
                continue

            # Initialize strain entry if needed
            if strain not in summary:
                summary[strain] = {}

            # Create summary entry for this condition
            condition_summary = {
                "biomass": result_data.get("biomass", 0),
                "active_reactions": result_data.get("active_reactions", 0),
                "off_reactions": result_data.get("off_reactions", []),
                "on_reactions": result_data.get("on_reactions", []),
                "solution_status": result_data.get("solution_status", "unknown"),
                "dgoa_flux": result_data.get("dgoa_flux", 0)  # Include for both conditions for consistency
            }

            # Add error if present
            if "error" in result_data:
                condition_summary["error"] = result_data["error"]

            # Store in summary
            summary[strain][dgoa_status] = condition_summary

        self.logger.info(f"Created summary for {len(summary)} strains")

        return summary

    def export_expression_flux_to_excel(self, results_dict, output_file, base_model):
        """
        Export expression flux analysis results to multi-sheet Excel workbook.

        Creates a summary sheet with overview statistics and individual sheets
        for each condition with detailed flux data.

        Args:
            results_dict (dict): Results dictionary with keys like "ACN2586_with_dgoa"
            output_file (str): Output filename (relative to notebook folder or absolute path)
            base_model (MSModelUtil): Base metabolic model for extracting reaction metadata

        Returns:
            str: Full path to created Excel file
        """
        try:
            # Ensure output is in nboutput directory
            if not output_file.startswith("/"):
                output_file = f"{self.notebook_folder}/nboutput/{output_file}"

            # Create summary data for the Summary sheet
            summary_rows = []
            for condition_key, result_data in results_dict.items():
                # Parse condition key
                if "_with_dgoa" in condition_key:
                    strain = condition_key.replace("_with_dgoa", "")
                    dgoa_status = "with_dgoa"
                elif "_without_dgoa" in condition_key:
                    strain = condition_key.replace("_without_dgoa", "")
                    dgoa_status = "without_dgoa"
                else:
                    continue

                row = {
                    "Strain": strain,
                    "DGOA_Status": dgoa_status,
                    "Biomass_Flux": result_data.get("biomass", 0),
                    "Active_Reactions_Count": result_data.get("active_reactions", 0),
                    "Off_Reactions_Count": len(result_data.get("off_reactions", [])),
                    "On_Reactions_Count": len(result_data.get("on_reactions", [])),
                    "DGOA_Flux": result_data.get("dgoa_flux", "N/A"),
                    "Solution_Status": result_data.get("solution_status", "unknown")
                }
                summary_rows.append(row)

            summary_df = pd.DataFrame(summary_rows)

            # Create Excel writer
            with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
                # Write summary sheet
                summary_df.to_excel(writer, sheet_name="Summary", index=False)
                self.logger.info(f"Wrote Summary sheet with {len(summary_rows)} rows")

                # Write individual condition sheets
                for condition_key, result_data in results_dict.items():
                    fluxes = result_data.get("fluxes", {})
                    if not fluxes:
                        self.logger.warning(f"No fluxes for {condition_key}, skipping sheet")
                        continue

                    # Create sheet name (max 31 chars for Excel)
                    sheet_name = condition_key.replace("_", " ").title()
                    if len(sheet_name) > 31:
                        # Abbreviate if too long
                        sheet_name = sheet_name.replace("With Dgoa", "W_DGOA").replace("Without Dgoa", "WO_DGOA")
                    if len(sheet_name) > 31:
                        sheet_name = sheet_name[:31]

                    # Create rows for this condition
                    condition_rows = []
                    for rxn_id, flux_value in fluxes.items():
                        # Get reaction metadata from model
                        rxn_name = "N/A"
                        gene_association = "N/A"
                        try:
                            rxn = base_model.model.reactions.get_by_id(rxn_id)
                            rxn_name = rxn.name if rxn.name else rxn_id
                            # Get gene association
                            genes = [str(g) for g in rxn.genes if not str(g).startswith("mRNA_")]
                            gene_association = "; ".join(genes) if genes else "spontaneous"
                        except:
                            # Reaction not in model (e.g., DgoA added separately)
                            pass

                        condition_rows.append({
                            "Reaction_ID": rxn_id,
                            "Reaction_Name": rxn_name,
                            "Flux": flux_value,
                            "Gene_Association": gene_association
                        })

                    # Sort by absolute flux value (largest first)
                    condition_rows.sort(key=lambda x: abs(x["Flux"]), reverse=True)

                    # Create DataFrame and write to sheet
                    condition_df = pd.DataFrame(condition_rows)
                    condition_df.to_excel(writer, sheet_name=sheet_name, index=False)
                    self.logger.info(f"Wrote sheet '{sheet_name}' with {len(condition_rows)} reactions")

            self.logger.info(f"Excel file created: {output_file}")
            return output_file

        except Exception as e:
            self.logger.error(f"Error creating Excel file: {str(e)}")
            raise

    def generate_all_escher_maps(self, results_dict, models_dict, map_name="Core", base_model_file=None):
        """
        Generate Escher metabolic maps for all conditions.

        Creates interactive HTML visualizations showing flux distributions
        for each strain and DGOA condition.

        Args:
            results_dict (dict): Results dictionary with keys like "ACN2586_with_dgoa"
            models_dict (dict): Dictionary mapping condition keys to model objects
                (used to extract reaction metadata for the maps)
            map_name (str): Name of Escher map to use (default: "Core")
            base_model_file (str): Path to base model file to use when models_dict is None
                (default: "data/FullyTranslatedPublishedModel.json")

        Returns:
            list: List of successfully created file paths
        """
        created_files = []
        failed_conditions = []
        
        # Load base model if models_dict is None and base_model_file is provided
        base_model = None
        if models_dict is None:
            if base_model_file is None:
                # Try default model file
                default_model_file = "data/FullyTranslatedPublishedModel.json"
                if os.path.exists(default_model_file):
                    base_model_file = default_model_file
                else:
                    # Try alternative
                    alt_model_file = "data/TranslatedPublishedModel.json"
                    if os.path.exists(alt_model_file):
                        base_model_file = alt_model_file
            
            if base_model_file and os.path.exists(base_model_file):
                try:
                    base_model = MSModelUtil.from_cobrapy(base_model_file)
                    self.logger.info(f"Loaded base model from {base_model_file} for Escher maps")
                except Exception as e:
                    self.logger.warning(f"Could not load base model from {base_model_file}: {e}")

        for condition_key, result_data in results_dict.items():
            try:
                # Get fluxes for this condition
                fluxes = result_data.get("fluxes", {})
                if not fluxes:
                    self.logger.warning(f"No fluxes for {condition_key}, skipping Escher map")
                    failed_conditions.append((condition_key, "No flux data"))
                    continue

                # Parse condition key for title and filename
                if "_with_dgoa" in condition_key:
                    strain = condition_key.replace("_with_dgoa", "")
                    dgoa_status = "with_dgoa"
                    title = f"{strain} - With DGOA - Pyruvate Media"
                elif "_without_dgoa" in condition_key:
                    strain = condition_key.replace("_without_dgoa", "")
                    dgoa_status = "without_dgoa"
                    title = f"{strain} - Without DGOA - Pyruvate Media"
                else:
                    self.logger.warning(f"Unexpected condition key format: {condition_key}")
                    failed_conditions.append((condition_key, "Invalid key format"))
                    continue

                # Create output filename
                output_file = f"nboutput/escher_{strain}_{dgoa_status}.html"

                # Get model for this condition (if available)
                model = models_dict.get(condition_key) if models_dict else base_model

                # Generate Escher map using inherited create_map_html method
                if model:
                    # Use model object if available
                    cobra_model = model.model if hasattr(model, 'model') else model
                    result_path = self.create_map_html(
                        model=cobra_model,
                        flux_solution=fluxes,
                        map_name=map_name,
                        output_file=output_file,
                        title=title
                    )
                else:
                    # Create map without model (just flux data) - this may fail if map requires model
                    # Try to use map_name only
                    try:
                        result_path = self.create_map_html(
                            model=None,
                            flux_solution=fluxes,
                            map_name=map_name,
                            output_file=output_file,
                            title=title
                        )
                    except Exception as e:
                        error_msg = f"Cannot create map without model: {str(e)}"
                        self.logger.error(error_msg)
                        failed_conditions.append((condition_key, error_msg))
                        continue

                created_files.append(result_path)
                self.logger.info(f"Created Escher map: {output_file}")

            except Exception as e:
                error_msg = f"Error creating Escher map for {condition_key}: {str(e)}"
                self.logger.error(error_msg)
                failed_conditions.append((condition_key, str(e)))
                # Continue to next condition rather than failing completely
                continue

        # Log summary
        self.logger.info(f"Successfully created {len(created_files)} Escher maps")
        if failed_conditions:
            self.logger.warning(f"Failed to create maps for {len(failed_conditions)} conditions:")
            for cond, reason in failed_conditions:
                self.logger.warning(f"  - {cond}: {reason}")

        return created_files

    def run_expression_flux_analysis(self, strains, proteomics_file, model_file, media_id,
                                     knockout_dahp=True, biomass_fraction=0.25,
                                     default_coef=0.01, activation_threshold=None,
                                     deactivation_threshold=0.000001, minimal_flux=0.001):
        """
        Main orchestrator for expression-constrained flux analysis pipeline.

        Runs complete analysis: loads data, averages replicates, processes all conditions,
        caches intermediate results, and returns comprehensive results.

        Args:
            strains (list): List of strain names (e.g., ["ACN2586", "ACN2821", ...])
            proteomics_file (str): Path to proteomics Excel file
            model_file (str): Path to metabolic model JSON file
            media_id (str): KBase media ID (e.g., "KBaseMedia/Carbon-Pyruvic-Acid")
            knockout_dahp (bool): Whether to knock out DAHP synthase (default: True)
            biomass_fraction (float): Minimum fraction of optimal biomass (default: 0.25)
            default_coef (float): Default coefficient for expression fitting (default: 0.01)
            activation_threshold (float): Threshold for activation (default: None)
            deactivation_threshold (float): Threshold for deactivation (default: 0.000001)
            minimal_flux (float): Minimal flux for on_on reactions (default: 0.001)

        Returns:
            dict: Complete results dictionary with keys for each condition (strain_dgoa_status)
        """
        self.logger.info("=" * 80)
        self.logger.info("Starting Expression-Constrained Flux Analysis Pipeline")
        self.logger.info("=" * 80)

        # Step 1: Load base model
        self.logger.info(f"Step 1: Loading metabolic model from {model_file}")
        base_model = MSModelUtil.from_cobrapy(model_file)
        base_model.util = self
        self.logger.info(f"  Model loaded: {len(base_model.model.reactions)} reactions, {len(base_model.model.metabolites)} metabolites")

        # Step 2: Load media
        self.logger.info(f"Step 2: Loading media: {media_id}")
        media = self.get_media(media_id)
        self.logger.info(f"  Media loaded")

        

        # Cache averaged expression
        

        # Step 5: Process all strain × DGOA conditions
        self.logger.info(f"Step 5: Processing {len(strains)} strains × 2 DGOA variants = {len(strains) * 2} conditions")
        results = {}
        models_dict = {}

        dgoa_variants = [True, False]  # with_dgoa, without_dgoa

        for i, strain in enumerate(strains, 1):
            self.logger.info(f"\n  [{i}/{len(strains)}] Processing strain: {strain}")

            for with_dgoa in dgoa_variants:
                dgoa_status = "with_dgoa" if with_dgoa else "without_dgoa"
                condition_key = f"{strain}_{dgoa_status}"

                self.logger.info(f"    Processing {condition_key}...")

                try:
                    # Process this condition
                    result = self.process_strain_with_expression(
                        strain=strain,
                        expression_data=averaged_expression,
                        base_model=base_model,
                        media=media,
                        with_dgoa=with_dgoa,
                        knockout_dahp=knockout_dahp,
                        biomass_fraction=biomass_fraction,
                        default_coef=default_coef,
                        activation_threshold=activation_threshold,
                        deactivation_threshold=deactivation_threshold,
                        minimal_flux=minimal_flux
                    )

                    results[condition_key] = result

                    # Log result summary
                    if result.get("solution_status") == "optimal":
                        self.logger.info(f"      ✓ Biomass={result['biomass']:.4f}, Active={result['active_reactions']}, Status={result['solution_status']}")
                    else:
                        self.logger.warning(f"      ✗ Status={result['solution_status']}, Biomass={result['biomass']:.4f}")

                    # Cache individual condition result
                    self.save(f"expr_flux_{condition_key}", result)

                except Exception as e:
                    self.logger.error(f"      ✗ Error: {str(e)}")
                    results[condition_key] = {
                        "fluxes": {},
                        "biomass": 0,
                        "active_reactions": 0,
                        "off_reactions": [],
                        "on_reactions": [],
                        "dgoa_flux": None,
                        "solution_status": "error",
                        "error": str(e)
                    }

        # Step 6: Create summary
        self.logger.info(f"\nStep 6: Creating analysis summary")
        summary = self.create_expression_flux_summary(results)
        self.save("expression_flux_summary", summary)
        self.logger.info(f"  Summary created and cached")

        # Step 7: Cache complete results
        self.logger.info(f"Step 7: Caching complete results")
        # Save only the non-flux data to avoid huge JSON files
        results_metadata = {}
        for key, val in results.items():
            results_metadata[key] = {k: v for k, v in val.items() if k != "fluxes"}
        self.save("expression_flux_results_metadata", results_metadata)
        self.logger.info(f"  Results metadata cached")

        # Final summary
        self.logger.info("\n" + "=" * 80)
        self.logger.info("Pipeline Complete!")
        self.logger.info(f"  Total conditions processed: {len(results)}")
        successful = sum(1 for r in results.values() if r.get("solution_status") == "optimal")
        self.logger.info(f"  Successful optimizations: {successful}/{len(results)}")
        self.logger.info(f"  Failed optimizations: {len(results) - successful}/{len(results)}")
        self.logger.info("=" * 80)

        return results

# Initialize the NotebookUtil instance
util = NotebookUtil() 