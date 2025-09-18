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
sys.path = [base_dir+"/KBUtilLib/src",base_dir+"/KBUtilLib/src",base_dir+"/KBUtilLib/src/kbutillib/dependencies/cobrakbase",base_dir+"/KBUtilLib/src/kbutillib/dependencies/ModelSEEDpy"] + sys.path

# Import utilities with error handling
from kbutillib import MSFBAUtils, AICurationUtils, NotebookUtils, KBGenomeUtils

import hashlib
import pandas as pd
from pandas import DataFrame, read_csv, concat, set_option
from cobrakbase.core.kbasefba import FBAModel
import cobra
from cobra import Reaction, Metabolite
from cobra.flux_analysis import pfba
from modelseedpy import AnnotationOntology, MSPackageManager, MSMedia, MSModelUtil, MSBuilder, MSATPCorrection, MSGapfill, MSGrowthPhenotype, MSGrowthPhenotypes, ModelSEEDBiochem, MSExpression
from modelseedpy.core.msprobability import MSProbability
from modelseedpy.core.annotationontology import convert_to_search_role, split_role
from modelseedpy.core.mstemplate import MSTemplateBuilder
from modelseedpy.core.msgenome import normalize_role
from modelseedpy.core.msensemble import MSEnsemble
from modelseedpy.community.mscommunity import MSCommunity
from modelseedpy.helpers import get_template

# Define the base classes based on what's available
class NotebookUtil(MSFBAUtils, AICurationUtils, NotebookUtils, KBGenomeUtils):
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
        model.set_media(media)
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
        cpd00236 = model.model.metabolites.get_by_id('cpd00236_c0')
        cpd00020 = model.model.metabolites.get_by_id('cpd00020_c0')
        cpd02857 = model.model.metabolites.get_by_id('cpd02857_c0')

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
    
# Initialize the NotebookUtil instance
util = NotebookUtil() 