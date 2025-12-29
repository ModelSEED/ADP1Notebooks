# ADP1 Research Project - Codebase Documentation

**Date**: 2025-10-22
**Research Type**: Codebase Structure Analysis
**Project**: Acinetobacter baylyi ADP1 Genome and Metabolic Model Analysis

---

## Executive Summary

This project contains Jupyter notebooks and supporting code for analyzing Acinetobacter baylyi ADP1 genome data, metabolic models, and proteomics. The research focuses on:

1. **Metabolic model comparison** - Comparing published models with ModelSEED-generated models
2. **DGOA enzyme analysis** - Analyzing the DgoA enzyme pathway and its role in aromatic amino acid biosynthesis
3. **Gene annotation analysis** - Integrating multiple annotation sources (RAST, DRAM, Snekmer, GLM4EC)
4. **Proteomics integration** - Incorporating UGA proteomics data with genomic and metabolic model data
5. **AI-assisted curation** - Using AI to validate gene-reaction associations and reaction directionality

---

## Project Structure

```
ADP1Research/
├── ADP1Notebooks/          # Symlink to ~/projects/ADP1Notebooks
│   ├── notebooks/          # Main analysis notebooks
│   │   ├── ADP1ModelReview.ipynb
│   │   ├── ADP1DGOAAnalysis.ipynb
│   │   ├── ADP1AIGeneAssociation.ipynb
│   │   ├── ADP1AnnotationAnalysis.ipynb
│   │   ├── ADP1AIModelCompare.ipynb
│   │   ├── ADP1ExpressionAnalysis.ipynb
│   │   ├── util.py         # Core utility functions
│   │   ├── vibe.md         # Brainstorming template
│   │   ├── data/           # Input data files
│   │   ├── datacache/      # Cached processed data
│   │   └── nboutput/       # Output files
│   ├── docs/               # Documentation
│   └── README.md
└── .claude/                # Claude Code configuration
```

---

## Core Notebooks

### 1. ADP1ModelReview.ipynb

**Purpose**: Standardize and compare the published ADP1 model with the ModelSEED model.

**Key Operations**:
- Loads published model from XML: `data/PublishedModel.XML`
- Loads KBase genome annotation: `Acinetobacter_baylyi_ADP1_RPG_Snekmer_LA_Pfam` (workspace 219173)
- Translates gene IDs between published model and KBase annotation
- Performs model standardization and comparison
- Uses AI to review model reactions for errors in directionality, definition, or stoichiometry
- Evaluates impact of directionality changes on flux solutions

**Key Outputs**:
- `gene_translation.json` - Gene ID translation mapping
- `match_stats.json` - Reaction matching statistics
- `model_match_stats.json` - Comprehensive model comparison stats
- `mdl_rxn_hash.json` - Parsed reaction stoichiometry

**ModelSEED Model**: `179225/Abaylyi_ADP1_RASTMS2_OMEGGA_Abaylyi_ADP1_RAST.mdlMS2_OMEGGA_iAbaylyi_Carbon_Succinic.gf`

---

### 2. ADP1DGOAAnalysis.ipynb

**Purpose**: Analyze the DgoA enzyme and its impact on metabolic flux under different conditions.

**Background**: The DgoA enzyme provides an alternative pathway for aromatic amino acid biosynthesis when the native DAHP synthase (rxn01332) is knocked out.

**Key Operations**:
- Adds DGOA reaction to the model: `cpd00236 + cpd00020 -> cpd02857`
- Simulates wild-type (rxn01332 active, DGOA inactive)
- Simulates DGOA mutant (rxn01332 inactive, DGOA active)
- Uses MOMA (Minimization of Metabolic Adjustment) to predict metabolic adaptations
- Analyzes flux differences between conditions
- Tests multiple knockout combinations to capture different adaptation strategies

**Media Condition**: Pyruvate media (`KBaseMedia/Carbon-Pyruvic-Acid`)

**Key Outputs**:
- `DGOA_transition_analysis.xlsx` - Comprehensive reaction flux comparison
- `ms_fva.json` - Flux variability analysis for ModelSEED model
- `ms_active_fluxes.json` - Active flux values
- `pubmod_dgoa_fva.json` - Published model FVA with DGOA
- `pubmod_dgoa_fluxes.json` - Published model fluxes with DGOA
- `escher_map.html` - Metabolic map visualization

**Strain Data**:
- ACN2586: Initial construct with dgoA* and native DAHP synthesis deleted
- ACN2821: Evolved strain with single copy of dgoA*
- ACN3015: Single copy of dgoA after evolution
- ACN3468: Multiple copies of dgoA*
- ACN3471: Multiple copies of dgoA
- ACN3474: Multiple copies of dgoA, partially evolved
- ACN3477: Multiple copies of dgoA, more evolved
- ADP1: Wild-type strain

---

### 3. ADP1AIGeneAssociation.ipynb

**Purpose**: Use AI to evaluate gene-reaction associations in the metabolic model.

**Key Operations**:
- Loads published model and reaction data
- Loads enhanced gene annotation data
- For each reaction, evaluates whether associated genes match the reaction chemistry
- Uses `util.evaluate_reaction_gene_association()` with AI
- Removes proteomics data before AI analysis (to focus on functional annotation)

**Key Outputs**:
- `rxn_gene_mapping.json` - AI-validated gene-reaction associations

---

### 4. ADP1AnnotationAnalysis.ipynb

**Purpose**: Integrate multiple gene annotation sources and proteomics data into a unified dataset.

**Key Operations**:

1. **Pull KBase genome annotations**:
   - Genome reference: `219173/Acinetobacter_baylyi_ADP1_RPG_Snekmer_LA_Pfam`
   - Uses AnnotationOntology API
   - Saves to `gene_term_hash.json`

2. **Add ModelSEED model data**:
   - Simulates model on pyruvate media
   - Adds reaction associations and flux data
   - Label: "modelseed"

3. **Add published model data**:
   - Model: `179225/iAbaylyi_ADP1`
   - Adds reaction associations and flux data
   - Label: "published"

4. **Add proteomics data**:
   - Source: `data/UGA_Proteomics_May2025_Report.xlsx`
   - Integrates Tukey test results (group comparisons)
   - Adds imputed proteome values for ACN2586 and ACN2821

5. **Add protein cluster data**:
   - Source: `data/Protein_HC_Cluster_Assignments.csv`
   - Source: `data/ClusterSizes.txt`

6. **Translate term IDs to names**:
   - Uses dictionaries: EC, KO, PF, SSO
   - Creates human-readable annotations

7. **Generate spreadsheets**:
   - `gene_annotations.xlsx` - Comprehensive gene annotation table
   - `data/ADP1Genes.csv` - Gene coordinates with annotations

**Key Outputs**:
- `gene_term_hash.json` - Raw annotation data by source
- `gene_term_hash_named.json` - Annotations with term names
- `genedata_enhanced.json` - Enhanced gene data with proteomics
- `gene_annotations.xlsx` - Master gene annotation spreadsheet
- `ADP1Genes.csv` - Gene table for GFF generation

**Annotation Sources**:
- RAST: Subsystem-based annotation
- GLM4EC: EC number prediction via AI
- DRAM: KEGG orthology annotation
- Snekmer: Pfam domain annotation

---

### 5. ADP1AIModelCompare.ipynb

**Purpose**: Use AI to evaluate reaction equivalence between published model and ModelSEED reactions.

**Key Operations**:
- Loads published model from XML
- Loads reaction comparison data
- For each published reaction, evaluates equivalence with ModelSEED matches
- Uses `util.evaluate_reaction_equivalence()` with AI
- Caches AI responses to avoid redundant API calls

**Key Outputs**:
- `rxn_mapping.json` - Reaction equivalence evaluations

---

### 6. ADP1ExpressionAnalysis.ipynb

**Purpose**: Analyze gene expression and proteomics data (notebook exists but content not fully explored).

---

## Core Utility Module: util.py

### Location
`ADP1Notebooks/notebooks/util.py`

### Class: NotebookUtil

**Inheritance**: `MSFBAUtils, AICurationUtils, NotebookUtils, KBGenomeUtils, EscherUtils`

**Dependencies**:
- KBUtilLib (ModelSEED utilities)
- cobrakbase (KBase-specific COBRApy extensions)
- modelseedpy (ModelSEED Python package)
- cobra (Constraint-based reconstruction and analysis)
- pandas (Data manipulation)
- escher (Metabolic map visualization)

**Key Configuration**:
```python
notebook_folder = script_dir
name = "ADP1NotebookUtils"
user = "chenry"
retries = 5
proxy_port = None
```

### Key Methods

#### 1. `get_model_and_simulate(model_id, media_id)`
Loads a KBase model and simulates it with specified media.

**Process**:
1. Loads model from KBase
2. Checks for cpd00236 (galacturonate) - if present, adds DGOA reaction
3. Sets media conditions
4. Runs pFBA (parsimonious flux balance analysis)
5. For each active reaction, performs knockout analysis
6. Calculates knockout impact ratio

**Returns**: Dictionary with model, reaction flux data, and solution

#### 2. `add_model_data_to_annotations(model, annotations, rxn_flux_data, label)`
Adds metabolic model information to gene annotations.

**Parameters**:
- `model`: FBAModel object
- `annotations`: Gene annotation dictionary
- `rxn_flux_data`: Reaction flux values
- `label`: Identifier for this model (e.g., "modelseed", "published")

**Process**:
- For each reaction in model:
  - Creates reaction string with equation and names
  - For each gene associated with reaction:
    - Adds reaction to gene's annotation
    - Adds flux data (flux value and KO ratio)

#### 3. `add_dgoa_reaction(model)`
Adds the DgoA reaction to a metabolic model.

**Reaction**: `cpd00236 + cpd00020 -> cpd02857`
- cpd00236: D-galacturonate
- cpd00020: Pyruvate
- cpd02857: 3-deoxy-D-arabino-hept-2-ulosonate 7-phosphate

**Properties**:
- ID: 'DgoA'
- Gene rule: 'DgoA'
- Bounds: 0 to 1000 (irreversible, forward only)

#### 4. `reaction_equation_with_names(rxn)`
Formats reaction equations using metabolite names instead of IDs.

**Output Format**: `"2 ATP + glucose <=> 2 ADP + glucose-6-phosphate"`

#### 5. `find_significant_differences(flux1, flux2, threshold=0.01)`
Identifies reactions with significantly different fluxes between two conditions.

**Returns**: List of dictionaries with reaction ID, fluxes, and difference

#### 6. `reactions_to_genes(reaction_list, model)`
Converts a list of reactions to their associated genes.

**Process**:
- For each reaction, extracts gene list
- Filters out mRNA features
- Creates reaction-to-gene mapping

**Returns**: Gene list and reaction-gene mapping dictionary

---

## Data Files

### Input Data (`notebooks/data/`)

1. **PublishedModel.XML** - Original published ADP1 metabolic model
2. **PublishedModel.json** - JSON version of published model
3. **TranslatedPublishedModel.json** - Published model with ModelSEED compound IDs
4. **UGA_Proteomics_May2025_Report.xlsx** - Proteomics data with multiple sheets:
   - Imputed: Imputed protein abundance values
   - Tukey: Statistical comparisons between strains
5. **MutantGrowthRatesData.xls** - Growth phenotype data
6. **ADP1Genes.csv** - Gene coordinate and annotation table
7. **ClusterSizes.txt** - Protein cluster size data
8. **Protein_HC_Cluster_Assignments.csv** - Protein hierarchical clustering
9. **EllenGeneCalls.tsv** - Manual gene coordinate curation

### Cached Data (`notebooks/datacache/`)

All cached data stored as JSON files:

**Genome & Annotation**:
- `genedata.json` - Gene data from model comparison
- `genedata_enhanced.json` - Gene data with proteomics
- `gene_term_hash.json` - Raw gene-term associations
- `gene_term_hash_named.json` - Gene terms with readable names
- `gene_translation.json` - Gene ID mappings

**Reaction Data**:
- `reactiondata.json` - Reaction comparison data
- `rxn_hash.json` - Reaction stoichiometry hash
- `rxn_id_hash.json` - Reaction ID mappings
- `rxn_mapping.json` - AI-evaluated reaction equivalences
- `rxn_gene_mapping.json` - AI-evaluated gene associations
- `mdl_rxn_hash.json` - Model reaction hash

**Flux Data**:
- `ms_fluxes.json` - ModelSEED model flux solutions
- `ms_fva.json` - ModelSEED FVA results
- `pubmod_fluxes.json` - Published model fluxes
- `pubmod_fva.json` - Published model FVA
- `pubmod_active_fluxes.json` - Active fluxes in published model
- `pubmod_dgoa_fluxes.json` - Published model with DGOA
- `pubmod_dgoa_fva.json` - Published model DGOA FVA

**Comparison Data**:
- `match_stats.json` - Reaction matching statistics
- `model_match_stats.json` - Model comparison statistics

### Output Data (`notebooks/nboutput/`)

- `flux_distribution.csv` - Gene-reaction flux distribution
- `escher_map.html` - Interactive metabolic map visualization

---

## Key External Dependencies

### KBase Workspace References

1. **Genome**: `219173/Acinetobacter_baylyi_ADP1_RPG_Snekmer_LA_Pfam`
2. **ModelSEED Model**: `179225/Abaylyi_ADP1_RASTMS2_OMEGGA_Abaylyi_ADP1_RAST.mdlMS2_OMEGGA_iAbaylyi_Carbon_Succinic.gf`
3. **Published Model**: `179225/iAbaylyi_ADP1`
4. **Media**: `KBaseMedia/Carbon-Pyruvic-Acid`

### External Libraries

- **KBUtilLib**: Custom KBase utility library (path: `/home/chenry/MyEnvs/modelseed_cplex`)
  - MSFBAUtils: FBA modeling utilities
  - AICurationUtils: AI-assisted curation
  - NotebookUtils: Jupyter notebook helpers
  - KBGenomeUtils: Genome data utilities
  - EscherUtils: Metabolic map visualization

- **cobrakbase**: KBase extensions for COBRApy (version 0.4.0)
  - FBAModel class for KBase models

- **modelseedpy**: ModelSEED Python package (version 0.4.2)
  - AnnotationOntology: Genome annotation handling
  - MSPackageManager: Package management
  - MSMedia: Media management
  - MSModelUtil: Model utilities
  - MSBuilder: Model building
  - MSATPCorrection: ATP correction
  - MSGapfill: Gap-filling
  - MSGrowthPhenotype(s): Phenotype simulation
  - ModelSEEDBiochem: Biochemistry database
  - MSExpression: Expression data integration

- **cobra**: Constraint-based modeling (COBRApy)
  - Reaction, Metabolite classes
  - pfba: Parsimonious FBA
  - moma: Minimization of Metabolic Adjustment

---

## Research Workflow

### Typical Analysis Pipeline

1. **Load and translate model**:
   - Load published model from XML
   - Get KBase genome annotation
   - Translate gene IDs
   - Compare with ModelSEED model

2. **Simulate metabolic conditions**:
   - Set media (pyruvate)
   - Run pFBA
   - Perform FVA (flux variability analysis)
   - Add DGOA reaction if needed

3. **Integrate annotations**:
   - Pull annotations from KBase
   - Add model associations
   - Add flux data
   - Add proteomics
   - Translate term IDs to names

4. **AI-assisted curation**:
   - Evaluate gene-reaction associations
   - Check reaction directionality
   - Validate reaction equivalences

5. **Generate outputs**:
   - Create spreadsheets
   - Generate visualizations
   - Export Escher maps

---

## AI Integration

The project uses AI extensively through `AICurationUtils`:

1. **Reaction Directionality**: `ai_analysis_of_model_reactions()`
   - Analyzes biochemical directionality
   - Compares with model directionality
   - Flags conflicts

2. **Gene Associations**: `evaluate_reaction_gene_association()`
   - Validates gene-reaction links
   - Uses gene annotations
   - Returns confidence scores

3. **Reaction Equivalence**: `evaluate_reaction_equivalence()`
   - Compares reactions from different sources
   - Evaluates stoichiometric equivalence
   - Considers compartmentalization

4. **Caching**: AI responses are cached to avoid redundant API calls
   - Location: Managed by KBUtilLib
   - Cache key: Based on function name and parameters

---

## Notable Analysis Features

### 1. DGOA Transition Analysis

The notebooks perform sophisticated analysis of metabolic adaptation:
- Wild-type baseline with native DAHP synthase
- DGOA mutant with alternative pathway
- MOMA-based prediction of metabolic rewiring
- Multiple knockout combinations to explore adaptation space

### 2. Multi-Source Annotation Integration

Combines 4+ annotation sources:
- RAST (subsystem-based)
- DRAM (KEGG-based)
- Snekmer (domain-based)
- GLM4EC (AI-predicted EC numbers)
- Proteomics (experimental)

### 3. Model Comparison

Systematic comparison of two models:
- Published model (literature)
- ModelSEED model (automated reconstruction)
- Gene matching
- Reaction matching
- Flux comparison

---

## Python Helper Scripts

Located in `notebooks/`:

1. **debug_gene_matching.py** - Debugging gene ID matching
2. **test_fix.py** - Testing model fixes
3. **check_gene_annotations.py** - Validating annotations
4. **test_with_genome.py** - Testing with genome data

---

## Development Notes

### Environment Setup

The notebooks expect:
- Jupyter environment
- KBase authentication token at `~/.kbase/token`
- KBUtilLib installed at specific path
- ModelSEED database locally available
- Python 3.9
- GLPK solver (warnings indicate lack of optimality tolerance support)

### Common Warnings

1. **"Section 'ModelSEEDBiochem' not found in config"** - Expected, uses default
2. **"Conditional Formatting extension is not supported"** - From openpyxl, safe to ignore
3. **"NUMEXPR_MAX_THREADS not set"** - NumExpr defaults to 8 threads

### Solver Notes

Uses GLPK solver by default, which doesn't support setting optimality tolerance.

---

## Research Questions Addressed

1. **How does the published ADP1 model compare to the ModelSEED-generated model?**
   - Gene coverage comparison
   - Reaction matching and equivalence
   - Flux solution differences

2. **What are the metabolic consequences of replacing DAHP synthase with DgoA?**
   - Flux redistribution analysis
   - Identification of compensatory pathways
   - Growth impact assessment

3. **Which genes are differentially expressed between ADP1 strains?**
   - Proteomics data integration
   - Statistical comparison (Tukey tests)
   - Clustering analysis

4. **How accurate are automated gene annotations?**
   - Multi-source comparison
   - AI-based validation
   - Subsystem assignment

5. **Which reactions have uncertain directionality?**
   - AI analysis of biochemistry
   - Comparison with model constraints
   - Impact on flux solutions

---

## File Format Notes

### Cached JSON Structure

**gene_term_hash_named.json**:
```json
{
  "gene_id": {
    "annotation_source": {
      "term_id": "term_name"
    },
    "ACN2586_vs_ACN2821": {"value": "0.5"},
    "modelseed_flux": ["rxn_id:flux:fva"],
    "cluster": "cluster_id",
    "cluster_size": 10
  }
}
```

**rxn_flux_data**:
```json
{
  "rxn_id": {
    "flux": 1.23,
    "ko_ratio": 0.85
  }
}
```

---

## Future Directions

Based on notebook structure and data:

1. Expression analysis (ADP1ExpressionAnalysis.ipynb needs exploration)
2. Additional strain comparisons (ACN3015, ACN3468, ACN3471, ACN3474, ACN3477)
3. Escher map refinement (multiple attempts visible in DGOA notebook)
4. GFF file generation and editing
5. Integration with BioCyc gene coordinates

---

## Contact & Attribution

**User**: chenry
**ModelSEED Database**: Used for compound and reaction biochemistry
**KBase**: Used for genome data and model storage
**AI Gateway**: ArgoGatewayClient (dev environment)

---

## Conclusion

This is a comprehensive systems biology research project combining:
- Metabolic modeling (FBA, FVA, MOMA)
- Genome annotation (multi-source integration)
- Proteomics (experimental validation)
- AI-assisted curation (reaction and gene validation)
- Strain evolution analysis (multiple ADP1 variants)

The primary focus is understanding how Acinetobacter baylyi ADP1 can adapt to loss of native DAHP synthase through the DgoA alternative pathway, with implications for metabolic engineering and synthetic biology.
