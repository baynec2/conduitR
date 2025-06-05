# Importing necessary packages
import os
import glob
import pandas as pd
import shutil
import pdb  # Add this for debugging
import logging
from datetime import datetime

# Load configuration from command line
# Usage: snakemake --configfile path/to/config.yaml --sdm apptainer
if not config.get("experiment_directory"):
    raise ValueError("Please provide experiment_directory in config file")

# Get experiment directory from config
EXPERIMENT_DIR = config["experiment_directory"]

# Print the expected config file path
expected_config_path = os.path.join(EXPERIMENT_DIR, "config/run_diann.cfg")

# Define a function to get rule-specific log paths
def get_log_path(rule_name):
    return os.path.join(EXPERIMENT_DIR, f"logs/{rule_name}.log")

# Define a function to setup logging for a rule
def setup_rule_logging(log_file):
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='[%(asctime)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    return logging.getLogger()

# Extracting the method from the config file
if not config.get("method"):
    raise ValueError("Please provide 'method' in config file")
# Defining the allowed methods. 
ALLOWED_METHODS = [
    "ncbi_taxonomy_id" # Will uncomment methods as they become supported,
   # "uniprot_proteome_id",
   # "proteotyping",
   # "MAG",
   # "metagenomic_profiling",
   # "16S"
]
# Checking that the method is allowed.   
METHOD = config["method"]
if METHOD not in ALLOWED_METHODS:
    raise ValueError(f"Method '{METHOD}' not allowed. Must be one of: {', '.join(ALLOWED_METHODS)}")

# Defining the .raw samples that are inputs to the workflow using config
if not config.get("raw_files"):
    raise ValueError("Please provide raw_files pattern in config file")
raw_files_pattern = config["raw_files"]

# Get the full paths to all raw files by joining experiment dir with pattern and then using glob
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/raw_files/*.raw"))
# Get just the base filenames without extension for SAMPLES
SAMPLES = []  # Initialize empty list
for filepath in RAW_FILEPATHS:
    # Get just the filename without path and without extension
    basename = os.path.basename(filepath)  # gets 'example2.raw'
    name_without_ext = os.path.splitext(basename)[0]  # gets 'example2'
    SAMPLES.append(name_without_ext)

# Print found samples for debugging
print(f"Found samples: {SAMPLES}")

# Checking to make sure that the raw file names match those in sample_annotation.txt
if not config.get("sample_annotation"):
    raise ValueError("Please provide sample_annotation file in config file")
sample_annotation = config["sample_annotation"]

# Read sample annotation file and get expected file names
try:
    sample_df = pd.read_csv(os.path.join(EXPERIMENT_DIR, sample_annotation), sep='\t')
    if 'file' not in sample_df.columns:
        raise ValueError("sample_annotation file must contain a 'file' column")
    expected_files = set(sample_df['file'].values)
except Exception as e:
    raise ValueError(f"Error reading sample annotation file: {str(e)}")

# Get actual raw file names (just the base names)
actual_files = set(SAMPLES)  # SAMPLES already contains just the base names

# Check for mismatches
missing_in_annotation = actual_files - expected_files
missing_in_raw = expected_files - actual_files

if missing_in_annotation:
    raise ValueError(f"Raw files found but not in sample annotation: {', '.join(missing_in_annotation)}")
if missing_in_raw:
    raise ValueError(f"Files in sample annotation but no matching raw files: {', '.join(missing_in_raw)}")

# Defining all of the output files
rule all:
    input:
        # Config files must be created first
        os.path.join(EXPERIMENT_DIR, "config/generate_diann_spectral_library.cfg"),
        os.path.join(EXPERIMENT_DIR, "config/run_diann.cfg"),
        os.path.join(EXPERIMENT_DIR,"logs/apptainer_checked"),
        "apptainer/conduitR.sif",
        "apptainer/diann2.1.0.sif",
        # Rest of the workflow outputs
        expand(os.path.join(EXPERIMENT_DIR, "output/00_database_resources/{file}"), 
               file=[
                   "00_database.fasta",
                   "00_proteome_ids.txt",
                   "01_taxonomy.txt",
                   "02_protein_info.txt",
                   "03_taxonomic_tree_of_database.pdf",
                   "04_database.predicted.speclib",
                   "README.md"
               ]),
        expand(os.path.join(EXPERIMENT_DIR, "output/00_database_resources/detected_protein_resources/{file}"),
               file=[
                   "00_detected_protein_info.txt",
                   "01_detected_protein.fasta",
                   "02_uniprot_annotated_protein_info.txt",
                   "03_go_annotations.txt",
                   "04_subcellular_locations.txt",
                   "05_kegg_annotations.txt"
               ]),
        expand(os.path.join(EXPERIMENT_DIR, "input/raw_files/{sample}.raw"), sample=SAMPLES),
        expand(os.path.join(EXPERIMENT_DIR, "output/01_diann_output/report.{suffix}.tsv"), suffix=["pr_matrix", "pg_matrix"]),
        expand(os.path.join(EXPERIMENT_DIR, "output/03_output_files/{level}_matrix.tsv"), 
               level=[
                   "domain", "kingdom", "phylum", "class", "order", 
                   "family", "genus", "species", "go", "go_taxa", "subcellular_locations",
                   "protein_group", "precursor", "peptide"
               ]),
        expand(os.path.join(EXPERIMENT_DIR, "output/03_output_files/{metric}.tsv"),
               metric=[
                   "database_taxonomy", "database_metrics",
                   "detected_protein_taxonomy", "detected_protein_metrics"
               ]),
        os.path.join(EXPERIMENT_DIR, "output/03_output_files/qf.rds"),
        os.path.join(EXPERIMENT_DIR, "output/03_output_files/conduit_output.rds")


# Handle DIANN spectral library config
rule setup_diann_spectral_library_config:
    output:
        sl_config_file = os.path.join(EXPERIMENT_DIR, "config/generate_diann_spectral_library.cfg")
    params:
        default_config = "config/generate_diann_spectral_library.cfg",
        sel_config = lambda wildcards: config.get("generate_diann_spectral_library_config")
    log: get_log_path("setup_diann_spectral_library_config")
    run:
        logger = setup_rule_logging(log[0])
        logger.info("Starting setup_diann_spectral_library_config rule")
        
        if params.sel_config == params.default_config:
            shell("cp {params.default_config} {output.sl_config_file}")
            logger.info("Copied default spectral library config to experiment directory")
        else:
            shell("cp {params.sel_config} {output.sl_config_file}")
            logger.info("Copied custom spectral library config to experiment directory")
            
        logger.info("setup_diann_spectral_library_config rule completed")

# Handle DIANN run config
rule setup_diann_run_config:
    output:
        diann_run_config_file = os.path.join(EXPERIMENT_DIR, "config/run_diann.cfg")
    params:
        default_config = "config/run_diann.cfg",
        sel_config = lambda wildcards: config.get("run_diann_config")
    log: get_log_path("setup_diann_run_config")
    run:
        logger = setup_rule_logging(log[0])
        logger.info("Starting setup_diann_run_config rule")
        
        if params.sel_config == params.default_config:
            shell("cp {params.default_config} {output.diann_run_config_file}")
            logger.info("Copied default run config to experiment directory")
        else:
            shell("cp {params.sel_config} {output.diann_run_config_file}")
            logger.info("Copied custom run config to experiment directory")
            
        logger.info("setup_diann_run_config rule completed")

################################################################################
# Making sure the system has all of the necessary software
################################################################################
# Checking if apptainer is installed, if using apptainer
rule check_apptainer:
    output:
        touch(os.path.join(EXPERIMENT_DIR,"logs/apptainer_checked"))
    log:
        get_log_path("check_apptainer")
    run:
        logger = setup_rule_logging(log[0])
        logger.info("Starting Apptainer check...")
        try:
            shell("apptainer --version")
            logger.info(f"Apptainer was found sucessfully")
        except Exception as e:
            logger.warning("Apptainer is not installed or not found in PATH.")
            logger.warning("Please install Apptainer before running this workflow if apptainer is desired.")
            logger.warning("Installation instructions: https://apptainer.org/docs/admin/main/installation.html")
        logger.info("check_apptainer rule completed successfully.")


# Build conduitR apptainer 
rule build_conduitR_apptainer:
    input:
        "apptainer/conduitR.def"
    output:
        "apptainer/conduitR.sif"
    log: get_log_path("build_conduitR_apptainer")
    run:
        logger = setup_rule_logging(log[0])
        logger.info("Starting build_conduitR_apptainer rule")
        try:
            logger.info("Attempting to build conduitR apptainer...")
            shell("apptainer build {output} {input}")
            logger.info("conduitR apptainer successfully built")
        except Exception as e:
            logger.warning(f"Unable to build conduitR apptainer: {str(e)}")
            logger.warning("If apptainer is desired, please ensure it is properly installed.")
            # Add a blank .sif file to the output directory to allow the workflow to continue to run
            shell("touch {output}")
        logger.info("build_conduitR_apptainer rule completed")

# Build diann apptainer 
rule build_diann_apptainer:
    input:
        "apptainer/diann2.1.0.def"
    output:
        "apptainer/diann2.1.0.sif"
    log: get_log_path("build_diann_apptainer")
    run:
        logger = setup_rule_logging(log[0])
        try:
            logger.info("Attempting to build diann2.1.0 apptainer...")
            shell("apptainer build {output} {input}")
            logger.info("diann2.1.0 apptainer successfully built")
        except Exception as e:
            logger.warning(f"Unable to build diann2.1.0 apptainer: {str(e)}")
            logger.warning("If apptainer is desired, please ensure it is properly installed.")
            # Add a blank .sif file to the output directory to allow the workflow to continue to run
            shell("touch {output}")
        logger.info("build_diann_apptainer rule completed")
            
################################################################################
# 00 Getting Database Resources (Defining the Search Space)
################################################################################        
# Downloading Fasta for each Organism ID
rule get_fasta:
    input:
      os.path.join(EXPERIMENT_DIR, "input/organisms.txt")
    output:
      # File containing proteome ids. Allows user to see what proteomes for each organism were downloaded.
      # Also contains NA values for NCBI taxon ids that are not in the Uniprot database.
      os.path.join(EXPERIMENT_DIR, "input/00_database_resources/00_proteome_ids.txt"),
      os.path.join(EXPERIMENT_DIR, "input/00_database_resources/00_database.fasta")
    log: get_log_path("get_fasta")
    container:
      "apptainer/conduitR.sif"
    script:
      "scripts/00_get_database_resources/00_get_fasta.R"
      
# Generating Full Taxonomy Information for each Organism ID
rule get_taxonomy:
    input:
        os.path.join(EXPERIMENT_DIR, "input/organisms.txt"),
        # File containing results of get fasta. Includes info about what taxa ids had cooresponding proteome ids.
        os.path.join(EXPERIMENT_DIR, "input/00_database_resources/00_proteome_ids.txt")
    output:
        os.path.join(EXPERIMENT_DIR, "input/00_database_resources/01_taxonomy.txt")
    log: get_log_path("get_taxonomy")
    container: "apptainer/conduitR.sif"
    script:
      "scripts/00_get_database_resources/01_get_taxonomy.R"
      
# Taking the data from fasta and taxonomy file, putting it in fasta.
rule get_protein_info_from_fasta:
    input:
      database_fasta=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/00_database.fasta"),
      taxonomy_txt=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/01_taxonomy.txt")
    output:
      os.path.join(EXPERIMENT_DIR, "input/00_database_resources/02_protein_info.txt")
    log: get_log_path("get_protein_info_from_fasta")
    container: "apptainer/conduitR.sif"
    script:
      "scripts/00_get_database_resources/02_get_protein_info_from_fasta.R"
      
# Plotting a taxonomic tree containing the taxonomy used in experiment
rule plot_taxonomic_tree:
    input:os.path.join(EXPERIMENT_DIR, "input/00_database_resources/01_taxonomy.txt")
    output:os.path.join(EXPERIMENT_DIR, "input/00_database_resources/03_taxonomic_tree_of_database.pdf")
    log: get_log_path("plot_taxonomic_tree")
    container:"apptainer/conduitR.sif"
    script:
     "scripts/00_get_database_resources/03_plot_taxonomic_tree.R"

# Creating a DataBase ReadMe          
rule make_database_resources_readme:
    input:os.path.join(EXPERIMENT_DIR, "input/00_database_resources/02_protein_info.txt")
    output:
      md =os.path.join(EXPERIMENT_DIR, "input/00_database_resources/README.md"),
      html = os.path.join(EXPERIMENT_DIR, "input/00_database_resources/README.html")
    log: get_log_path("make_database_resources_readme")
    container: "apptainer/conduitR.sif"
    script:
     "scripts/00_get_database_resources/04_make_database_resources_readme.R"

#################################################################################
# Generating Spectral Library
#################################################################################
rule generate_diann_spectral_library:
    input:
        fasta = os.path.join(EXPERIMENT_DIR, "input/00_database_resources/00_database.fasta"),
        config_file = os.path.join(EXPERIMENT_DIR, "config/generate_diann_spectral_library.cfg")
    output:
        os.path.join(EXPERIMENT_DIR, "input/00_database_resources/04_database.predicted.speclib")
    log: get_log_path("generate_diann_spectral_library")
    container:
        "apptainer/diann2.1.0.sif"
    params:
        config_file = lambda wildcards, input: input.config_file
    shell:
        """
        diann --cfg {{params.config_file}} \
        --fasta {{input.fasta}} \
        --out-lib {{output}}
        """

################################################################################
# 01 Running DIANN
################################################################################
rule run_diann:
    input:
        raw_files = RAW_FILEPATHS,
        spectral_library = os.path.join(EXPERIMENT_DIR, "input/00_database_resources/04_database.predicted.speclib"),
        fasta = os.path.join(EXPERIMENT_DIR, "input/00_database_resources/00_database.fasta"),
        config_file = os.path.join(EXPERIMENT_DIR, "config/run_diann.cfg")
    output:
        out = os.path.join(EXPERIMENT_DIR, "output/01_diann_output/"),
        report_pr_matrix = os.path.join(EXPERIMENT_DIR, "output/01_diann_output/report.pr_matrix.tsv"),
        report_pg_matrix = os.path.join(EXPERIMENT_DIR, "output/01_diann_output/report.pg_matrix.tsv")
    log: get_log_path("run_diann")
    container:
        "apptainer/diann2.1.0.sif"
    params:
        config_file = lambda wildcards, input: input.config_file
    shell:
        """
        diann --cfg {{params.config_file}} \
        --fasta {{input.fasta}} \
        --out {{output.out}} \
        --f {{input.raw_files}} \
        --lib {{input.spectral_library}} \
        --threads {{threads}} --verbose 1 
        """

################################################################################
# 02 Annotating Detected Proteins 
################################################################################
# Extracting all proteins that were detected
rule extract_detected_proteins:
  input:
    protein_info_df=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/02_protein_info.txt"),
    protein_info_fasta =os.path.join(EXPERIMENT_DIR, "input/00_database_resources/00_database.fasta"),
    report_pg_matrix= os.path.join(EXPERIMENT_DIR, "output/01_diann_output/report.pg_matrix.tsv")
  output:
    detected_protein_info_df = os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/00_detected_protein_info.txt"),
    detected_protein_info_fasta = os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/01_detected_protein.fasta")
  log: get_log_path("extract_detected_proteins")
  container: "apptainer/conduitR.sif"
  script:
    "scripts/02_get_detected_proteins_annotation/00_extract_detected_proteins.R"
  
# Get detected protein information from Uniprot
rule get_annotations_from_uniprot:
  input:
    detected_protein_info = os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/00_detected_protein_info.txt")
  output:
    uniprot_annotated_protein_info = os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt")
  log: get_log_path("get_annotations_from_uniprot")
  container: "apptainer/conduitR.sif"
  script:
    "scripts/02_get_detected_proteins_annotation/01_get_annotations_from_uniprot.R"
    
# Extracting GO infromation from data frame
rule extract_go_info:
    input:
     uniprot_annotated_protein_info=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt")
    output:
     go_annotations=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/03_go_annotations.txt")
    log: get_log_path("extract_go_info")
    container: "apptainer/conduitR.sif"
    script:
     "scripts/02_get_detected_proteins_annotation/02_extract_go_info.R" 

# Extracting cellular location information
rule extract_cellular_location_info:
    input:
      uniprot_annotated_protein_info=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt")
    output:
      subcellular_locations=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/04_subcellular_locations.txt")
    log: get_log_path("extract_cellular_location_info")
    container: "apptainer/conduitR.sif"
    script:
     "scripts/02_get_detected_proteins_annotation/03_extract_subcellular_locations.R" 
     
# Getting all Kegg infromation 
rule get_kegg_info:
    input:
      uniprot_annotated_protein_info=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt")
    output:
       kegg_annotations=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/05_kegg_annotations.txt")
    log: get_log_path("get_kegg_info")
    container: "apptainer/conduitR.sif"
    script:
     "scripts/02_get_detected_proteins_annotation/04_get_kegg_info.R" 
     
# Getting all Kegg infromation 
rule extract_detected_taxonomy:
    input:
      uniprot_annotated_protein_info=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt")
    output:
       detected_taxonomy=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/06_detected_taxonomy.txt")
    log: get_log_path("extract_detected_taxonomy")
    container: "apptainer/conduitR.sif"
    script:
     "scripts/02_get_detected_proteins_annotation/05_extract_detected_taxonomy.R"      
     
# rule get_cazyme_info:
#     input:
#      uniprot_annotated_protein_info="user_input/00_database_resources/04_detected_protein_info.txt"
#     output:
#      go_annotations="user_input/00_database_resources/05_go_annotations.txt",
#      subcellular_locations="user_input/00_database_resources/06_subcellular_locations.txt"
#     script:
#      "scripts/00_get_database_resources/03_extract_go_and_cellular_location_info.R" 
#      
# rule get_pfam_info:
#     input:
#      uniprot_annotated_protein_info="user_input/00_database_resources/04_detected_protein_info.txt"
#     output:
#      go_annotations="user_input/00_database_resources/05_go_annotations.txt",
#      subcellular_locations="user_input/00_database_resources/06_subcellular_locations.txt"
#     script:
#      "scripts/00_get_database_resources/03_extract_go_and_cellular_location_info.R" 
#      
# rule interpro_info:
#     input:
#      uniprot_annotated_protein_info="user_input/00_database_resources/04_detected_protein_info.txt"
#     output:
#      go_annotations="user_input/00_database_resources/05_go_annotations.txt",
#      subcellular_locations="user_input/00_database_resources/06_subcellular_locations.txt"
#     script:
#      "scripts/00_get_database_resources/03_extract_go_and_cellular_location_info.R" 

# # Making detected protein Resources ReadME
# rule make_detected_protein_resources_readme:
#     input:
#       detected_protein_info="user_input/00_database_resources/detected_protein_resources/00_detected_protein_info.txt"
#     output:
#        kegg_annotations="user_input/00_database_resources/detected_protein_resources/05_kegg_annotations.txt"
#     script:
#      "scripts/00_get_database_resources/04_get_kegg_info.R" 
#     

################################################################################
# 03 processing matrices
################################################################################
rule process_taxonomic_matrices:
  input:
    report_pr_matrix=os.path.join(EXPERIMENT_DIR, "output/01_diann_output/report.pr_matrix.tsv"),
    protein_info=os.path.join(EXPERIMENT_DIR, "input/00_database_resources/02_protein_info.txt")
  output:
    domain_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/domain_matrix.tsv"),
    kingdom_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/kingdom_matrix.tsv"),
    phylum_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/phylum_matrix.tsv"),
    class_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/class_matrix.tsv"),
    order_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/order_matrix.tsv"),
    family_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/family_matrix.tsv"),
    genus_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/genus_matrix.tsv"),
    species_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/species_matrix.tsv")
  log: get_log_path("process_taxonomic_matrices")
  container: "apptainer/conduitR.sif"
  script: "scripts/03_processing_matrices/00_processing_taxonomic_matrices.R"
    
rule process_go_matrix:
  input:
    report_pg_matrix=os.path.join(EXPERIMENT_DIR, "output/01_diann_output/report.pg_matrix.tsv"),
    go_annotations = os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/03_go_annotations.txt")
  output: 
    go_annotations_matrix = os.path.join(EXPERIMENT_DIR, "output/03_output_files/go_matrix.tsv"),
    go_annotations_taxa_matrix = os.path.join(EXPERIMENT_DIR, "output/03_output_files/go_taxa_matrix.tsv")
  log: get_log_path("process_go_matrix")
  container: "apptainer/conduitR.sif"
  script: "scripts/03_processing_matrices/01_process_go_matrix.R"
  
rule process_subcellular_locations_matrix:
  input:
    report_pg_matrix=os.path.join(EXPERIMENT_DIR, "output/01_diann_output/report.pg_matrix.tsv"),
    subcellular_locations = os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/04_subcellular_locations.txt")
  output: 
    subcellular_locations_matrix = os.path.join(EXPERIMENT_DIR, "output/03_output_files/subcellular_locations_matrix.tsv")
  log: get_log_path("process_subcellular_locations_matrix")
  container: "apptainer/conduitR.sif"
  script:"scripts/03_processing_matrices/02_process_subcellular_locations_matrix.R"
  
# rule process_kegg_matrices:
#   input:
#     report_pg_matrix="output/03_diann_output/report.pg_matrix.tsv"
#   output: 
#     kegg_matrix = "output/05_output_files/subcellular_locations_matrix.tsv"
#   script:"scripts/06_processing_matrices/03_process_kegg_matrix.R"  
# 

rule process_diann_matrices:
  input:
    report_pg_matrix=os.path.join(EXPERIMENT_DIR, "output/01_diann_output/report.pg_matrix.tsv"),
    report_pr_matrix =os.path.join(EXPERIMENT_DIR, "output/01_diann_output/report.pr_matrix.tsv")
  output: 
    protein_group_matrix = os.path.join(EXPERIMENT_DIR, "output/03_output_files/protein_group_matrix.tsv"),
    precursor_matrix = os.path.join(EXPERIMENT_DIR, "output/03_output_files/precursor_matrix.tsv"),
    peptide_matrix = os.path.join(EXPERIMENT_DIR, "output/03_output_files/peptide_matrix.tsv")
  log: get_log_path("process_diann_matrices")
  container: "apptainer/conduitR.sif"
  script:"scripts/03_processing_matrices/04_process_diann_matrices.R"
    
    
rule move_database_resources:
    input:
        expand(os.path.join(EXPERIMENT_DIR, "input/00_database_resources/{filename}"), 
               filename=[
                   "00_database.fasta",
                   "00_proteome_ids.txt",
                   "01_taxonomy.txt",
                   "02_protein_info.txt",
                   "03_taxonomic_tree_of_database.pdf",
                   "04_database.predicted.speclib",
                   "README.md",
                   "README.html"
               ]),
        expand(os.path.join(EXPERIMENT_DIR, "input/00_database_resources/detected_protein_resources/{filename}"), 
               filename=[
                   "00_detected_protein_info.txt",
                   "01_detected_protein.fasta",
                   "02_uniprot_annotated_protein_info.txt",
                   "03_go_annotations.txt",
                   "04_subcellular_locations.txt",
                   "05_kegg_annotations.txt"
               ]),
    output:
        expand(os.path.join(EXPERIMENT_DIR, "output/00_database_resources/{filename}"), 
               filename=[
                   "00_database.fasta",
                   "00_proteome_ids.txt",
                   "01_taxonomy.txt",
                   "02_protein_info.txt",
                   "03_taxonomic_tree_of_database.pdf",
                   "04_database.predicted.speclib",
                   "README.md",
                   "README.html"
               ]),
        expand(os.path.join(EXPERIMENT_DIR, "output/00_database_resources/detected_protein_resources/{filename}"), 
               filename=[
                   "00_detected_protein_info.txt",
                   "01_detected_protein.fasta",
                   "02_uniprot_annotated_protein_info.txt",
                   "03_go_annotations.txt",
                   "04_subcellular_locations.txt",
                   "05_kegg_annotations.txt"
               ]),
    log: get_log_path("move_database_resources")
    shell:
        """
        mkdir -p output/00_database_resources/detected_protein_resources
        cp -u -r input/00_database_resources/* output/00_database_resources/
        cp -u -r input/00_database_resources/detected_protein_resources/* output/00_database_resources/detected_protein_resources/
        """
################################################################################
# 04 Ingestion into R
################################################################################
# Prepare a QFeatures object with all the data
rule prepare_qf:
  input:
    annotation = os.path.join(EXPERIMENT_DIR, "input/sample_annotation.txt"),
    precursor_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/precursor_matrix.tsv"),
    peptide_matrix = os.path.join(EXPERIMENT_DIR, "output/03_output_files/peptide_matrix.tsv"),
    protein_group_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/protein_group_matrix.tsv"),
    domain_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/domain_matrix.tsv"),
    kingdom_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/kingdom_matrix.tsv"),
    phylum_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/phylum_matrix.tsv"),
    class_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/class_matrix.tsv"),
    order_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/order_matrix.tsv"),
    family_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/family_matrix.tsv"),
    genus_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/genus_matrix.tsv"),
    species_matrix=os.path.join(EXPERIMENT_DIR, "output/03_output_files/species_matrix.tsv"),
    go_matrix = os.path.join(EXPERIMENT_DIR, "output/03_output_files/go_matrix.tsv"),
    go_taxa_matrix = os.path.join(EXPERIMENT_DIR, "output/03_output_files/go_taxa_matrix.tsv"),
    subcellular_locations_matrix = os.path.join(EXPERIMENT_DIR, "output/03_output_files/subcellular_locations_matrix.tsv")
  output:
    qf= os.path.join(EXPERIMENT_DIR, "output/03_output_files/qf.rds")
  log: get_log_path("prepare_qf")
  container: "apptainer/conduitR.sif"
  script:
    "scripts/04_ingestion_into_R/00_prepare_qf.R"
    
# prepare metrics    
rule extract_metrics:
  input:
    protein_info = os.path.join(EXPERIMENT_DIR, "output/00_database_resources/02_protein_info.txt"),
    detected_protein_info = os.path.join(EXPERIMENT_DIR, "output/00_database_resources/detected_protein_resources/00_detected_protein_info.txt")
  output:
    database_taxonomy=os.path.join(EXPERIMENT_DIR, "output/03_output_files/database_taxonomy.tsv"),
    database_metrics=os.path.join(EXPERIMENT_DIR, "output/03_output_files/database_metrics.tsv"),
    detected_protein_taxonomy = os.path.join(EXPERIMENT_DIR, "output/03_output_files/detected_protein_taxonomy.tsv"),
    detected_protein_metrics = os.path.join(EXPERIMENT_DIR, "output/03_output_files/detected_protein_metrics.tsv"),
    combined_metrics = os.path.join(EXPERIMENT_DIR, "output/03_output_files/combined_metrics.tsv")
  log: get_log_path("extract_metrics")
  container: "apptainer/conduitR.sif"
  script:
    "scripts/04_ingestion_into_R/01_extract_metrics.R"
    
# Prepare a conduit object
# Experiment metadata = date of experiment
# Database taxonomy = n proteins belonging to each taxa in database
# database metrics = n_proteins, n_organism_types, domain, kingdom, phylum, class, order, family, genus,species
# detected protein taxonomy = n proteins detected belonging to each taxa
# detected protein metrics = n_protein_groups,n_uniquely_ided_proteins,n_precursors,n_peptides n_organism_types, domain, kingdom, phylum, class, order, family, genus,species
rule create_conduit:
  input:
    qf = os.path.join(EXPERIMENT_DIR, "output/03_output_files/qf.rds"),
    database_taxonomy = os.path.join(EXPERIMENT_DIR, "output/03_output_files/database_taxonomy.tsv"),
    database_metrics = os.path.join(EXPERIMENT_DIR, "output/03_output_files/database_metrics.tsv/"),
    detected_protein_taxonomy = os.path.join(EXPERIMENT_DIR, "output/03_output_files/detected_protein_taxonomy.tsv"),
    detected_protein_metrics = os.path.join(EXPERIMENT_DIR, "output/03_output_files/detected_protein_metrics.tsv")
  output:
    conduit_obj = os.path.join(EXPERIMENT_DIR, "output/03_output_files/conduit_output.rds")
  log: get_log_path("create_conduit")
  container: "apptainer/conduitR.sif"
    # From this point- the conduit object contains all the stuff you need. 
    # Load that into Conduit GUI.
  script:
    "scripts/04_ingestion_into_R/02_create_conduit_obj.R"