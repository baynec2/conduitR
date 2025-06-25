import glob
import os
EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/raw_files/*.raw"))
#################################################################################
# Generating Spectral Library
#################################################################################
rule generate_diann_spectral_library:
    input:
        fasta = os.path.join(EXPERIMENT_DIR,"input/database_resources/database.fasta"),
        config_file = os.path.join(EXPERIMENT_DIR,"config/generate_diann_spectral_library.cfg")
    output:
        os.path.join(EXPERIMENT_DIR,"input/database_resources/database.predicted.speclib")
    log: os.path.join(EXPERIMENT_DIR,"logs/diann/generate_diann_spectral_library.log")
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
# Running DIANN
################################################################################
rule run_diann:
    input:
        raw_files = RAW_FILEPATHS,
        spectral_library = os.path.join(EXPERIMENT_DIR,"input/database_resources/database.predicted.speclib"),
        fasta = os.path.join(EXPERIMENT_DIR,"input/database_resources/database.fasta"),
        config_file = os.path.join(EXPERIMENT_DIR,"config/run_diann.cfg")
    output:
        out = os.path.join(EXPERIMENT_DIR,"output/diann_output/"),
        report_pr_matrix = os.path.join(EXPERIMENT_DIR,"output/diann_output/report.pr_matrix.tsv"),
        report_pg_matrix = os.path.join(EXPERIMENT_DIR,"output/diann_output/report.pg_matrix.tsv")
    log: os.path.join(EXPERIMENT_DIR,"logs/diann/run_diann.log")
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
# Extracting Detected Proteins
################################################################################
rule extract_detected_proteins:
  input:
    protein_info_df=os.path.join(EXPERIMENT_DIR,"input/database_resources/protein_info.txt"),
    protein_info_fasta =os.path.join(EXPERIMENT_DIR,"input/database_resources/database.fasta"),
    report_pg_matrix=os.path.join(EXPERIMENT_DIR,"output/diann_output/report.pg_matrix.tsv")
  output:
    detected_protein_info_df = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/detected_protein_info.txt"),
    detected_protein_info_fasta = os.path.join(EXPERIMENT_DIR,"input/database_resources/detected_protein_resources/detected_protein.fasta")
  log: os.path.join(EXPERIMENT_DIR,"logs/diann/extract_detected_proteins.log")
  container: "apptainer/conduitR.sif"
  script:
    "scripts/extract_detected_proteins.R"