import glob
import os
EXPERIMENT_DIR = os.path.join("experiments",config["experiment"])
RAW_FILEPATHS = glob.glob(os.path.join(EXPERIMENT_DIR, "input/raw_files/*.raw"))
################################################################################
# Building the Apptainer Image for Umgap
################################################################################
# This is needed to run the script generating the . See https://github.com/unipept/unipept-database/issues/75
rule build_umgap_apptainer:
    input:
        apptainer/umgap.def
    output:
        apptainer/umgap.sif
    log: os.path.join(EXPERIMENT_DIR,"logs/search_space/build_umgap_apptainer.log")
    shell:
        """
        apptainer build umgap.sif umgap.def
        """
################################################################################
# Generating the Sequence Index
################################################################################
# This is needed to generate the file containing all peptides in TREMBL and SWISSPROT
# And their LCAs. See https://github.com/unipept/unipept-database/issues/75
rule build_sequence_index:
    output:
        resources/databases/sequences.tsv.lz4
    container: "apptainer/umgap.sif"
    log: os.path.join(EXPERIMENT_DIR,"logs/search_space/generate_first_pass_proteotyping_db.log")
    shell:
        """
        modules/search_space/proteotyping/scripts/unipept-database/scripts/generate_umgap_tables.sh tryptic \
        --output-dir resources/databases\
        --database-sources  'swissprot,trembl' \
        --temp-dir resources/databases/temp/ \
        --min-peptide-length 5 \
        --max-peptide-length 50 
        >> {log} 2>&1
        """
################################################################################
# Determining Version of Uniprotkb used to generate the sequence index. 
################################################################################
# The user should be able to know what version of uniprotkb was used to generate
# the sequence index. Furthermore, the user should be identified if the version
# of the sequence index is up to date with the current version of uniprotkb or not. 
rule timestamp_sequence_index:
    output:
        resources/databases/sequences.tsv.lz4.timestamp
    log: os.path.join(EXPERIMENT_DIR,"logs/search_space/timestamp_sequence_index.log")
    shell:
        """
        # Determine the version of uniprotkb used to generate the sequence index
        
        # Determine the current version of uniprotkb

        # Determine if the version used is current or not

        """
################################################################################
# Generating the Proteotyping Fasta DB
################################################################################
# From the file containing the peptides and their LCAS, we can generate a fasta file
# only containing proteotypic peptides for a given taxonomy. We can subsequently use 
# this as a first pass database to identify taxa that are likely present in the experiment.
rule build_proteotyping_fasta_db:
    input:
        resources/databases/sequences.tsv.lz4
    output:
        resources/databases/first_pass_proteotyping.fasta
    container: "apptainer/conduitR.sif"
    log: os.path.join(EXPERIMENT_DIR,"logs/search_space/generate_first_pass_proteotyping_db.log")
    script: 
        "scripts/generate_proteotyping_fasta_db.R"


################################################################################
# Generating the First Pass Spectral Library
################################################################################
# We can generate a spectral library from the proteotyping fasta db to use with
# Diann. This can then be used to perform the first pass search identifying 
# species that are likely present in the experiment. 

rule generate_first_pass_spectral_library:
    input: 
        fasta = "resources/databases/first_pass_database.fasta",
        config_file = "config/proteotyping_firstpass_diann_spectral_library.cfg"
    output: "resources/databases/first_pass_database.predicted.speclib"
    container: "apptainer/diann.sif"
    log: os.path.join(EXPERIMENT_DIR,"logs/search_space/generate_first_pass_spectral_library.log")
    threads: workflow.cores
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --threads {threads} \
        --out-lib resources/databases/first_pass_proteotyping_db >> {log} 2>&1
        """

################################################################################
# Performing the First Pass Search
################################################################################
# Searching our first pass spectral library with Diann
rule perform_first_pass_search:
    input:
        raw_files_dir = os.path.join(EXPERIMENT_DIR,"input/raw_files"),
        spectral_library = "resources/databases/first_pass_proteotyping_db.predicted.speclib",
        fasta = "resources/databases/first_pass_proteotyping_db.predicted.speclib",
        config_file = os.path.join(EXPERIMENT_DIR,"config/run_diann.cfg")
    output:
        report_pr_matrix =os.path.join(EXPERIMENT_DIR,"resources/database/first_pass_diann.pr_matrix.tsv"),
        report_pg_matrix = "resources/database/first_pass_diann.pg_matrix.tsv",
    log: os.path.join(EXPERIMENT_DIR,"logs/diann/run_diann.log")
    container:
        "apptainer/diann2.1.0.sif"
    threads: workflow.cores 
    shell:
        """
        diann --cfg {input.config_file} \
        --fasta {input.fasta} \
        --out  resources/databases/first_pass_diann \
        --dir {input.raw_files_dir} \
        --lib {input.spectral_library} \
        --threads {threads} --verbose 1 >> {log} 2>&1
        """

################################################################################
# Scoring the likelihood of species being present based on the first pass search results
################################################################################
# We aim to score the likelihood that a species is truly present based on the number of detected
# species-specific proteotypic peptides. To do this, we calculate the fraction of that species’ 
# known proteotypic peptides that were observed in the experiment, and compare this to the number 
# we would expect to appear by chance under a 1% FDR.

rule score_species_likelihood:
    input:
        report_pr_matrix = "resources/database/first_pass_diann.pr_matrix.tsv"
    output:
        "resources/database/first_pass_diann.pr_matrix.tsv"
    shell:

        # Calculate the expected false discoveries per species
        # expected false discoveries per species ≈ 
        # total hits × FDR × (species-specific peptide count / total DB peptides)
        
        will need:
	•	P_detected = number of peptides observed for the species
	•	P_species = number of proteotypic peptides for the species in the database
	•	P_total = total number of proteotypic peptides in the full database
	•	N_total_detected = total number of peptides detected (passing FDR threshold)
	•	FDR = your false discovery rate (e.g. 0.01)

    E_species = N_total_detected * FDR * (P_species / P_total)
    score = P_detected / (E_species + ε)

################################################################################
# Generating a database based on first pass search results
################################################################################
# We can generate a database based on the first pass search results. This will
# be used to identify taxa that are likely present in the experiment.

rule generate_first_pass_identified_taxa:
    input:
        species_scores = "resources/database/first_pass_diann.pr_matrix.tsv"
    output:
        os.path.join(EXPERIMENT_DIR,"input/database_resources/database.fasta")
    shell:
        """
        # Define threshold for species to be considered present. 
        
        # Filter the species scores to only include species with a score greater than 0.5

        # Get taxaproteome ids corresponding to these species

        # Write ncbi_taxonomy_ids to a file. 
        ncbi_taxonomy_ids  
        """

# After this, we can just plug the first pass identified taxa into the the ncbi_taxonomy workflow to perform the 
# main search!


   