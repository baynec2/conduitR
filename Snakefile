# Importing os
import os
# Defining the .raw samples that are inputs to the workflow. 
SAMPLES = [os.path.splitext(f)[0] for f in os.listdir("user_input/raw_files/") if f.endswith(".raw")]
# Defining all of the output files
rule all:
    input:
        expand("output/00_database_resources/{file}", 
               file=[
                   "00_database.fasta",
                   "01_taxonomy.txt",
                   "02_protein_info.txt",
                   "03_taxonomic_tree_of_database.pdf",
                   "README.md"
               ]),
        expand("output/00_database_resources/detected_protein_resources/{file}",
               file=[
                   "00_detected_protein_info.txt",
                   "01_detected_protein.fasta",
                   "02_uniprot_annotated_protein_info.txt",
                   "03_go_annotations.txt",
                   "04_subcellular_locations.txt",
                   "05_kegg_annotations.txt"
               ]),
        expand("user_input/raw_files/{sample}.raw", sample=SAMPLES),
        "output/02_speclib_gen/database.predicted.speclib",
        expand("output/03_diann_output/report.{suffix}.tsv", suffix=["pr_matrix", "pg_matrix"]),
        expand("output/05_output_files/{level}_matrix.tsv", 
               level=[
                   "superkingdom", "kingdom", "phylum", "class", "order", 
                   "family", "genus", "species", "go", "go_taxa", "subcellular_locations",
                   "protein_group", "precursor", "peptide"
               ]),
        expand("output/05_output_files/{metric}.tsv",
               metric=[
                   "database_taxonomy", "database_metrics",
                   "detected_protein_taxonomy", "detected_protein_metrics"
               ]),
        "output/05_output_files/qf.rds",
        "output/05_output_files/conduit_output.rds"
################################################################################
# Making sure the system has all of the necessary software
################################################################################
# Need apptainer
rule check_apptainer:
    output:
        "log/log.txt"
    shell:
        """
        if ! command -v apptainer &> /dev/null; then
            echo "apptainer not found. Attempting to install apptainer..." >&2
            sudo apt-get update -y && sudo apt-get install -y apptainer
            if ! command -v apptainer &> /dev/null; then
                echo "Apptainer installation failed. Please install apptainer manually." >&2
                exit 1
            fi
        fi
        apptainer --version > {output}
        """
# Build diann apptainer 
rule build_conduitR_apptainer:
    input:
        "apptainer/conduitR.def"
    output:
        log="log/log.txt",
        sif="apptainer/conduitR.sif"
    shell:
        """
        apptainer build {output.sif} {input}
        echo "conduitR apptainer successfully built" > {output.log}
        """
# We no longer need this since Diann2.1 supports .raw files naitively.
# thermorawfileparser container can be downloaded from 
#https://quay.io/repository/biocontainers/thermorawfileparser?tab=tags
# rule build_thermorawfileparser_apptainer:
#     output:
#         log="log/log.txt",
#         sif="apptainer/thermorawfileparser.sif"
#     shell:
#         """
#         apptainer build {output.sif} "docker://quay.io/biocontainers/thermorawfileparser:1.4.5--h05cac1d_1"
#         echo "ThermoRawFileParser Apptainer successfully built" > {output.log}
#         """

# Build diann apptainer 
rule build_diann_apptainer:
    output:
        "log/log.txt",
        "apptainer/diann2.1.0.sif"
    shell:
        """
        apptainer build apptainer/diann2.1.0.sif apptainer/diann2.1.0.def
        echo "diann2.1.0 apptainer successfully built" > {output}
        """
################################################################################
# 00 Getting Database Resources (Defining the Search Space)
################################################################################        
# Downloading Fasta for each Organism ID
rule get_fasta:
    input:
      "user_input/organisms.txt"
    output:
      "user_input/00_database_resources/00_database.fasta"
    container:
      "apptainer/conduitR.sif"
    script:
      "scripts/00_get_database_resources/00_get_fasta.R"
      
# Generating Full Taxonomy Information for each Organism ID
rule get_taxonomy:
    input:
        "user_input/organisms.txt"
    output:
        "user_input/00_database_resources/01_taxonomy.txt"
    container: "apptainer/conduitR.sif"
    script:
      "scripts/00_get_database_resources/01_get_taxonomy.R"
      
# Taking the data from fasta and taxonomy file, putting it in fasta.
rule get_protein_info_from_fasta:
    input:
      database_fasta="user_input/00_database_resources/00_database.fasta",
      taxonomy_txt="user_input/00_database_resources/01_taxonomy.txt"
    output:
      "user_input/00_database_resources/02_protein_info.txt"
    container: "apptainer/conduitR.sif"
    script:
      "scripts/00_get_database_resources/02_get_protein_info_from_fasta.R"
      
# Plotting a taxonomic tree containing the taxonomy used in experiment
rule plot_taxonomic_tree:
    input:"user_input/00_database_resources/01_taxonomy.txt"
    output:"user_input/00_database_resources/03_taxonomic_tree_of_database.pdf"
    container:"apptainer/conduitR.sif"
    script:
     "scripts/00_get_database_resources/03_plot_taxonomic_tree.R"

# Creating a DataBase ReadMe          
rule make_database_resources_readme:
    input:"user_input/00_database_resources/02_protein_info.txt"
    output:
      md ="user_input/00_database_resources/README.md",
      html = "user_input/00_database_resources/README.html"
    container: "apptainer/conduitR.sif"
    script:
     "scripts/00_get_database_resources/04_make_database_resources_readme.R"
# ################################################################################
# # 01 Converting .raw files to MzmL
# ################################################################################
# rule convert_raw_to_mzml:
#     input:
#         expand("user_input/raw_files/{sample}.raw", sample=SAMPLES)
#     output:
#         expand("output/01_file_conversion/{sample}.mzML", sample=SAMPLES)
#     container: "apptainer/thermorawfileparser.sif"
#     shell:
#         """
#         # Run ThermoRawFileParser using container-mapped paths
#         thermorawfileparser --input {input} --output {output}
#         """

#################################################################################
# 02 Generating Spectral Library
#################################################################################
rule generate_diann_spectral_library:
    input:
        fasta = "user_input/00_database_resources/00_database.fasta"
    output:
        "output/02_speclib_gen/database.predicted.speclib"
    container:
        "apptainer/diann2.1.0.sif"
    params:
        config_file = "config/generate_diann_spectral_library.cfg"
    shell:
        """
        diann --cfg {params.config_file} \
        --fasta {input.fasta} \
        --out-lib {output}
        """  
################################################################################
# 03 Running DIANN
################################################################################
rule run_diann:
    input:
        raw_files = expand("user_input/raw_files/{sample}.raw", sample=SAMPLES),
        spectral_library = "output/02_speclib_gen/database.predicted.speclib",
        fasta = "user_input/00_database_resources/00_database.fasta"
    output:
        out = "output/03_diann_output/",
        report_pr_matrix = "output/03_diann_output/report.pr_matrix.tsv",
        report_pg_matrix = "output/03_diann_output/report.pg_matrix.tsv"
    container:
        "apptainer/diann2.1.0.sif"
    params:
        config_file = "config/run_diann.cfg"  # Config file for DIA-NN
    shell:
        """
        diann --cfg {params.config_file} \
        --fasta {input.fasta} \
        --out "output/03_diann_output" \
        --f {input.raw_files} \
        --lib {input.spectral_library} \
        --threads {threads} --verbose 1 
        """
################################################################################
# 04 Cleaning up MZML Files.
################################################################################
# Files from the Astral take a ton of space. It is not feasible to store what is
# essentially two copies of the information. As such, we will deleted the MZML
# files. Ideally, future versions of diann will be able to operate directly on
# .raw files on linux, circumventing the need for the file conversion. 

# rule remove_mzml_files:
#     input:
#         "output/01_file_conversion/mzml_files/{sample}.mzML"  # Ensure input files are specified if necessary
#     output:
#         "log/log.txt"
#     shell:
#         """
#         rm -r output/01_file_conversion/mzml_files/
#         echo "MZML files removed" > {output}
#         """
        
################################################################################
# 05 Annotating Detected Proteins 
################################################################################
# Extracting all proteins that were detected
rule extract_detected_proteins:
  input:
    protein_info_df="user_input/00_database_resources/02_protein_info.txt",
    protein_info_fasta ="user_input/00_database_resources/00_database.fasta",
    report_pg_matrix= "output/03_diann_output/report.pg_matrix.tsv"
  output:
    detected_protein_info_df = "user_input/00_database_resources/detected_protein_resources/00_detected_protein_info.txt",
    detected_protein_info_fasta = "user_input/00_database_resources/detected_protein_resources/01_detected_protein.fasta"
  container: "apptainer/conduitR.sif"
  script:
    "scripts/05_get_detected_proteins_annotation/00_extract_detected_proteins.R"
  
# Get detected protein information from Uniprot
rule get_annotations_from_uniprot:
  input:
    detected_protein_info = "user_input/00_database_resources/detected_protein_resources/00_detected_protein_info.txt"
  output:
    uniprot_annotated_protein_info = "user_input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt"
  container: "apptainer/conduitR.sif"
  script:
    "scripts/05_get_detected_proteins_annotation/01_get_annotations_from_uniprot.R"
    
# Extracting GO infromation from data frame
rule extract_go_info:
    input:
     uniprot_annotated_protein_info="user_input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt"
    output:
     go_annotations="user_input/00_database_resources/detected_protein_resources/03_go_annotations.txt"
    container: "apptainer/conduitR.sif"
    script:
     "scripts/05_get_detected_proteins_annotation/02_extract_go_info.R" 

# Extracting cellular location information
rule extract_cellular_location_info:
    input:
      uniprot_annotated_protein_info="user_input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt"
    output:
      subcellular_locations="user_input/00_database_resources/detected_protein_resources/04_subcellular_locations.txt"
    container: "apptainer/conduitR.sif"
    script:
     "scripts/05_get_detected_proteins_annotation/03_extract_subcellular_locations.R" 
     
# Getting all Kegg infromation 
rule get_kegg_info:
    input:
      uniprot_annotated_protein_info="user_input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt"
    output:
       kegg_annotations="user_input/00_database_resources/detected_protein_resources/05_kegg_annotations.txt"
    container: "apptainer/conduitR.sif"
    script:
     "scripts/05_get_detected_proteins_annotation/04_get_kegg_info.R" 
     
# Getting all Kegg infromation 
rule extract_detected_taxonomy:
    input:
      uniprot_annotated_protein_info="user_input/00_database_resources/detected_protein_resources/02_uniprot_annotated_protein_info.txt"
    output:
       detected_taxonomy="user_input/00_database_resources/detected_protein_resources/06_detected_taxonomy.txt"
    container: "apptainer/conduitR.sif"
    script:
     "scripts/05_get_detected_proteins_annotation/05_extract_detected_taxonomy.R"      
     
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
# 06 processing matrices
################################################################################
rule process_taxonomic_matrices:
  input:
    report_pr_matrix="output/03_diann_output/report.pr_matrix.tsv",
    protein_info="user_input/00_database_resources/02_protein_info.txt"
  output:
    superkingdom_matrix="output/05_output_files/superkingdom_matrix.tsv",
    kingdom_matrix="output/05_output_files/kingdom_matrix.tsv",
    phylum_matrix="output/05_output_files/phylum_matrix.tsv",
    class_matrix="output/05_output_files/class_matrix.tsv",
    order_matrix="output/05_output_files/order_matrix.tsv",
    family_matrix="output/05_output_files/family_matrix.tsv",
    genus_matrix="output/05_output_files/genus_matrix.tsv",
    species_matrix="output/05_output_files/species_matrix.tsv"
  container: "apptainer/conduitR.sif"
  script: "scripts/06_processing_matrices/00_processing_taxonomic_matrices.R"
    
rule process_go_matrix:
  input:
    report_pg_matrix="output/03_diann_output/report.pg_matrix.tsv",
    go_annotations = "user_input/00_database_resources/detected_protein_resources/03_go_annotations.txt"
  output: 
    go_annotations_matrix = "output/05_output_files/go_matrix.tsv",
    go_annotations_taxa_matrix = "output/05_output_files/go_taxa_matrix.tsv"
  container: "apptainer/conduitR.sif"
  script: "scripts/06_processing_matrices/01_process_go_matrix.R"
  
rule process_subcellular_locations_matrix:
  input:
    report_pg_matrix="output/03_diann_output/report.pg_matrix.tsv",
    subcellular_locations = "user_input/00_database_resources/detected_protein_resources/04_subcellular_locations.txt"
  output: 
    subcellular_locations_matrix = "output/05_output_files/subcellular_locations_matrix.tsv"
  container: "apptainer/conduitR.sif"
  script:"scripts/06_processing_matrices/02_process_subcellular_locations_matrix.R"
  
# rule process_kegg_matrices:
#   input:
#     report_pg_matrix="output/03_diann_output/report.pg_matrix.tsv"
#   output: 
#     kegg_matrix = "output/05_output_files/subcellular_locations_matrix.tsv"
#   script:"scripts/06_processing_matrices/03_process_kegg_matrix.R"  
# 

rule process_diann_matrices:
  input:
    report_pg_matrix="output/03_diann_output/report.pg_matrix.tsv",
    report_pr_matrix ="output/03_diann_output/report.pr_matrix.tsv"
  output: 
    protein_group_matrix = "output/05_output_files/protein_group_matrix.tsv",
    precursor_matrix = "output/05_output_files/precursor_matrix.tsv",
    peptide_matrix = "output/05_output_files/peptide_matrix.tsv"
  container: "apptainer/conduitR.sif"
  script:"scripts/06_processing_matrices/04_process_diann_matrices.R"
    
    
rule move_database_resources:
    input:
        expand("user_input/00_database_resources/{filename}", 
               filename=[
                   "00_database.fasta",
                   "01_taxonomy.txt",
                   "02_protein_info.txt",
                   "03_taxonomic_tree_of_database.pdf",
                   "README.md",
                   "README.html"
               ]),
        expand("user_input/00_database_resources/detected_protein_resources/{filename}", 
               filename=[
                   "00_detected_protein_info.txt",
                   "01_detected_protein.fasta",
                   "02_uniprot_annotated_protein_info.txt",
                   "03_go_annotations.txt",
                   "04_subcellular_locations.txt",
                   "05_kegg_annotations.txt"
               ]),
    output:
        expand("output/00_database_resources/{filename}", 
               filename=[
                   "00_database.fasta",
                   "01_taxonomy.txt",
                   "02_protein_info.txt",
                   "03_taxonomic_tree_of_database.pdf",
                   "README.md",
                   "README.html"
               ]),
        expand("output/00_database_resources/detected_protein_resources/{filename}", 
               filename=[
                   "00_detected_protein_info.txt",
                   "01_detected_protein.fasta",
                   "02_uniprot_annotated_protein_info.txt",
                   "03_go_annotations.txt",
                   "04_subcellular_locations.txt",
                   "05_kegg_annotations.txt"
               ]),
    shell:
        """
        mkdir -p output/00_database_resources/detected_protein_resources
        cp -u -r user_input/00_database_resources/* output/00_database_resources/
        cp -u -r user_input/00_database_resources/detected_protein_resources/* output/00_database_resources/detected_protein_resources/
        """
# rule cp_diann_output:
#   # Now that all of the database stuff is finished processing, let's move to output. 
#     input:
#       # Database Resources
#         report_pr_matrix="output/03_diann_output/report.pr_matrix.tsv",
#         report_pg_matrix = "output/03_diann_output/report.pg_matrix.tsv"
#       # Detected Protein Annotations
#     output:
#         precursor_matrix = "output/05_output_files/precursor_matrix.tsv",
#         protein_group_matrix = "output/05_output_files/protein_group_matrix.tsv"
#     run:
#         shell("cp output/03_diann_output/report.pr_matrix.tsv output/05_output_files/precursor_matrix.tsv")
#         shell("cp output/03_diann_output/report.pg_matrix.tsv output/05_output_files/protein_group_matrix.tsv")

################################################################################
# 07 Ingestion into R
################################################################################
# Prepare a QFeatures object with all the data
rule prepare_qf:
  input:
    annotation = "user_input/sample_annotation.txt",
    precursor_matrix="output/05_output_files/precursor_matrix.tsv",
    peptide_matrix = "output/05_output_files/peptide_matrix.tsv",
    protein_group_matrix="output/05_output_files/protein_group_matrix.tsv",
    superkingdom_matrix="output/05_output_files/superkingdom_matrix.tsv",
    kingdom_matrix="output/05_output_files/kingdom_matrix.tsv",
    phylum_matrix="output/05_output_files/phylum_matrix.tsv",
    class_matrix="output/05_output_files/class_matrix.tsv",
    order_matrix="output/05_output_files/order_matrix.tsv",
    family_matrix="output/05_output_files/family_matrix.tsv",
    genus_matrix="output/05_output_files/genus_matrix.tsv",
    species_matrix="output/05_output_files/species_matrix.tsv",
    go_matrix = "output/05_output_files/go_matrix.tsv",
    go_taxa_matrix = "output/05_output_files/go_taxa_matrix.tsv",
    subcellular_locations_matrix = "output/05_output_files/subcellular_locations_matrix.tsv"
  output:
    qf= "output/05_output_files/qf.rds"
  container: "apptainer/conduitR.sif"
  script:
    "scripts/07_ingestion_into_R/00_prepare_qf.R"
  
    
# prepare metrics    
rule extract_metrics:
  input:
    protein_info = "output/00_database_resources/02_protein_info.txt",
    detected_protein_info = "output/00_database_resources/detected_protein_resources/00_detected_protein_info.txt"
  output:
    database_taxonomy="output/05_output_files/database_taxonomy.tsv",
    database_metrics="output/05_output_files/database_metrics.tsv",
    detected_protein_taxonomy = "output/05_output_files/detected_protein_taxonomy.tsv",
    detected_protein_metrics = "output/05_output_files/detected_protein_metrics.tsv",
    combined_metrics = "output/05_output_files/combined_metrics.tsv"
  container: "apptainer/conduitR.sif"
  script:
    "scripts/07_ingestion_into_R/01_extract_metrics.R"
    
# Prepare a conduit object
# Experiment metadata = date of experiment
# Database taxonomy = n proteins belonging to each taxa in database
# database metrics = n_proteins, n_organism_types, superkingdom, kingdom, phylum, class, order, family, genus,species
# detected protein taxonomy = n proteins detected belonging to each taxa
# detected protein metrics = n_protein_groups,n_uniquely_ided_proteins,n_precursors,n_peptides n_organism_types, superkingdom, kingdom, phylum, class, order, family, genus,species
rule create_conduit:
  input:
    qf = "output/05_output_files/qf.rds",
    database_taxonomy = "output/05_output_files/database_taxonomy.tsv",
    database_metrics = "output/05_output_files/database_metrics.tsv/",
    detected_protein_taxonomy = "output/05_output_files/detected_protein_taxonomy.tsv",
    detected_protein_metrics = "output/05_output_files/detected_protein_metrics.tsv"
  output:
    conduit_obj = "output/05_output_files/conduit_output.rds"
  container: "apptainer/conduitR.sif"
    # From this point- the conduit object contains all the stuff you need. 
    #Load that into Conduit GUI.
  script:
    "scripts/07_ingestion_into_R/02_create_conduit_obj.R"
    
    
