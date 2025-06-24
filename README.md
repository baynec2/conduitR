# Conduit: A Modular Metaproteomics Analysis Platform

Conduit is a scalable and modular workflow management system for metaproteomics data analysis, designed to seamlessly integrate with metagenomic data if avalible. Built using Snakemake, it provides a robust pipeline for processing Data Independent Acquisition (DIA) mass spectrometry data with particular emphasis on metaproteomics applications.


## Features

Note this is currently a work in progress and all features are not yet avalible. 

- **Search Space Definition**: Choose the best way to define the search space for *your* experiment.
    - **User Defined Taxa**: Know what taxa are in your sample? Great! Just give Conduit the ncbi taxa ids, we will do the rest.
    - **User Definied Proteomes**: Know the specific proteome ids in your sample? Also amazing, we can work with those too. 
    - **Reference Based Metagenomic**: Don't know what taxa are in your sample, but were able to get metagenomic sequencing data? No problem! We can use it to define the search space.
    - **MAGS**: Have metagenomic sequencing data, but reference based metagenomic data not complicated enough for you? No problem, we will make MAGs and use those to defne the search space.
    - **16S**: Don't know what is in your sample, can't get your hands on any metagenomic sequencing data, but have 16S sequencing data? Definitely the worst option, but better than nothing. We will do our best to use it.
- **Metadata Handling**: Easily associate your sample data with the infromation about what each sample actually is. Less time spent on data wrangling means more time to figure out what it all means.
- **DIA-NN Integration**: Automated, spectral library free processing of DIA data.
- **Taxonomic Analysis**: Multi-level taxonomic classification of proteins
- **Functional Analysis**: GO term, KEGG pathway, and subcellular location annotation
- **R Integration**: Direct integration with R to enable easy statistical analysis, plotting, and more. 
- **Inegration with a Dedicated GUI**: Want to explore your data quickly without having to write any code? Conduit-GUI has you covered. 
- **Container Support**: Full containerization via Apptainer. Run it on any reasonable linux machine with ease.
- **Scalable**: Conduit runs on a single machine or scales to HPC clusters. Tip: metaproteomics loves compute — the more, the better.
- **Open Source**: Conduit is open source, allowing users to customize the pipeline to their specific needs.

## Dependencies

### Required Software
- Snakemake (≥7.0.0)
- conduitR
- Python libraries
   * glob
   * pandas 
   * shutil
   * pdb  
   * logging
   * datetime
- DIA-NN (2.1.0)
- Apptainer/Singularity (≥1.1.0)

Note: All other dependencies (R, Python, DIA-NN, etc.) are automatically handled through Apptainer containers. You only need to ensure that Snakemake and Apptainer are installed on your *linux* system.

If not using linux or apptainer, you will need to install the dependencies manually. Will add more details here soon.

### Hardware Requirements
- 64GB+ RAM recommended.
- Multi-core processor (8+ cores recommended)
- Storage space for raw data and results. Raw files are usually around ~7 GB each, so it gets to Tb scale very quickly.

## Quick Start

1. Clone the repository:
```bash
git clone https://github.com/baynec2/conduit.git
cd conduit
```

2. Set up an experiment directory following the example structure:
```
experiments/your_experiment/
├── input/
│   ├── raw_files/        # Your .raw files
│   ├── organisms.txt     # Organism IDs
│   └── sample_annotation.txt
├── config/
│   └── snakemake.yaml    # Experiment configuration
└── output/               # Results will be here
```

3. Run the workflow:
```bash
snakemake --configfile experiments/your_experiment/config/snakemake.yaml --use-apptainer
```
4. Run [Conduit-GUI](https://github.com/baynec2/conduit-GUI) to explore your data. Or use [conduitR](https://github.com/baynec2/conduitR) if you prefer a more bespoke experience. Up to you.

5. Repeat for each experiment you want to run (keeping a reproducible record of your experiments).

## Project Structure

```
conduit/
├── modules/                          # Snakemake modules. Each module should have snakefile and associated scripts. 
│   ├── search_space/                 # Defining the search space
        ├── user_specified            # with ncbi_taxonomy and proteome_id 
        ├── proteotyping              # with species specific peptides
        ├── metagenomic_profiling     # with metagenomic profiling (reference based_)
        ├── mags                      # with metagenome assembled genome
        └── 16s                       # with 16S data
│   ├── diann/                        # Identification and Quantification with DIA-NN
│   ├── annotation/                   # Protein / taxonomic annotation
│   |── matrices/                     # Matrix processing
|   └── r_integration                 # Integration into R
├── config/                           # Default configurations (copied to each experiment config if not altered)
├── experiments/                      # Experiment directories
│   └── example/                      # Example experiment
|       |── config                    # Configuration files for specific experiment
|       |── logs                      # Log files for each rule. 
        └─── input                     # Inputs into the workflow
            |── sample_annotation.txt # Sample annotation file matching names of .raw MS files. 
            |── organisms.txt         # Text file NCBI organism IDs defining the search space (if using ncbi_taxonomy)
            |── proteome_id.txt       # Text file containing proteome ids (if using proteome id workflow)
            |── fastq                 # Fastq input for metagenomic profiling or 16S
            |── mags                  # MAG inputs.
            └── raw_files             # Thermo Mass Spectrometry .raw files. 

└── tests/                            # Test suite
```
## Configuration

The main configuration is done through `snakemake.yaml`. Key parameters include:

```yaml
experiment_directory: "experiments/example"
method: "ncbi_taxonomy_id"
raw_files: "experiments/example/input/raw_files/*.raw"
sample_annotation: "input/sample_annotation.txt"
```
See the example configuration in `experiments/example/config/` for a complete template.

## Output Structure

The workflow generates several key outputs:

1. **Database Resources**
   - Protein FASTA files
   - Taxonomy information

2. **DIA-NN Results**
   - Spectral library
   - Protein/peptide quantification

3. **Detected Protein Annotations**
   - GO term annotations
   - KEGG pathway annotations
   - Subcellular location annotations

4. **Processed Matrices**
   - Taxonomic level matrices
   - Functional annotation matrices
   - Protein-level matrices
   - Peptide-level matrices
   - Precursor-level matrices

5. **R Objects**
   - QFeatures object
   - Analysis metrics
   - Conduit object (accepted by Conduit-GUI/conduitR)

## Acknowledgments

Conduit would not be possible without the great work of many people making great software that it relies on. These are as follows:

- Snakemake developers
- DIA-NN developers
- R Bioconductor community
- R tidyverse community
- R 
- [Other acknowledgments]  

