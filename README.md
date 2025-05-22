# Conduit

Conduit is a Snakemake workflow that provides a path through a mountain of
metaproteomics data produced by the Astral mass spectrometer. 
It facilitates a comprehensive analysis of DIA metaproteomics 
data using multiple strategies to define the search space
(generate a database). 

## Search Space Definition Methods

The workflow supports the following methods to define the search space:

1. **User-defined Organisms**: Using NCBI taxonomy IDs or UniProt proteome IDs
2. **Proteotyping**: Uses a first pass search with species-specific tryptic peptides to define the search space.
3. **Metagenomics**:
   - Reference-based via ?
   - Assembly-based (MAGs) via ?
4. **16S rRNA Analysis**

The workflow produces files that are intended to be used with the 
Conduit-GUI, enabling users to visualize and draw biological 
conclusions from the data.

Conduit-GUI can be found at the following links:
Github: https://github.com/baynec2/conduit-GUI
Deployment: https://gonzalezlab.shinyapps.io/conduit-GUI/

## Table of Contents
- [Installation](#installation)
- [Quick Start Guide](#quick-start-guide)
- [Input Requirements](#input-requirements)
- [Workflow Overview](#workflow-overview)
- [Configuration](#configuration)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)
- [Dependencies](#dependencies)

## Installation

### Prerequisites
- Linux-based operating system (Ubuntu 20.04 or later recommended)
  - Note: While other operating systems can be used, Apptainer containers won't work, requiring manual dependency management (not recommended)
- Python 3.8 or later
- Minimum 100GB free disk space
- Internet connection for downloading dependencies
- Recommended: High RAM capacity for handling large search spaces

### Required Software Installation

1. **Snakemake**
   ```bash
   pip install snakemake
   ```

2. **Apptainer**
   ```bash
   sudo apt-get update
   sudo apt-get install -y apptainer
   ```

3. **Conduit**
   ```bash
   git clone "https://github.com/baynec2/conduit"
   cd conduit
   ```

## Quick Start Guide

### 1. Prepare Configuration

Create a new experiment directory and prepare the following configuration:

```bash
mkdir -p experiments/your_experiment/config
```

Create `experiments/your_experiment/config/snakemake_config.yaml` with the following options:

```yaml
raw_files: "path/to/raw_files"
method: "your_chosen_method"  # See available methods below
fastq_files: "path/to/fastq/files"
output: "path/to/output/direcotry"
```

Available methods:
- `ncbi_taxonomy_id`: Use NCBI taxonomy IDs
- `uniprot_proteome_id`: Use UniProt proteome IDs
- `proteotyping`: Use species-specific tryptic peptides
- `MAG`: Use Metagenome-Assembled Genomes
- `metagenomic_profiling`: Use metagenomic profiling
- `16S`: Use 16S rRNA gene sequences

### 2. Prepare Input Files

Organize your input files in the following structure:

```
experiments/your_experiment/
├── config/
│   └── snakemake_config.yaml
├── input/
│   ├── raw/           # Place .raw files here
│   └── fastq/         # Place fastq files here (if needed)
├── sample_annotation.txt
└── organisms.txt      # If using taxonomy/proteome IDs
```

### 3. Run the Workflow

```bash
snakemake --use-singularity --cores 4 --configfile experiments/your_experiment/config/snakemake_config.yaml
```

### 4. Monitor Progress

- Check progress in `experiments/your_experiment/log/`
- Find output files in `experiments/your_experiment/output/`

## Input Requirements

### Essential Files

1. **Configuration File**
   - Location: `experiments/your_experiment/config/snakemake_config.yaml`
   - Required parameters: `samples` and `method`

2. **Raw Files**
   - Location: `experiments/your_experiment/input/raw/`
   - Format: Thermo Scientific .raw files
   - Naming: Avoid special characters in filenames

3. **Sample Metadata**
   - File: `experiments/your_experiment/input/sample_annotation.txt`
   - Format: Tab-delimited text file
   - Required columns:
     - `file`: Raw file name (without .raw extension)
   - Optional columns: treatment, timepoint, replicate, etc.
   - Example:
     ```
     file    treatment    timepoint    replicate
     sample1 control      0h           A
     sample2 control      0h           B
     sample3 treatment    24h          A
     sample4 treatment    24h          B
     ```

### Method-Specific Requirements

#### For NCBI Taxonomy ID Method
- File: `experiments/your_experiment/input/organisms.txt`
- Format: Tab-delimited text file
- Required columns: `organism_id`, `organism_type`
- Example:
  ```
  organism_id    organism_type
  9606           host
  83333          microbiome
  818            microbiome
  3847           diet
  ```

#### For UniProt Proteome ID Method
- File: `experiments/your_experiment/input/organisms.txt`
- Format: Tab-delimited text file
- Required columns: `proteome_id`, `organism_type`
- Example:
  ```
  proteome_id    organism_type
  UP000005640    host
  UP000000625    microbiome
  UP000095541    microbiome
  UP000008827    diet
  ```

#### For MAG, Metagenomic Profiling, and 16S Methods
- Place fastq files in `experiments/your_experiment/input/fastq/`
- For 16S method: Place fasta files in `experiments/your_experiment/input/fasta/`

### Optional Configuration

You can customize DIA-NN settings by creating:
- `experiments/your_experiment/config/generate_diann_spectral_library.cfg`
- `experiments/your_experiment/config/run_diann.cfg`

If not specified, default configurations will be copied from `config/` to ensure reproducibility.

## Workflow Overview

The workflow consists of several main stages:

1. **Database Resource Generation**
   - Downloads and processes protein databases
   - Generates taxonomy information
   - Creates spectral libraries

2. **Data Processing**
   - Processes raw files
   - Performs DIA analysis
   - Generates protein identifications

3. **Output Generation**
   - Creates various output matrices
   - Generates visualization files
   - Produces summary statistics

## Configuration

### Main Configuration Files
- `config/generate_diann_spectral_library.cfg`: Spectral library generation settings
- `config/run_diann.cfg`: DIA-NN analysis parameters

### Customizing the Workflow
1. Edit configuration files in the `config/` directory
2. Modify Snakemake parameters in the command line
3. Adjust resource allocation in the Snakefile

## Output Files

### Main Output Directories
- `output/00_database_resources/`: Database and taxonomy files
- `output/02_speclib_gen/`: Spectral library files
- `output/03_diann_output/`: DIA-NN analysis results
- `output/05_output_files/`: Final output matrices and visualizations

### Key Output Files
- `conduit_output.rds`: Main output object for Conduit-GUI
- `*_matrix.tsv`: Various quantification matrices
- `qf.rds`: QFeatures object for R analysis

## Troubleshooting

### Common Issues
1. **Container Build Failures**
   - Ensure you have the correct DIA-NN license
   - Check apptainer installation
   - Verify sufficient disk space

2. **Memory Issues**
   - Increase available RAM
   - Reduce number of parallel processes
   - Check for memory leaks in logs

3. **Input File Problems**
   - Verify file formats
   - Check file permissions
   - Ensure correct file naming

### Getting Help
- Check the [GitHub Issues](https://github.com/baynec2/conduit/issues)
- Contact the development team
- Review the logs in the `log/` directory

## Contributing

We welcome contributions! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request
4. Follow our coding standards

## License

MIT License

Copyright (c) 2025 Charlie Bayne

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Dependencies

Conduit relies on the following key software:

1. **DIA-NN 2.1.0**
   - Performs DIA analysis
   - Handles spectrum matching and protein identification
   - Requires a valid license

3. **conduitR**
   - R package for search space definition and data analysis
   - Dependencies include:
     - QFeatures
     - SummarizedExperiment
     - ggplot2
     - metacoder
     - dplyr
     - tidyr
     - And other R packages
