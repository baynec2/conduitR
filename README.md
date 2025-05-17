# Conduit

Conduit is a snakemake workflow that provides a path through a mountain of
metaproteomics data produced by the Astral mass spectrometer. 
It is intended to facillitate a comprehensive analysis of DIA metaproteomics 
data using a number of strategies to define the search space
(generate a database). 

These are as follows:

1. User defined by providing organism_ids.  
2. Proteotyping
3. Metagenomics
  * Reference based via Kraken2
  * Assembly based


It produces files that are intended to be used in combination with the 
Conduit-GUI allowing the user the ability to visualize and draw biological 
concliusions from the data.

Github: https://github.com/baynec2/conduit
Deployement: https://gonzalezlab.shinyapps.io/conduit/


# Instructions

Will add later. Basic idea is below.

1. Need to have snakemake installed.
2. Need to install apptainer on system.
3. Need to build docker image for DIANN, as license does not allow for the distribution of prebuilt image

# Dependancies

Conduit relies heavily on preexisting software to function. These are as follows:

1. Thermorawfileparser - used to convert .raw to mzml to be compatible with linux
2. DIA-NN - used to perform the computational intense proteomics analysis. Does the 
important stuff. Figures out what specra beclong to 
3. conduitR - Used to define search space, pull data from APIs, make plots, etc. 
  * This itself depends on numerous R packages 
  * QFeatures,SummarizedExperiment, ggplot2, metacoder, dplyr, tidyr,etc.
