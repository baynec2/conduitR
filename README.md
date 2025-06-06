
<!-- README.md is generated from README.Rmd. Please edit that file -->

# conduitR <img src="man/figures/logo.png" align="right" height="139" />

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

## Overview

conduitR is an R package that provides tools for metaproteomics analysis
using conduit and conduit-GUI. It offers a comprehensive suite of
functions for processing, analyzing, and visualizing metaproteomics
data. These functions are designed to be used in conjunction with the
[conduit](https://github.com/baynec2/conduit) and
[conduit-GUI](https://github.com/baynec2/conduit-GUI) packages, or
independently for custom analyses.

## Features

- **Getting Data**
  - Download FASTA files.
  - Download taxonomy data from NCBI
- **Data Processing**
  - Create custom protein sequence databases
  - Associate sample proteomics data and metadata using QFeatures
    integration.
  - Transformation, Imputation, and Normalization, of proteomics data.
- **Statistical Analysis**
  - Perform LIMMA differential expression analysis
  - Support for various machine learning models (LASSO, Random Forest,
    XGBoost)
  - Feature importance analysis
- **Visualization**
  - Interactive heatmaps with annotations

  - PCA biplots with customizable aesthetics

  - Taxonomic barplots

  - Taxonomic heat trees

  - Volcano plots for differential expression

  - KEGG pathway visualization

  - Feature-specific plots

  - Custom plot aesthetics

  - **Miscellaneous**

  - Logging functions.

## Installation

You can install the development version of conduitR from GitHub:

`r, eval = FALSE # install.packages("devtools") devtools::install_github("baynec2/conduitR")`

## Usage

This is a work in progress, will update in the future.
