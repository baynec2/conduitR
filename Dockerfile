FROM rocker/bioconductor:latest

# System dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libgit2-dev \
    libhdf5-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libtiff-dev \
    libgdal-dev \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Bioconductor packages
RUN Rscript -e "BiocManager::install(c( \
    'QFeatures', \
    'Biostrings', \
    'SummarizedExperiment', \
    'S4Vectors', \
    'limma', \
    'MsCoreUtils', \
    'ComplexHeatmap', \
    'PCAtools', \
    'sechm', \
    'clusterProfiler' \
), ask = FALSE, update = FALSE)"

# GitHub packages
RUN Rscript -e "remotes::install_github('grunwaldlab/metacoder')"

# CRAN packages
RUN Rscript -e "install.packages(c( \
    'future', \
    'furrr', \
    'purrr', \
    'tibble', \
    'httr2', \
    'KEGGREST', \
    'rentrez', \
    'XML', \
    'tidyr', \
    'dplyr', \
    'readr', \
    'stringr', \
    'ggplot2', \
    'ggprism', \
    'ggraph', \
    'ggkegg', \
    'tidygraph', \
    'cowplot', \
    'viridis', \
    'heatmaply', \
    'sunburstR', \
    'd3r', \
    'arrow', \
    'jsonlite', \
    'glue', \
    'rlang', \
    'assertthat', \
    'httr', \
    'Matrix', \
    'xgboost', \
    'parsnip', \
    'recipes', \
    'rsample', \
    'workflows', \
    'tune', \
    'dials', \
    'hardhat', \
    'yardstick', \
    'testthat', \
    'withr', \
    'imputeLCMD' \
), repos = 'https://cloud.r-project.org')"

# Install conduitR from the local source
COPY . /pkg
RUN Rscript -e "remotes::install_local('/pkg', dependencies = FALSE)"

CMD ["R"]
