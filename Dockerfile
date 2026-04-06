FROM rocker/bioconductor:3.22

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

# Install pak for fast parallel package installation
RUN Rscript -e "install.packages('pak', repos = 'https://r-lib.github.io/p/pak/stable/')"

# Install all R dependencies (CRAN, Bioconductor, GitHub) in one pass
RUN Rscript -e "pak::pak(c( \
    'bioc::QFeatures', \
    'bioc::Biostrings', \
    'bioc::SummarizedExperiment', \
    'bioc::S4Vectors', \
    'bioc::limma', \
    'bioc::MsCoreUtils', \
    'bioc::ComplexHeatmap', \
    'bioc::PCAtools', \
    'bioc::sechm', \
    'bioc::clusterProfiler', \
    'grunwaldlab/metacoder', \
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
))"

# Install conduitR from the local source
COPY . /pkg
RUN Rscript -e "pak::pak('local::/pkg')"

CMD ["R"]
