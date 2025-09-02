FROM posit/r-base:4.5-jammy

# ---- all system deps in one cached layer ----
RUN apt-get update && apt-get install -y --no-install-recommends \
    ## basic toolchain --------------------------------------------------------
    build-essential        \
    gfortran               \
    cmake                  \
    pkg-config             \
    ## compression / archives -------------------------------------------------
    zlib1g-dev             \
    libbz2-dev             \
    ## networking / crypto ----------------------------------------------------
    libcurl4-openssl-dev   \
    libssl-dev             \
    ## XML / HTML parsing -----------------------------------------------------
    libxml2-dev            \
    ## graphics & text rendering (ragg, systemfonts, gdtools, svglite, etc.) --
    libcairo2-dev          \
    libpng-dev             \
    libjpeg-dev            \
    libtiff5-dev           \
    libfontconfig1-dev     \
    libfreetype6-dev       \
    libharfbuzz-dev        \
    libfribidi-dev         \
    ## image processing -------------------------------------------------------
    libmagick++-dev        \
    ## geospatial stack (sf, s2, units) ---------------------------------------
    libgdal-dev            \
    libgeos-dev            \
    libproj-dev            \
    libudunits2-dev        \
    ## optimisation / linear algebra (igraph, glmnet, etc.) -------------------
    libblas-dev            \
    liblapack-dev          \
    ## misc dependencies sometimes pulled in ----------------------------------
    libgit2-dev            \
    libssh2-1-dev          \
    libglpk-dev            \
    ## reticulate wants a Python header ---------------------------------------
    python3                \
    python3-dev            \
 && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"

WORKDIR /fc_changes
COPY renv.lock renv.lock

RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

RUN R -e "renv::restore()"

# ---- add Quarto CLI here (small, fast) ----
RUN wget -qO /tmp/quarto.deb \
        https://quarto.org/download/latest/quarto-linux-amd64.deb && \
    apt-get update && \
    apt-get install -y /tmp/quarto.deb && \
    rm /tmp/quarto.deb && \
    rm -rf /var/lib/apt/lists/*

COPY src/ src/
COPY paper/ paper/
COPY data/atlas_data/ data/atlas_data/
COPY data/gradients/ data/gradients/
COPY data/processed_and_cleaned/ data/processed_and_cleaned/


