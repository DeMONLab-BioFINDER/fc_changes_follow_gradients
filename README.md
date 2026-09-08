# Different functional connectivity gradients reflect aging and Alzheimer's disease

This repository contains the analysis code accompanying the manuscript **“Different functional connectivity gradients reflect aging and Alzheimer’s disease.”**

Preprint: [https://doi.org/10.1101/2025.05.22.655469](https://doi.org/10.1101/2025.05.22.655469)

The original participant data cannot be shared. Synthetic BioFINDER and ADNI datasets are included so that the structure of the workflow can be inspected and the analysis code can be exercised without access to restricted data.

The synthetic results are not scientifically meaningful and will not reproduce the numerical results reported in the manuscript. The repository should be understood as a reproducible record of the analysis workflow, with synthetic data provided for testing.

## Reproducibility status

The supported environment is R with `renv`. Running `renv::restore()` will however also install some python dependencies for parts of the workflow that you normally would not need to use but that is kept for provenance. 


The main processing and analyis script is called `main.R`. 

The intended execution pattern is:

1. Run `renv::restore()` to download all packages used in this repo.
2. Run the analysis with `FROM_START=TRUE` first time you run it. This will create connectomes from synthetic fMRI timeseries, calculate derivatives and clean the synthetic datasets just as the real datasets were processed. 
3. The derivatives and cleaned datasets will be saved locally and after that you can run the script with `FROM_START=FALSE` for later runs that reuse the generated connectomes and cached results.

Running everything from start may take one to two hours or more depending on the system. Only running analyses and producing figures will be faster. However, the brain figures do take some time to render. 

Generated connectomes, processed data, figure files, and most other outputs are intentionally excluded from version control. Consequently, `FROM_START=FALSE` is not expected to work immediately after a fresh clone.

The full pipeline has been successfully run from start to finish on a fresh clone. However, it has not yet been tested on a machine other than the one on which the code was developed.

## Requirements

- R 4.6.0 
- Quarto, for rendering `paper/fc_changes_paper.qmd`
- System libraries required by packages such as `sf` and `magick`
- Sufficient memory for arrays containing subject-level 1000 × 1000 connectomes

The analysis has primarily been developed on Ubuntu. Other operating systems may require different system-library installation steps.

Python is not required for the standard synthetic workflow. `requirements.txt` records packages used by optional Python-based processing.

## Installation

Clone the repository and start R from its root directory. The project `.Rprofile` activates `renv` automatically.

Restore the recorded R environment:

```r
renv::restore()
```

If `renv` is not already installed, install it first:

```r
install.packages("renv")
renv::restore()
```

## Running the synthetic workflow

All commands should be run from the repository root.

For the first run:

```bash
FROM_START=TRUE \
CREATE_BRAIN_PERMUTATIONS=FALSE \
EXTRACT_TIMESERIES=FALSE \
REAL_DATA=FALSE \
Rscript src/main.R
```

The equivalent commands from an interactive R session are:

```r
Sys.setenv(
  FROM_START = "TRUE",
  CREATE_BRAIN_PERMUTATIONS = "FALSE",
  EXTRACT_TIMESERIES = "FALSE",
  REAL_DATA = "FALSE"
)
source("src/main.R")
```

After a successful from-scratch run, cached results can be reused with:

```bash
FROM_START=FALSE REAL_DATA=FALSE Rscript src/main.R
```


### Configuration variables

| Variable | Default | Purpose |
| --- | --- | --- |
| `FROM_START` | `FALSE` | When `TRUE`, calculate connectomes, connectivity derivatives, gradients, and expensive window analyses. When `FALSE`, reuse locally generated results. |
| `CREATE_BRAIN_PERMUTATIONS` | `FALSE` | Generate a new set of cortical rotations for spatial permutation tests. The repository already includes the rotations used by the workflow. |
| `EXTRACT_TIMESERIES` | `FALSE` | Extract parcel time series from raw NIfTI images using Python and Nilearn. Raw images are not included, so this should remain `FALSE` for the supplied synthetic data. |
| `REAL_DATA` | `FALSE` | Select restricted BioFINDER and ADNI source-data paths. These data are not distributed; leave this `FALSE` outside the authorized analysis environment. |

## Outputs

A from-scratch run creates or updates the following local outputs:

- `data/bf_src_data_synthetic/connectomes/`: BioFINDER synthetic connectomes.
- `data/adni_src_data_synthetic/connectomes/`: baseline ADNI synthetic connectomes.
- `data/processed_and_cleaned/`: cleaned data, connectivity derivatives, gradients, and window-analysis caches.
- `paper/figures/`: main and extended-data figures, primarily as PDF files.
- `paper/suppfig_original/jpg/`: supplementary figures embedded by the manuscript.
- `paper/figures_source_data/`: CSV source data written alongside the figures.
- `paper/tables/`: tables consumed by the Quarto manuscript.
- `paper/fc_changes_paper.docx`: the rendered manuscript.

These outputs are generally ignored by Git to keep the archived source repository smaller and to avoid mixing generated synthetic results with the published results.

### Manuscript-rendering limitations

The methods figure at `paper/figures/conceptual_plot/methods_plot_grayed.pdf` was assembled manually in Inkscape from components generated during development. The final assembled PDF is retained because it cannot be recreated automatically by the current scripts.

The parcel-wise tau requires data that I have not yet synthesized. Therefore, the synthetic workflow does not generate `paper/figures/nodal_tau.pdf`, although the manuscript currently references that figure. 

After all required figures and tables are available, the manuscript can be rendered separately with:

```bash
Rscript src/render_paper.R
```

## Repository structure

```text
.
├── data/
│   ├── adni_src_data_synthetic/    # synthetic ADNI metadata and time series
│   ├── bf_src_data_synthetic/      # synthetic BioFINDER metadata and time series
│   ├── atlas_data/                 # atlas geometry, labels, rotations, and reference data
│   └── gradients/                  # reference functional-gradient data
├── paper/
│   ├── _extensions/                # local Quarto/Pandoc extensions
│   ├── fc_changes_paper.qmd        # manuscript source
│   ├── references.bib              # bibliography
│   └── custom_ref_new.docx         # Word reference template
├── renv/                           # renv activation and settings
├── src/                            # analysis and visualization code
├── renv.lock                       # pinned R package environment
└── requirements.txt                # optional Python dependencies
```

## Source-code guide

### Main workflow

#### `src/main.R`

The entry point for the analysis. It performs the following broad stages in sequence:

1. Reads configuration variables and prepares output directories.
2. Loads the Schaefer-1000 atlas, Yeo network labels, cortical geometry, and spin-test permutations (if you don't want to calculate them yourself).
3. Reads and cleans the synthetic BioFINDER metadata.
4. Constructs BioFINDER connectomes from parcel time series.
5. Calculates nodal connectivity strength and a variety of connectivity-similarity/affinity measures.
6. Derives functional gradients and evaluates alternative gradient parameters.
7. Reads baseline ADNI metadata, matches scan-level time series to subject-level motion files, performs motion scrubbing, and calculates replication-cohort connectomes and derivatives.
8. Runs cross-sectional, nonlinear, longitudinal, cognition, mediation, and sensitivity analyses. See the paper for more information. 
9. Writes figures, figure source data (only statistics), and descriptive tables.
10. Invokes `src/render_paper.R`.

`main.R` is intentionally a sequential analysis script rather than an R package or workflow-manager project. Run it from the repository root because paths are relative to that location.

## Scripts

`src/util.R`

Contains the main numerical and modelling helpers, including:

- nodal strength, within-network, and between-network connectivity; connectivity similarity/affinity;
- vectorised parcel-wise linear and mixed-effects regression;
- extraction of nodal model estimates (getting t-values from the vectorised models);
- parcelwise GAM related functions;
- main function for the sensitivity table of the supplementary.

Several functions rely on atlas objects created near the beginning of `main.R`, so this file is not designed as a standalone library.

`src/util_gradients.R`

Implements gradient construction and alignment. It contains the diffusion-map implementation, PCA/diffusion gradient estimation, component reordering, sign alignment to reference gradients, and optional visualization.

### Visualization helpers

`src/plot_grad_rels.R`

Fits or receives parcel-wise models and builds the cortical maps and scatterplots used to show relationships between model-effect maps and functional gradients.

`src/plot_gams.R`

Builds the nonlinear-analysis figure from parcel-wise generalized additive model predictions and derivatives. The file contains the current plotting implementation and retained legacy wrappers.

`src/util_vis.R`

Composes higher-level manuscript figures from the lower-level plotting functions. This includes the main cross-sectional figures, longitudinal/window figures, network overlays, gradient comparisons, and shared layout helpers.

`src/plot_gradient_coverage.R`

Creates the figure describing the anatomical, network, and cognitive-term coverage of the principal gradients. It uses the precomputed NeuroQuery result stored in `data/atlas_data/schaefer1000_NQ_results.rds`, which you can get using code in `neurosynth.py`. 

`src/mass_mediation_src.R`

Contains the parcel-wise mediation function. This function makes it possible to run 1000 mediation analyses in a matter of seconds. 


## Citation

If you use this code, please cite:

Rittmo, J., Franzmeier, N., Strandberg, O., Chauveau, L., Satterthwaite, T. D., Wisse, L. E., Spotorno, N., Behjat, H. H., Dehsarvi, A., van Westen, D., Anijärv, T. E., Palmqvist, S., Janelidze, S., Stomrud, E., Ossenkoppele, R., Mattsson-Carlgren, N., Hansson, O., & Vogel, J. W. (2025). *Different functional connectivity gradients reflect aging and Alzheimer’s disease*. bioRxiv. [https://doi.org/10.1101/2025.05.22.655469](https://doi.org/10.1101/2025.05.22.655469)

## License

The code is released under the terms in `LICENSE`.
