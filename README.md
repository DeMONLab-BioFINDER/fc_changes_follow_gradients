## Functional Connectivity Changes in Aging and Alzheimer's Disease

This repository provides code for reproducing analyses and figures from the manuscript:

**"Age and Alzheimer’s disease affect functional connectivity along separate axes of functional brain organization"**

[Preprint](https://doi.org/10.1101/2025.05.22.655469)

Due to data-sharing restrictions, the repository uses synthetic data by default.


ℹ️ **NOTE**

The source code in the repo is for the final version of the manuscript, the code in the docker container produces the version of the manuscript
that is published as a preprint. You won't be able to run all of the source code in the current repo with the synthetic data. As soon as I have time
I will synthesize new data and update everything so that you can run everything from start to finish for the final version.

---

## Getting Started

You can either clone the repo as it is and try to set it up using in R using `renv::restore()` and then
source `src/main.R`, but the easiest way to reproduce the analyses and manuscript is likely to use a docker
container with the dependencies.

### Docker Image

First install docker https://www.docker.com/get-started/

The Docker image is available on Docker Hub:

```
docker pull jorittmo/fc-changes-image:latest
```

### Running the Analysis with Docker

To reproduce the paper using the synthetic data provided clone the repo to 
a directory of your choice and then run:

```
docker run --rm \
  -v /fc_changes_follow_gradients/data/adni_src_data_synthetic:/fc_changes/data/adni_src_data_synthetic \
  -v /fc_changes_follow_gradients/data/bf_src_data_synthetic:/fc_changes/data/bf_src_data_synthetic \
  -v /fc_changes_follow_gradients/docker_out:/fc_changes/paper/output \
  -e FROM_START=FALSE \
  jorittmo/fc-changes-image \
  Rscript src/main.R
```

This will run all analyses (without pre-calculating measures etc) and spit out the manuscript with figures usin synthetic data in `docker_out`.

Ensure the host mount paths are the absolute path to the repo (i.e. `~/fc_changes_follow_gradients`, if you put it in you home folder).

* `-v` mounts local directories into the container for data input and output.
* `-e FROM_START` sets if you want to do all pre-computation to get derivatives, gradients and such. Set to `TRUE` to recompute everything from scratch and to `FALSE` to only run analyses and models and generate figures.

If you want to try the image without running the full analysis pipeline which takes ~10-30 minutes depending on setup, and if you run it from start
you can just render the paper by doing:

```
docker run --rm \
  -v /fc_changes_follow_gradients/docker_out:/fc_changes/paper/output \
  -e FROM_START=FALSE \
  jorittmo/fc-changes-image \
  Rscript src/render_paper.R
```

The rendered paper (with synthetic data) should then be reachable in the folder `docker_out`.

---

## Configuration and Arguments

The script `main.R` recognizes the following environmental arguments:

| Argument                    | Description                                                                                  | Recommended Setting (Synthetic Data) |
| --------------------------- | -------------------------------------------------------------------------------------------- | ------------------------------------ |
| `FROM_START`                | Recompute everything (connectome calculation, FC derivatives, gradients) from scratch        | `FALSE` (data precomputed)           |
| `CREATE_BRAIN_PERMUTATIONS` | Generate brain permutations for spin tests  (takes quite long to run and is not necessary)   | `FALSE` (already included)           |
| `EXTRACT_TIMESERIES`        | Extract parcel time series from raw images  (should not be used for now)                     | `FALSE` (timeseries pre-extracted)   |
| `REAL_DATA`                 | Toggle between real and synthetic data  (Should not be used for now)                         | `FALSE` (synthetic data)             |

To override the defaults, set them as environmental variables when calling Docker (e.g., `-e FROM_START=TRUE`).

---

## Workflow Summary

The code in `main.R` runs the following major steps:

1. **Initial setup**

   * Set flags (`FROM_START`, `CREATE_BRAIN_PERMUTATIONS`, `EXTRACT_TIMESERIES`, `REAL_DATA`), create folders.

2. **Atlas Loading**

   * Load/parcellation (Schaefer, Yeo), build ROI-to-network lookup, gradient color palettes.

3. **(Optional) Permutations**

   * If requested, generate 1,000 spin/test parcellation rotations and save them.

4. **(Optional) Timeseries Extraction**

   * If requested, use Nilearn to read raw NIfTIs, extract 1,000-parcel time series, save as `.rds`.

5. **Load & Clean BioFINDER Data**

   * Read CSV → filter for motion & diagnosis/pathology → compute continuous pathology trajectories via SCORPIUS → finalize `biofinder_df`.

6. **(FROM\_START) Build BioFINDER Connectomes & Derivatives**

   * Read saved timeseries → compute 1000×1000 correlation matrices → save → read them back into a big 3D array → compute node-strength + two "affinity” measures via `fc_strength`/`get_affinity` → save those results.

7. **(FROM\_START) Compute BioFINDER Gradients & Parameter Sweep**

   * Load “Margulies” reference gradients (CSV or via Nilearn), build “healthy-young” average connectome → plot supplementary gradient-comparison figure and save → call `get_gradients` once at (threshold=0, method=pca; the gradients used) → then loop over several{threshold, method, affinity} combos to build a large dataframes with gradients derived over different values.

8. **Load & Clean ADNI Data**

   * Read CSV → rename columns → compute continuous `pathology_ad` via SCORPIUS → load FD CSVs (not computed in this workflow) → (if `FROM_START`) extract/scrub timeseries → build ADNI connectomes → filter by motion → save cleaned `adni_df`.

9. **(FROM\_START) Compute ADNI Connectivity Derivatives**

   * Read 3D array of ADNI connectomes → compute node-strength + two affinity variants → save.

10. **(FROM\_START) Compute ADNI Gradients & Merge**

    * Build “young, healthy” average ADNI connectome → call `get_gradients` once (PCA) → bind with BioFINDER gradients → loop over parameter grid to build `gradient_data_adni` → merge with BioFINDER → save a

11. **(Optional) Methods Figure** 
    
    * Run code to generate pngs used in the methods plot, this is now commented out in main.R, just uncomment our source `methods_figure.R` to generate the figures used in it. This plot is not reproducible as it was put together in inkscape. All other figures are generated.

12. **Main Figures**

    * **Figure 1**: `figure_one(...)` cross-sectional biofinder + ADNI (cognition is also generated but we take only the first figure panels), save as PNG
    * **Figure 2**: first calculate non-linear trajectories: `gam_pred_nodes(...)` → `plot_gams_v1(...)`, save as PNG
    * **Figure 3**: Build "longitudinal" df → `longitudinal_and_window_analysis(...)`, save main & supp PNG
    * **Figure 4**: `figure_one(...)` cognition panels only, save as PNG

13. **Supplementary Figures**

    * Many calls to `plot_gradient_relationships(...)`, `figure_one(...)`, `overlaid_main_results(...)` with different and/or measures (diagnosis, gradient 2 only, within/between-network affinity, nodal strength, no-threshold affinity, clinical interactions, healthy group interactions, pathology densities etc.). Each saved to a distinct file under `paper/figures/…`.

14. **Table 1 (Cross-Sectional)**

    * Build descriptive `finalfit` tables for BioFINDER & ADNI, row-bind, save as `CS_tbl1.rds`.

15. **Longitudinal Table 1 (BioFINDER)**

    * Compute baseline/follow-up summary for each BioFINDER subject → `summary_factorlist(...)`, save as `LT_tbl1.rds`.

16. **(FROM\_START) Gradient Sensitivity Analyses (Supplementary Table)**

    * If `FROM_START`, call `write_gradient_supp_table(params)`, save as `supp_table.rds`; otherwise just load.

17. **Quarto Rendering**

    * Source the script render_paper.R which calls quarto to render the manuscript and then moves the rendered paper into paper/output

    
---

## Troubleshooting

If you encounter any issues, please open an issue on GitHub, and we'll assist you as quickly as possible.

The non-container version of this code has only been tested locally on ubuntu 22.04. If you have dependency troubles I recommend
the containerized version. 

The code is not as well documented as I would have wanted at this time, so please raise an issue or contact me directly for further inquiries.

---

## Citation

If you use this code, please cite:

Rittmo, J., Franzmeier, N., Strandberg, O., Chauveau, L., Satterthwaite, T. D., Wisse, L. E., Spotorno, N., Behjat, H. H., Dehsarvi, A., Westen, D. van, Anijärv, T. E., Palmqvist, S., Janelidze, S., Stomrud, E., Ossenkoppele, R., Mattsson-Carlgren, N., Hansson, O., & Vogel, J. W. (2025). Age and Alzheimer’s disease affect functional connectivity along separate axes of functional brain organization (p. 2025.05.22.655469). bioRxiv. https://doi.org/10.1101/2025.05.22.655469

---

