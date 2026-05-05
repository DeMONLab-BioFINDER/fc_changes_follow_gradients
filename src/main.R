 ######  ######## ######## ##     ## ########  
##    ## ##          ##    ##     ## ##     ## 
##       ##          ##    ##     ## ##     ## 
 ######  ######      ##    ##     ## ########  
      ## ##          ##    ##     ## ##        
##    ## ##          ##    ##     ## ##        
 ######  ########    ##     #######  ##           
 

##########
# Packages
##########

# To install all dependencies, for _the whole_ project uncomment and run the snippet below
# renv::restore()
# However, you only need the python dependencies if you want to extract timeseries from images

library(tidyverse)
library(broom)
library(SCORPIUS)
library(conflicted)
library(magick)
library(sf)
library(readxl)
library(ggpmisc)
library(cowplot)
library(pbapply)
conflicts_prefer(ggpp::annotate)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::lag)


# These arguments are set as environmental variables when running the docker image
# if you are running the code on your own machine, set them manually

from_start <- as.logical(Sys.getenv("FROM_START", "FALSE"))
create_brain_permutations <- as.logical(Sys.getenv("CREATE_BRAIN_PERMUTATIONS", "FALSE"))
extract_timeseries <- as.logical(Sys.getenv("EXTRACT_TIMESERIES", "FALSE"))
real_data <- as.logical(Sys.getenv("REAL_DATA", "TRUE"))

figure_path <- "paper/figures"
dir.create(figure_path, showWarnings = FALSE)
source_figure_path <- "paper/figures_source_data"
dir.create(source_figure_path, showWarnings = FALSE)
table_path <- "paper/tables"
dir.create(table_path, showWarnings = FALSE)

data_dir <- ifelse(real_data, "data/bf_src_data",  "data/bf_src_data_synthetic")
dir.create(data_dir, showWarnings = FALSE)

ts_dir <- file.path(data_dir, "timeseries")
connectome_dir <- file.path(data_dir, "connectomes")
df_file <- file.path(data_dir, "biofinder_df.csv")

data_dir_adni <- ifelse(real_data, "data/adni_src_data",  "data/adni_src_data_synthetic")
dir.create(data_dir_adni, showWarnings = FALSE)

ts_dir_adni <- file.path(data_dir_adni, "timeseries")
connectome_dir_adni <- file.path(data_dir_adni, "connectomes")
df_file_adni <- file.path(data_dir_adni, "adni_df.csv")

clean_dir <- "data/processed_and_cleaned"


##############
# Atlas setup
##############

atlas_dir <- "data/atlas_data"

schaef1k <- readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds"))
schaef1k <- readRDS(file.path(atlas_dir, "schaef_ggseg2.rds"))
yeo7 <- readRDS(file.path(atlas_dir, "Yeo2011_7_geometry.rds"))

# Yeo 7 network colors and names
net_names <- data.frame(name = c('Vis', 'SomMot', 'DorsAttn','SalVentAttn','Limbic', 'Cont', 'Default'),
                        col = c("#781286", "#4682B4", "#00760E", "#C43AFA", "#c7cc7a", "#E69422", "#CD3E4E"), 
                        label = c(1:7))

# How they map to the Schaefer 1000 parcels
yeo_msk <- as.numeric(read_lines(file.path(atlas_dir, "org_mask_yeo_1000.txt")))
# Schaefer 1000 parcel names
rois <- readLines(file.path(atlas_dir, "Schaefer2018_1000Parcels_order.txt"))

# Everything put together
roi_data <- data.frame(region = rois, 
                       yeo_label = yeo_msk) |> inner_join(net_names, join_by(yeo_label == label))

# Colors for the different ends of the gradients
gradient_cols <- data.frame(gradient1 = c("#3F596D", "#D38A4E"), 
                            gradient2 = c("#4682B4", "#781286"),  
                            gradient3 = c("#8A6081", "#738518"))


if (create_brain_permutations) {
  # Create the brain permutations to use for spintest, long compute time
  
  # File containing an array of functions to facillitate calculations
  source("src/util.R")
  
  # Coords downloaded from https://github.com/ThomasYeoLab/CBIG/tree/master/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal/Parcellations/MNI/Centroid_coordinates
  coords <- read_csv(file.path(atlas_dir, "Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv"))
  
  perms <- rotate.parcellation(coord.l = coords[1:500, 3:5] |> as.matrix(), 
                               coord.r = coords[501:1000, 3:5] |> as.matrix(), 
                               nrot=1000, 
                               method='hungarian')
  write_rds(perms, file.path(atlas_dir, "permutations_1000_hungarian.rds"))
} 


##        #######     ###    ########     ########     ###    ########    ###    
##       ##     ##   ## ##   ##     ##    ##     ##   ## ##      ##      ## ##   
##       ##     ##  ##   ##  ##     ##    ##     ##  ##   ##     ##     ##   ##  
##       ##     ## ##     ## ##     ##    ##     ## ##     ##    ##    ##     ## 
##       ##     ## ######### ##     ##    ##     ## #########    ##    ######### 
##       ##     ## ##     ## ##     ##    ##     ## ##     ##    ##    ##     ## 
########  #######  ##     ## ########     ########  ##     ##    ##    ##     ##

###########################################
# Images
# -------
# This should only be run if you have
# nifti images that you are extracting
# the timeseries from. In this repo
# we start the analysis with the timeseries
# already extracted.
###########################################


if (extract_timeseries){
  
  require(reticulate)
  nilearn <- import("nilearn", convert = FALSE)
  
  schaefer <- nilearn$datasets$fetch_atlas_schaefer_2018(n_rois = 1000L, data_dir = "data/atlas_data/")
  atlas_file <- schaefer$maps
  
  masker <- nilearn$maskers$NiftiLabelsMasker(
    labels_img = atlas_file,
    smoothing_fwhm = as.integer(6),
    standardize = FALSE,
    verbose = 0L
  )
  
  dir.create(ts_dir, showWarnings = FALSE)
  
  # images should be in subject/session named folders with only the nifti image within
  img_dir <- file.path(data_dir, "images")
  nifti_images <- list.files(img_dir, recursive = TRUE)
  
  i <- 0
  pb <- txtProgressBar(min = i, max = length(nifti_images), style = 3)
  for (img in nifti_images){
    ts <- masker$fit_transform(file.path(img_dir, img))
    ts_r <- py_to_r(ts)
    image_file <- str_split_1(img, "/")[1]
    write_rds(ts_r, file.path(ts_dir, paste0(image_file, ".rds")))
    i <-  i + 1
    setTxtProgressBar(pb, value = i)
  }
close(pb)
  
}


###########################################
# Subject data
# ------------
# Here we're loading the subject data.
# The data that gets loaded is synthetic
# and has been generated from the original.
###########################################

sub_data <- read_csv(df_file, show_col_types = FALSE)

# This is the base subject data we will be working with
biofinder_df__ <- sub_data |> 
  mutate(
    fmri_date = as.Date(str_extract(csv_rsqa__index, "(?<=__).*"), format = "%Y%m%d"),
    visit = Visit,
    etiology = underlying_etiology_text_baseline_variable,
    image_file = csv_rsqa__index,
    tau_file = csv_tnic_sr_mr_fs__index,
    abnorm_ab = Abnormal_CSF_Ab42_Ab40_Ratio,
    sex = gender_baseline_variable,
    education = education_level_years_baseline_variable,
    diagnosis = diagnosis_baseline_variable,
    ab_ratio = CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022,
    a_syn = CSF_SynucleinAlpha_RTQuIC_2022,
    apoe4 = grepl("4", as.character(apoe_genotype_baseline_variable)),
    cho12 = tnic_cho_com_I_II,
    cho34 = tnic_cho_com_III_IV,
    cho56 = tnic_cho_com_V_VI,
    motion_filter = (rsqa__fd_max < 3 & rsqa__MeanFD < 0.3)
  ) |> 
  mutate(fmri_bl = fmri_date == min(fmri_date),
         has_longitudinal = n()>1, .by = "sid") |> 
  select(sid,  
         fmri_bl,
         fmri_date,
         visit,
         etiology,
         image_file,
         tau_file,
         has_longitudinal,
         age,
         education,
         diagnosis, 
         mPACC_v1, 
         apoe4,
         ab_ratio,
         abnorm_ab,
         a_syn,
         sex, 
         rsqa__MeanFD,
         rsqa__fd_max,
         motion_filter,
         cho12,
         cho34,
         cho56
  ) |> 
  group_by(sid) |> 
  filter(
    (any(abnorm_ab == 1, na.rm = TRUE) | !all(diagnosis == "MCI" & abnorm_ab == 0, na.rm = TRUE)) &
      any(diagnosis %in% c("AD", "MCI", "SCD", "Normal") | is.na(diagnosis), na.rm = TRUE)
  ) |>
  filter(any(diagnosis != "MSA", na.rm = TRUE)) |>
  filter(any(diagnosis != "PD", na.rm = TRUE)) |>
  filter(any(diagnosis != "PPA_NOS", na.rm = TRUE)) |>
  arrange(sid, fmri_date) |> 
  fill(diagnosis, .direction = "down") |> 
  ungroup() |> 
  filter(!(diagnosis %in% c("MCI", "AD") & abnorm_ab == 0) | is.na(diagnosis) | is.na(abnorm_ab))
  


# This is to get a dataframe with all subjects before filtering on motion
biof_motion_unfilt <- biofinder_df__ |> group_by(sid) |> 
  fill(sex, .direction = "downup") |>
  arrange(fmri_date) |>
  mutate(
    age_known = if_else(!is.na(age), age, NA_real_),
    date_known = if_else(!is.na(age), fmri_date, as.Date(NA))
  ) |>
  fill(age_known, date_known, .direction = "down") |>
  mutate(
    time_diff_years = as.numeric(difftime(fmri_date, date_known, units = "days")) / 365.25,
    age_estimated = age_known + time_diff_years,
    age = if_else(is.na(age), age_estimated, age)
  ) |>
  ungroup()

biofinder_df__ <- biofinder_df__ |> filter(motion_filter)

patvars <- biofinder_df__ |> filter(fmri_bl) |> 
  dplyr::select(sid,
                ab_ratio,
                cho12,
                cho34,
                cho56
  ) |>  drop_na()

withr::with_seed(12345, {
  space <- reduce_dimensionality(base::as.matrix(patvars |> select(-sid)), "euclidean")
  trajectory <- infer_trajectory(space)
})
traject <- trajectory$time
#biofinder_df__ |> left_join(patvars |> mutate(pathology_ad = traject) |> select(sid, pathology_ad)) |> ggplot(aes(pathology_ad, fill = diagnosis)) +geom_density(alpha = 0.5)
biofinder_df__ <- biofinder_df__ |> left_join(patvars |> mutate(pathology_ad = traject) |> select(sid, pathology_ad)) 

patvars_tau <- biofinder_df__ |> 
  dplyr::select(image_file,
                cho12,
                cho34,
                cho56
  ) |>  drop_na()
withr::with_seed(123456, {
  space <- reduce_dimensionality(base::as.matrix(patvars_tau |> select(-image_file)), "euclidean")
  trajectory <- infer_trajectory(space)
})
traject <- trajectory$time
#biofinder_df__ |> left_join(patvars_tau |> mutate(tau_path = traject) |> select(image_file, tau_path)) |> ggplot(aes(tau_path, fill = diagnosis)) +geom_density(alpha = 0.5)
biofinder_df <- biofinder_df__ |> left_join(patvars_tau |> mutate(tau_pathology = traject) |> select(image_file, tau_pathology)) |> 
  group_by(sid) |> 
  fill(sex, .direction = "downup") |>
  arrange(fmri_date) |>
  mutate(
    age_known = if_else(!is.na(age), age, NA_real_),
    date_known = if_else(!is.na(age), fmri_date, as.Date(NA))
  ) |>
  fill(age_known, date_known, .direction = "down") |>
  mutate(
    time_diff_years = as.numeric(difftime(fmri_date, date_known, units = "days")) / 365.25,
    age_estimated = age_known + time_diff_years,
    age = if_else(is.na(age), age_estimated, age)
  ) |>
  ungroup()
rm(biofinder_df__, space, trajectory, patvars_tau, patvars)


########  ########  ########          ######     ###    ##        ######  
##     ## ##     ## ##               ##    ##   ## ##   ##       ##    ## 
##     ## ##     ## ##               ##        ##   ##  ##       ##       
########  ########  ######   ####### ##       ##     ## ##       ##       
##        ##   ##   ##               ##       ######### ##       ##       
##        ##    ##  ##               ##    ## ##     ## ##       ##    ## 
##        ##     ## ########          ######  ##     ## ########  ######  

#########################################################
## Here we do some pre-calculations to get the 
## derivatives which we use to model later on
#########################################################

#########################################################
## Extract connectomes
#########################################################


if (from_start) {
  
  ts_files <- list.files(ts_dir, pattern = "\\.rds$", full.names = TRUE)
  # The below is much faster if using the .rds fileformat
  print("Reading timeseries")
  ts_list <- pblapply(ts_files, read_rds)

  names(ts_list) <- tools::file_path_sans_ext(basename(ts_files))
  
  # Take only the timeseries from subjects of interest
  timeseries <- ts_list[biofinder_df$image_file]
  
  dir.create(connectome_dir, showWarnings = FALSE)
  already_proc_conn <- list.files(connectome_dir) |> str_remove_all(".rds")
  proc_conn <- names(timeseries)[!(names(timeseries) %in% already_proc_conn)]
  
  print("Writing connectomes")
  pb = txtProgressBar(min = 0, 
                      max = ifelse(length(proc_conn) != 0, 
                                   length(proc_conn), 1),
                      initial = 0,
                      style = 3)
  i = 0
  for (img_file in proc_conn){
    tryCatch(
      {
        ts <- timeseries[[img_file]]
        ts <- scale(ts)
        colnames(ts) <- rois
        connectome <- cor(ts)
        if (sum(is.na(connectome)) > 0) stop("There are zero variance timeseries")
        write_rds(connectome, file.path(connectome_dir, paste0(img_file, ".rds")))
      },
      error = function(cond) {
        print(cond)
      },
      warning = function(cond) {
        print(cond)
      })
    i = i + 1
    setTxtProgressBar(pb,i)
  }
  close(pb)
}

# These are all subjects who successfully got calculated connectomes
success_vec <- list.files(connectome_dir) |> str_remove_all(".rds")

# Filter the data on that
biofinder_df <- biofinder_df |> filter(image_file %in% success_vec)
biof_motion_unfilt <- biof_motion_unfilt |> filter(image_file %in% success_vec)


print("Writing clean data")
dir.create(clean_dir, showWarnings = FALSE)
write_rds(biofinder_df, file.path(clean_dir, "biofinder_df.rds"))
write_rds(biof_motion_unfilt, file.path(clean_dir, "biofinder_motion_unfiltered.rds"))


if (from_start) {
  
  print("Reading connectomes")
  con_cube_bf <- array(dim = c(length(rois), length(rois), length(success_vec)), 
                       dimnames = list(rois, rois, success_vec))
  pb = txtProgressBar(min = 0, max = length(success_vec), style = 3)
  i = 0
  for (img_file in success_vec){
    tryCatch(
      {
        con_cube_bf[, , img_file] <- read_rds(file.path(connectome_dir, paste0(img_file, ".rds")))
      },
      error = function(cond) {
        print(cond)
      },
      warning = function(cond) {
        print(cond)
      })
    i = i + 1
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
}
#########################################################
## Calculate connectivity derivatives
#########################################################

source("src/util.R")

if (from_start) {
  bf_fc <- fc_strength(con_cube_bf, roi_names = rois, threshold = 0, replace = 0)
  bf_affinity <- get_affinity(con_cube_bf, 
                              similarity_method = "cosine", 
                              roi_names = rois,
                              threshold = 0.75, 
                              atlas = yeo_msk)
  bf_aff_nothres <- get_affinity(con_cube_bf, 
                                 similarity_method = "correlation",
                                 roi_names = rois, 
                                 threshold = 0, 
                                 atlas = yeo_msk)
  
  fc_measures_bf <- c(bf_fc, bf_affinity)
  write_rds(fc_measures_bf, file.path(clean_dir, "connectivity_derivatives_biofinder.rds"))
  write_rds(bf_aff_nothres, file.path(clean_dir, "bf_aff_no_thresh_corr.rds"))
} else {
  fc_measures_bf <- readRDS(file.path(clean_dir, "connectivity_derivatives_biofinder.rds"))
  bf_aff_nothres <- readRDS(file.path(clean_dir, "bf_aff_no_thresh_corr.rds"))
}

#########################################################
## Derive group gradients
########################################
#################

# Source file for calculations related to deriving gradients
source("src/util_gradients.R")
source("src/util_vis.R")

if (from_start) {
  
  gradient_dir <- "data/gradients"
  
  if (extract_timeseries) {
  # images have been accessed from https://neurovault.org/collections/1598/
  nifti_images <- list.files(gradient_dir, pattern = "\\.nii.gz$")
  
  require(reticulate)
  nilearn <- import("nilearn", convert = FALSE)
  schaefer <- nilearn$datasets$fetch_atlas_schaefer_2018(n_rois = 1000L, data_dir = "data/atlas_data/")
  atlas_file <- schaefer$maps
  masker <- nilearn$maskers$NiftiLabelsMasker(
    labels_img = atlas_file,
    smoothing_fwhm = as.integer(6),
    standardize = FALSE,
    verbose = 0L
  )
  
  marg_gradients <- matrix(ncol = length(nifti_images), nrow = length(rois))
  colnames(marg_gradients) <- c("gradient1", "gradient2", "gradient3", "gradient4", "gradient5")
  i <- 0
  pb <- txtProgressBar(min = i, max = length(nifti_images), style = 3)
  for (img in nifti_images){
    grad_val <- masker$fit_transform(file.path(gradient_dir, img))
    grad_val_r <- c(py_to_r(grad_val))
    gradient <- str_split_1(img, "\\.")[1]
    marg_gradients[, gradient] <- grad_val_r
    i <-  i + 1
    setTxtProgressBar(pb, value = i)
  }
  close(pb)
  marg_gradients <- as_tibble(marg_gradients)
  } else {
    marg_gradients <- read_csv("data/gradients/margulies_gradients.csv")
  }
  
  ###################################################
  # Create gradients from average healthy connectomes
  ###################################################
  
  healthy_young_connectomes <- con_cube_bf[, , biofinder_df |> filter(fmri_bl, age < 61, diagnosis == "Normal", abnorm_ab==0, !apoe4) |> pull(image_file)]
  average_connectome <- apply(healthy_young_connectomes, c(1, 2), mean)
  #write_rds(average_connectome, file.path(atlas_dir, "average_connectome_normalyoung.rds"))
  average_connectome = read_rds(file.path(atlas_dir, "average_connectome_normalyoung.rds"))
  rm(healthy_young_connectomes)
  
  # This takes some time, creates supplementary figure 9
  gradient_comparison <- plot_grads_over_params(connectome_list = list(biofinder=average_connectome))
  
  comp_plot <- ggdraw() +
    draw_plot(gradient_comparison[["biofinder"]], x = 0, y = 0, width = 1, height = 0.95) +
    draw_label(label = "PCA", x = 0.34, y = 0.985, size = 28) +
    draw_label(label = "DME", x = 0.785, y = 0.985, size = 28) +
    draw_label(label = "Threshold:", x = 0.34, y = 0.955, size = 18) +
    draw_label(label = "Threshold:", x = 0.785, y = 0.955, size = 18) +
    draw_line(x = c(0.15, 0.545), y = c(0.9725, 0.9725)) +
    draw_line(x = c(0.59, 0.985), y = c(0.9725, 0.9725)) 
  
  comp_plot <- comp_plot + theme(plot.background = element_rect(color = "black"))
  # This is max millimiter (of journals) in inches
  img_width = 180 / 25.4
  
  # Scale image to balance figure elements and fontsize, but need "shrink" image after
  # so the fontsize chosen should be large enough to shrink it by the scaling factor
  scaling_factor <-  3
  magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")
  p_name <- "gradient_param_comparison_bf.png"
  
  ggsave(file.path(figure_path, p_name), comp_plot , #patch_plots[["biofinder"]], 
         width = img_width*scaling_factor, height = img_width*scaling_factor*0.675, bg ="white")
  img <- magick::image_read(file.path(figure_path, p_name))
  img_resized <- magick::image_resize(img, magick_geom_scaling)
  magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)

  
  grad_list_bf <- get_gradients(connectome_ests = list(biofinder = average_connectome),
                                reference_gradients = marg_gradients,
                                n_gradients = c(1,2,3),
                                threshold = 0.0,
                                similarity_method = "cosine",
                                on_affinity = FALSE,
                                method = "pca",
                                visualize = FALSE)
  
  
  params <- expand_grid(method = c("pca", "diffusion"), affinity = c(FALSE, TRUE), sim_method = "cosine", threshold = c(0.0, 0.25, 0.5, 0.75)) |> 
    filter(!(method == "diffusion" & !affinity)) |> 
    filter(!(method == "pca" & affinity)) |> 
    mutate(sim_method = ifelse(!affinity, NA, sim_method))
  
  print("Calculating gradients over various parameters")
  gradient_data <- c()
  varexp_df <- c()
  for (i in 1:nrow(params)) {
    param_i <- params[i, ]
    grad_list <- get_gradients(connectome_ests = list(biofinder = average_connectome),
                               n_gradients = c(1,2,3),
                               threshold = param_i$threshold,
                               similarity_method = param_i$sim_method,
                               on_affinity = param_i$affinity,
                               method = param_i$method,
                               visualize = FALSE,
                               side_density = FALSE)
    
    gradient_data <- rbind(gradient_data, grad_list$gradients)
    varexp_df <- rbind(varexp_df, grad_list$varexp)
  }
  
  gradient_data <- gradient_data |> distinct()
  
}



# ########  ######## ########  ##       ####  ######     ###    ######## ####  #######  ##    ## 
# ##     ## ##       ##     ## ##        ##  ##    ##   ## ##      ##     ##  ##     ## ###   ## 
# ##     ## ##       ##     ## ##        ##  ##        ##   ##     ##     ##  ##     ## ####  ## 
# ########  ######   ########  ##        ##  ##       ##     ##    ##     ##  ##     ## ## ## ## 
# ##   ##   ##       ##        ##        ##  ##       #########    ##     ##  ##     ## ##  #### 
# ##    ##  ##       ##        ##        ##  ##    ## ##     ##    ##     ##  ##     ## ##   ### 
# ##     ## ######## ##        ######## ####  ######  ##     ##    ##    ####  #######  ##    ## 


adni_df <- read_csv(df_file_adni, show_col_types = FALSE)
adni_df_ <- adni_df |> select(ID, file_func, DX, age, sex, amyloid_status, centiloid,
                              alpha_syn = CSF_AlphaSyn_seeding,
                              education_yrs,
                              APOE4_alleles, EXAMDATE_func, 
                              ABETA42, ABETA40, contains("braak")) 
# This is just a workaround for ID handling, 
# which differs slightly from the synthetic 
if (real_data) {
  id_ses <- c()
  for (i in adni_df_$file_func) {
    id_ses <- c(id_ses, str_split_1(str_split_i(i, "/", 10), "_")[1:2] %>% paste0(collapse = "_") )
  }
  adni_df_ <- adni_df_ |> mutate(file_func = id_ses)
}

rm(adni_df)

adni_df__ <- adni_df_ |> 
  mutate(EXAMDATE_func = as.Date(EXAMDATE_func, format = "%m-%d-%y"),
         id_ses = file_func) |> 
  relocate(id_ses) |> 
  mutate(fmri_bl = EXAMDATE_func == min(EXAMDATE_func), .by = ID) |> 
  mutate(
    apoe4 = APOE4_alleles > 0,
    abnorm_ab = ifelse(amyloid_status == "Ab.neg", 0, 1),
    diagnosis = DX,
    fmri_date = EXAMDATE_func,
    braak1 = tau.SUVR.DK.braak1,
    braak34 = tau.SUVR.DK.braak34,
    braak56 = tau.SUVR.DK.braak56,
    ab_ratio = ABETA42/ABETA40,
    sex = ifelse(sex == "male", 0, 1)
  ) 

patvars <- adni_df__ |> 
  dplyr::select(id_ses,
                ab_ratio,
                braak1,
                braak34,
                braak56
  ) |>  drop_na()

withr::with_seed(123456, {
  space <- reduce_dimensionality(base::as.matrix(patvars |> select(-id_ses)), "euclidean")
  trajectory <- infer_trajectory(space)
})
traject <- trajectory$time

adni_df___ <- adni_df__ |> left_join(patvars |> mutate(pathology_ad = traject) |> select(id_ses, pathology_ad))

patvars_tau <- adni_df__ |> 
  dplyr::select(id_ses,
                braak1,
                braak34,
                braak56
  ) |>  drop_na()

withr::with_seed(123456, {
  space <- reduce_dimensionality(base::as.matrix(patvars_tau |> select(-id_ses)), "euclidean")
  trajectory_tau <- infer_trajectory(space)
})

traject_tau <- trajectory_tau$time
adni_df___ <- adni_df___ |> left_join(patvars_tau |> mutate(tau_pathology = traject_tau) |> select(id_ses, tau_pathology))

rm(adni_df_, adni_df__)

############
# Timeseries
############

fd_files <- list.files(ts_dir_adni, pattern = "displacement\\.csv")

adni_fd <- list()
for (file in fd_files) {
  sub_and_ses <- tools::file_path_sans_ext(file) |> str_remove_all("_framewise_displacement") |> 
    str_remove_all("_task-rest_bold_st_mcflirt_bp_det_24HMP8Phys_4mmFWHM")
  adni_fd[[sub_and_ses]] <- read_csv(file.path(ts_dir_adni, file), show_col_types = FALSE) 
  if ("V1" %in% names(adni_fd[[sub_and_ses]])) {
    adni_fd[[sub_and_ses]] <- adni_fd[[sub_and_ses]] |> 
      rename(displacement = V1)
  }
}

rsqa <- lapply(adni_fd, function(x) c(rsqa__MeanFD = x |> pull(displacement) |> mean(), rsqa__MaxFD = x |> pull(displacement) |> max()))
rsqa_fd <- do.call(rbind, rsqa) |> as_tibble(rownames = NA) |> rownames_to_column("id_ses")


if (from_start) {
  
  timeseries_files <- list.files(ts_dir_adni, pattern = "timeseries\\.csv")
  
  adni_timeseries <- pblapply(timeseries_files, function(file) {
    read_csv(file.path(ts_dir_adni, file), progress = FALSE, show_col_types = FALSE) |> 
      as.matrix()
  })
  
  names(adni_timeseries) <- timeseries_files |> tools::file_path_sans_ext() |>  str_remove_all("_timeseries")
  
  # This is to scrub timeseries
  set_false_window <- function(log_vec) {
    result <- log_vec
    for (i in seq_along(log_vec)) {
      if (!log_vec[i]) {
        if (i > 1) result[i - 1] <- FALSE
        result[i] <- FALSE
        if (i + 1 <= length(log_vec)) result[i + 1] <- FALSE
        if (i + 2 <= length(log_vec)) result[i + 2] <- FALSE
      }
    }
    return(result)
  }
  
  scrubbed_time_series <- list()
  for(file in names(adni_timeseries)) {
    
    fd <- c(0, adni_fd[[file]])
    # Setting the same threshold as originally done by Franzmeier
    fd_filt <- set_false_window(fd<0.5)
    scrubbed_time_series[[file]] <- adni_timeseries[[file]][fd_filt, ]
    
  }
  frame_length <- sapply(scrubbed_time_series, nrow)
  # Uncomment this line if running with real data
  # timeseries <- scrubbed_time_series[frame_length>100]
   timeseries <- scrubbed_time_series
}


##############
# Connectomes
##############
if (from_start) {
 
  dir.create(connectome_dir_adni, showWarnings = FALSE)
  already_proc_conn <- list.files(connectome_dir_adni) |> tools::file_path_sans_ext()
  proc_conn <- names(timeseries)[!(names(timeseries) %in% already_proc_conn)]
  
  print("Writing connectomes ADNI")
  pb = txtProgressBar(min = 0, max = ifelse(length(proc_conn) == 0, 1, length(proc_conn)), style = 3)
  i = 0
  for (img_file in proc_conn){
    tryCatch(
      {
        ts <- timeseries[[img_file]]
        ts <- scale(ts)
        colnames(ts) <- rois
        connectome <- cor(ts)
        if (sum(is.na(connectome)) > 0) stop("There are zero variance timeseries")
        write_rds(connectome, file.path(connectome_dir_adni, paste0(img_file, ".rds")))
      },
      error = function(cond) {
        print(cond)
      },
      warning = function(cond) {
        print(cond)
      })
    i = i + 1
    setTxtProgressBar(pb,i)
  }
  close(pb)
}


success_vec <- list.files(connectome_dir_adni) |> tools::file_path_sans_ext()
adni_df <- adni_df___ |> filter(id_ses %in% success_vec) |> 
  inner_join(rsqa_fd, join_by(id_ses==id_ses)) # Set the join variable depending on if longitudinal data is used

adni_df_unfilt <- adni_df |> mutate(motion_filter = (rsqa__MeanFD<0.3 & rsqa__MaxFD<3))
adni_df <- adni_df_unfilt |>  filter(motion_filter)


write_rds(adni_df, file.path(clean_dir, "adni_df.rds"))
write_rds(adni_df_unfilt, file.path(clean_dir,"adni_motion_unfilt.rds"))

if (from_start) {
  print("Reading connectomes")
  con_cube_adni <- array(dim = c(length(rois), 
                                 length(rois), 
                                 length(success_vec)), 
                         dimnames = list(rois, rois, success_vec))
  pb = txtProgressBar(min = 0, max = length(success_vec), initial = 0, style = 3)
  i = 0
  for (img_file in success_vec){
    tryCatch(
      {
        con_cube_adni[, , img_file] <- read_rds(file.path(connectome_dir_adni, paste0(img_file, ".rds")))
      },
      error = function(cond) {
        print(cond)
      },
      warning = function(cond) {
        print(cond)
      })
    i = i + 1
    setTxtProgressBar(pb,i)
  }
  close(pb)
}

#########################################################
## Calculate connectivity derivatives
#########################################################

source("src/util.R")

if (from_start) {
  print("Calculating FC derivatives")
  # Calculate connectivity strenght
  adni_fc <- fc_strength(con_cube_adni, 
                         roi_names = rois, 
                         threshold = 0, 
                         replace = 0)
  
  # Calculate connectivity affinity
  adni_affinity <- get_affinity(con_cube_adni, 
                                similarity_method = "cosine",
                                roi_names = rois, 
                                threshold = 0.75, 
                                atlas = yeo_msk)
  
  # Calculate affinity using different method
  adni_affinity_no_thresh <- get_affinity(con_cube_adni, 
                                          similarity_method = "correlation", 
                                          roi_names = rois, 
                                          threshold = 0, 
                                          atlas = yeo_msk)
  
  fc_measures_adni <- c(adni_fc, adni_affinity)
  write_rds(fc_measures_adni, file.path(clean_dir,"connectivity_derivatives_adni.rds"))
  write_rds(adni_affinity_no_thresh, file.path(clean_dir,"adni_affinity_no_thresh_corr.rds"))
}

fc_measures_adni <- readRDS(file.path(clean_dir, "connectivity_derivatives_adni.rds"))
adni_affinity_no_thresh <- readRDS(file.path(clean_dir, "adni_affinity_no_thresh_corr.rds"))
#########################################################
## Derive group gradients
#########################################################

source("src/util_gradients.R")

if (from_start) {
  # Create connectome from young healthy individuals
  yh_filt <- adni_df |> filter(fmri_bl, abnorm_ab==0, !apoe4, DX == "CN") |> pull(file_func)
  healthy_young_connectomes <- con_cube_adni[, , yh_filt]
  average_connectome_adni <- apply(healthy_young_connectomes, c(1, 2), mean)
  rm(healthy_young_connectomes)
  
  # if you want to compare ADNI gradients to margulies, uncomment:
  # gradient_comparison <- plot_grads_over_params(connectome_list = list(adni=average_connectome))
  
  grad_list_adni <- get_gradients(connectome_ests = list(adni = average_connectome_adni),
                                  reference_gradients = marg_gradients,
                                  n_gradients = c(1,2,3),
                                  threshold = 0.0,
                                  similarity_method = "cosine",
                                  on_affinity = FALSE,
                                  method = "pca",
                                  visualize = FALSE,
                                  side_density = TRUE)
  
  grad_df <- grad_list_bf$gradients |> bind_rows(grad_list_adni$gradients |> filter(study=="adni") )
  write_rds(grad_df, file.path(clean_dir, "grad_df.rds"))
  
  gradient_data_adni <- c()
  varexp_df <- c()
  for (i in 1:nrow(params)) {
    param_i <- params[i, ]
    grad_list <- get_gradients(connectome_ests = list(adni = average_connectome),
                               n_gradients = c(1,2,3),
                               threshold = param_i$threshold,
                               similarity_method = param_i$sim_method,
                               on_affinity = param_i$affinity,
                               method = param_i$method,
                               visualize = FALSE,
                               side_density = FALSE)
    
    gradient_data_adni <- rbind(gradient_data_adni, grad_list$gradients)
    varexp_df <- rbind(varexp_df, grad_list$varexp)
  }
  
  gradient_data_adni <- gradient_data_adni |> distinct()
  write_rds(rbind(gradient_data, gradient_data_adni) |> distinct(), file.path(clean_dir, "gradients_over_params.rds"))
}

gradient_over_params <- readRDS(file.path(clean_dir, "gradients_over_params.rds"))
grad_df <- readRDS(file.path(clean_dir, "grad_df.rds"))

##     ##  #######  ########  ######## ##       ##       #### ##    ##  ######   
###   ### ##     ## ##     ## ##       ##       ##        ##  ###   ## ##    ##  
#### #### ##     ## ##     ## ##       ##       ##        ##  ####  ## ##        
## ### ## ##     ## ##     ## ######   ##       ##        ##  ## ## ## ##   #### 
##     ## ##     ## ##     ## ##       ##       ##        ##  ##  #### ##    ##  
##     ## ##     ## ##     ## ##       ##       ##        ##  ##   ### ##    ##  
##     ##  #######  ########  ######## ######## ######## #### ##    ##  ######   

########################
# Methods plot 
########################

# Source the methods_figure.R script to generate all plots for the methods figure.
# However, the figures have then been manually put together in inkscape 

# source("src/methods_figure.R")

########################
# Main figures
########################

# This script contains an array of functions that creates the figures in the paper
source("src/util_vis.R")

# These parameters set the final figure width, and the scaling to use to 
# get elements and font sizes aligned
img_width <-  180 
scaling_factor <-  1

################
# FIGURE 1
################


fig_one <- figure_one(subject_data = biofinder_df,
                      subject_data_replication = adni_df |> filter(fmri_bl),
                      measures_list =list(nodal_affinity = fc_measures_bf$affinity), 
                      measures_list_replication = list(nodal_affinity = fc_measures_adni$affinity),
                      gradients_df = grad_df |> filter(study == "biofinder"),
                      gradients_df_replication = grad_df |> filter(study == "adni"),
                      shade = TRUE,
                      shade_a = 0.025,
                      shade_s = 0.05,
                      tag_size = 10,
                      p_size = 0.1,
                      r2_sizing1 = 0.9,
                      raster = TRUE,
                      ggrastr_dpi = 300,
                      selected_gradients = c(1, 3), 
                      empt_row_height = -0.15,
                      brain_title_size = rel(1),
                      brain_title_vjust = -1.5,
                      b_size = 7,
                      draw_size = 7,
                      plot_title_size = 0.9,
                      axes_title_size = 0.9,
                      ax_txt_size = rel(0.85),
                      boxed = TRUE,
                      split = TRUE,
                      fig = 1,
                      fig1_formula_bf = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"),
                      fig1_formula_ad = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"),
                      brain_plot_names_f1 = c(NA, "AD Pathology"))


p_name <- "figure1.pdf"
ggsave(file.path(figure_path, p_name), fig_one[[1]] |> pad_plot(),
       width = img_width*scaling_factor, height = img_width*0.5*scaling_factor, 
       units = "mm", dpi = 300, device = "pdf", bg = "transparent")

source_data <- fig_one$tmaps$bf |> mutate(study = "biofinder") |> 
  bind_rows(fig_one$tmaps$adni |> mutate(study = "adni")) |> 
  select(-n, -model_formula) |> 
  left_join(grad_df |> filter(study %in% c("adni", "biofinder")) |> 
              select(region, study, starts_with("gradient")))
write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))


###########################
# FIGURE 2 (non-linear)
###########################

source("src/util.R")
source("src/util_vis.R")

print("Running non-linear models")
  
gam_preds_affinity <- gam_pred_nodes(biofinder_df |> filter(fmri_bl),
                                     fc_matrix = fc_measures_bf$affinity,
                                     roi_names = rois, 
                                     id_var = "image_file",
                                     print_ticker = TRUE,
                                     model_formula = formula("FC ~ s(age) + s(pathology_ad) + sex + rsqa__MeanFD"))


source("src/util_vis.R")
nonlin_p <- plot_gams_v1(
  gam_predictions = gam_preds_affinity,
  grad_df = grad_df |> filter(study == "biofinder"),
  biofinder_data = biofinder_df |> filter(fmri_bl),
  b_size = 7,
  corr_label_size = 1.2,
  legend_key_height = 0.4,
  legend_key_width = 0.01,
  pathology_legend_width = 0.10,
  pathology_legend_height = 0.06,
  pathology_label_size = 1.75,
  net_legend_text_rel = 0.5,
  net_legend_point_size = 0.75,
  net_legend_width = 0.07,
  net_legend_height = 0.04,
  add_shade = TRUE,
  shade_alpha = 0.1,
  shade_size = 0.1,
  point_size = 0.075,
  point_alpha = 0.5,
  struct_linewidth = 0.2,
  plot_linewidth = 0.3,
  scatter_raster_dpi = 300,
  rasterize_brains = TRUE,
  rasterize_scatters = TRUE,
  canvas_text_size = 5
)


p_name <- "figure2.pdf"
ggsave(
  file.path(figure_path, p_name ),
  nonlin_p$plot |> pad_plot(),
  width = 88,
  height = 88*1.3,
  units = "mm",
  device = "pdf",
  bg = "white"
)

# Write source data
write_csv(nonlin_p$source_data$pathology_plot |> 
            mutate(value = as.vector(value)), file.path(source_figure_path, paste0("figure2_pathology_plot", ".csv")))
write_csv(nonlin_p$source_data$quartile_fcs_slopes, file.path(source_figure_path, paste0("figure2_quartile_fcs_slopes", ".csv")))
write_csv(nonlin_p$source_data$ventile_mean_predicted_fcs, file.path(source_figure_path, paste0("figure2_ventile_mean_predicted_fcs", ".csv")))
write_csv(nonlin_p$source_data$gradient_correlations, file.path(source_figure_path, paste0("figure2_gradient_correlations", ".csv")))


#######################################
# FIGURE 3 (Longitudinal and windowing)
#######################################

source("src/util_vis.R")

# Create the longitudinal dataframe
biofinder_df |> 
  drop_na(age, tau_pathology, sex, 
          rsqa__MeanFD) |> 
  group_by(sid) |> 
  arrange(sid, fmri_date) |> 
  mutate(
    diag_bl = first(diagnosis),
    ab_bl = first(abnorm_ab),
    age_bl = min(age),
    path_bl = first(tau_pathology),
    mPACC_bl = first(mPACC_v1),
    ageΔ  = age - age_bl,
    pathΔ = tau_pathology - path_bl,
    pathology_mc = tau_pathology - mean(tau_pathology),
    path = tau_pathology,
    mPACCΔ = mPACC_v1 - mPACC_bl,
    time = difftime(fmri_date, first(fmri_date), units = "days")/365
  ) |>  
  ungroup() |>   
  filter(n()>1, .by = "sid") -> long_bf_


roi_fc <- fc_measures_bf$affinity[, rois[30]] %>% 
  enframe("image_file", "FC")
reg_df <- long_bf_ %>% 
  inner_join(roi_fc, by = "image_file")

library(performance)
library(lme4)
fit <- lmer(formula("FC ~ time  + age_bl + path_bl + pathΔ + sex + rsqa__MeanFD + (1 | sid)"), data = reg_df)
check_collinearity(fit)

bf_longitudinal <-  plot_gradient_relationships(long_bf_ |> mutate(yearly_path = pathΔ/as.numeric(time)), 
                                                gradient_data = grad_df %>% filter(study=="biofinder"), 
                                                gradients = c(1, 3),
                                                empty_row_height = -0.1,
                                                gradient_colors = gradient_cols,
                                                base_size_ = 7,
                                                rasterize = TRUE,
                                                ggrastr_dpi = 300,
                                                add_shade = TRUE, 
                                                shade_alpha = 0.0075,
                                                shade_size = 0.001,
                                                r2_size = rel(3.8),
                                                list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                covariates = c("sex", "rsqa__MeanFD"),
                                                filter_criteria = quo(),
                                                show_networks = FALSE,
                                                plt_title = NULL,
                                                tag_prefix = "",
                                                tag_sep = "",
                                                layout_construction = "horizontal",
                                                include_gradient_plots = TRUE,
                                                right_term_side = FALSE,
                                                cache_runs = FALSE,
                                                longitudinal = TRUE,
                                                sub_id = "sid",
                                                longitudinal_formula = formula(paste0("FC ~ age_bl + path_bl + yearly_path + sex + rsqa__MeanFD + (1 | sid)")))


# Very case specific function for creating figure 3

longitudindal_figs <- longitudinal_and_window_analysis(
  long_df = long_bf_,
  b_size = 7,
  run_windowing = 
    ifelse(from_start & 
             !file.exists("data/processed_and_cleaned/sliding_window_age_res.rds") &
             !file.exists("data/processed_and_cleaned/sliding_window_path_res.rds"), 
           TRUE, FALSE) 
)

source("src/util_vis.R")
analysis_data <- longitudindal_figs$analysis_data
longitudindal_figs2 <- plot_longitudinal_and_window_analysis(
  analysis_data = analysis_data,
  b_size = 7,
  label_size = 10,
  title_size = 7,
  annotation_size = 3,
  strip_size_rel = 1.2,
  axis_title_rel = 1.2,
  subtitle_size = 5
)

img_width = 180
scaling_factor <- 1

p_name <- "figure3.pdf"
ggsave(file.path(figure_path, p_name), longitudindal_figs2[["main_fig"]] |> pad_plot(),
       width = img_width*scaling_factor, height = img_width*1*scaling_factor, units = "mm", dpi = 300, device = "pdf")

write_csv(analysis_data$bf_longitudinal$tmaps |> 
            inner_join(grad_df |> filter(study == "biofinder") |> 
                         select(region, starts_with("gradient"))), 
          file.path(source_figure_path, paste0("figure3_longitudinal_estimates", ".csv")))

write_csv(analysis_data$hist_map_data, 
          file.path(source_figure_path, paste0("figure3_2D_histogram", ".csv")))

write_csv(analysis_data$grad_cor_age, 
          file.path(source_figure_path, paste0("figure3_windowing_age", ".csv")))

write_csv(analysis_data$grad_cor_path, 
          file.path(source_figure_path, paste0("figure3_windowing_path", ".csv")))


p_name <- "supplementary_windowing.pdf"
ggsave(file.path(figure_path, p_name), longitudindal_figs[["supp_fig"]] |> pad_plot(),
       width = img_width*scaling_factor, height = img_width*0.5*scaling_factor, units = "mm", dpi = 300, device = "pdf", bg = "white")



######################
# FIGURE 4 (cognition)
######################

source("src/util_vis.R")

img_width = 180 
scaling_factor <-  1

fig4 <- figure_one(subject_data = biofinder_df,
                   subject_data_replication = adni_df |> filter(fmri_bl),
                   measures_list =list(nodal_affinity = fc_measures_bf$affinity), 
                   measures_list_replication = list(nodal_affinity = fc_measures_adni$affinity),
                   gradients_df = grad_df |> filter(study == "biofinder"),
                   gradients_df_replication = grad_df |> filter(study == "adni"),
                   selected_gradients = c(1, 3), 
                   empt_row_height = -0.25,
                   shade = TRUE,
                   shade_a = 0.1,
                   shade_s = 0.1,
                   tag_size = 9,
                   p_size = 0.01,
                   raster = TRUE,
                   ggrastr_dpi = 300,
                   b_size = 7,
                   draw_size = 7,
                   plot_title_size = 0.85,
                   axes_title_size = 0.85,
                   ax_txt_size = rel(0.75),
                   r2_sizing2 = 0.65,
                   boxed = FALSE,
                   split = TRUE,
                   fig = 2)

p_name <- "figure4.pdf"
ggsave(file.path(figure_path, p_name), fig4[[1]] |> pad_plot(), 
       width = img_width*scaling_factor, 
       height = img_width*0.4*scaling_factor,
       units = "mm", dpi = 300, 
       device = "pdf", bg = "transparent")


source_data <- fig4$tmaps$health |> mutate(plot = "healthy") |> 
  bind_rows(fig4$tmaps$clin |> mutate(plot = "clinical")) |> 
  select(-n, -model_formula) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient")))
write_csv(source_data, file.path(source_figure_path, paste0("figure4_estimates", ".csv")))



##############################################
# Cognition (imputation) and cognitive domains
##############################################

library(psych)
library(mice)


nic <- read_delim("stuff_for_revisions/nicola__20240326_143046.csv") 


cog <- nic |> select(1:5, mPACC_v1:fcsrt_immediate,
                     cognitive_test_date
)



cog_bl <- cog |> 
  filter(!is.na(cognitive_test_date)) |> 
  select(
    sid,
    cognitive_test_date,
    adas_delayed_word_recall,
    adas_immediate_word_recall_average,
    symbol_digit,
    trailmaking_a, 
    trailmaking_b,
    letter_s,
    animal_fluency,
    vosp_cube,
    bnt_15_2,
    aqt_color_form,
  ) |> 
  filter(cognitive_test_date == min(cognitive_test_date), .by = "sid") 

df <- cog_bl %>%
  mutate(cognitive_test_date = as.Date(cognitive_test_date))

baseline_dates <- df %>%
  group_by(sid) %>%
  summarise(
    baseline_date = min(cognitive_test_date, na.rm = TRUE),
    .groups = "drop"
  )

df <- df %>%
  left_join(baseline_dates, by = "sid") %>%
  mutate(
    days_from_baseline = as.numeric(cognitive_test_date - baseline_date)
  )


get_baseline_within_window <- function(date_diff, values, window_days = 365) {
  
  # index of baseline row
  i0 <- which(date_diff == 0)
  
  # if baseline is observed, use it
  if (length(i0) == 1 && !is.na(values[i0])) {
    return(values[i0])
  }
  
  # otherwise look within window
  candidates <- which(
    !is.na(values) &
      abs(date_diff) <= window_days
  )
  
  if (length(candidates) == 0) {
    return(NA_real_)
  }
  
  # take closest in time
  closest <- candidates[which.min(abs(date_diff[candidates]))]
  values[closest]
}


cog_vars <- c(
  "adas_delayed_word_recall",
  "adas_immediate_word_recall_average",
  "symbol_digit",
  "trailmaking_a",
  "trailmaking_b",
  "letter_s",
  "animal_fluency",
  "vosp_cube",
  "bnt_15_2",
  "aqt_color_form"
)


baseline_df <- df %>%
  group_by(sid) %>%
  summarise(
    baseline_date = first(baseline_date),
    
    across(
      all_of(cog_vars),
      ~ get_baseline_within_window(
        date_diff = days_from_baseline,
        values    = .x,
        window_days = 365
      ),
      .names = "{.col}"
    ),
    .groups = "drop"
  )

biofinder_cog <- baseline_df |> 
  inner_join(biofinder_df |> filter(fmri_bl), by = "sid") 

imp_df <- biofinder_cog %>%
  select(
    sid,
    age,
    sex,
    pathology_ad,
    education,
    diagnosis,
    mPACC_v1,
    any_of(cog_vars)
  ) |> 
  mutate(n_missing = rowSums(is.na(across(all_of(cog_vars))))) |> 
  filter(n_missing < 5) |> 
  select(-n_missing)


id_vars   <- "sid"

library(mice)


ini  <- mice(imp_df, maxit = 0, print = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix


meth[id_vars] <- ""
pred[id_vars, ] <- 0
pred[, id_vars] <- 0


no_impute_vars <- c(
  "sex",
  "diagnosis",
  "age",
  "education",
  "mPACC_v1",
  "pathology_ad"
)

pred[no_impute_vars, ] <- 0


cog_predictors <- c(
  "age",
  "sex",
  "education",
  "pathology_ad",
  "mPACC_v1"
)

pred[cog_vars, cog_predictors] <- 1
pred[cog_vars, cog_vars]       <- 1


diag(pred) <- 0


imp <- mice(
  imp_df,
  m      = 20,
  maxit  = 30,
  method = meth,
  predictorMatrix = pred,
  seed   = 42,
  print  = FALSE
)
imp_list <- complete(imp, action = "all")

ave_imp <- imp_list[[1]] |> column_to_rownames("sid") |> select(all_of(cog_vars))
ave_imp[, ] <- 0
for(i in seq_along(imp_list)) {
  ave_imp <- ave_imp + imp_list[[i]] |> column_to_rownames("sid") |> select(all_of(cog_vars))
}
ave_imp <- ave_imp/30
ave_imp <- ave_imp |> rownames_to_column("sid")


#####################
# Cognition modelling
#####################

cog_mat_cu <- ave_imp |> 
  select(
    sid,
    adas_delayed_word_recall,
    adas_immediate_word_recall_average,
    symbol_digit,
    trailmaking_a, 
    trailmaking_b,
    letter_s,
    animal_fluency,
    vosp_cube,
    bnt_15_2,
    aqt_color_form,
  ) |> inner_join(biofinder_cog |> select(-all_of(cog_vars))) |> 
  drop_na(pathology_ad) |>
  filter(diagnosis %in% c("Normal"), age <70 ) |> 
  select(
    sid,
    adas_delayed_word_recall,
    adas_immediate_word_recall_average,
    symbol_digit,
    trailmaking_a, 
    trailmaking_b,
    letter_s,
    animal_fluency,
    vosp_cube,
    bnt_15_2,
    aqt_color_form,
  ) |> 
  mutate(
    trailmaking_a = -trailmaking_a,
    trailmaking_b = -trailmaking_b,
    adas_delayed_word_recall = -adas_delayed_word_recall,
    adas_immediate_word_recall_average = -adas_immediate_word_recall_average,
    aqt_color_form = -aqt_color_form
  ) |> 
  column_to_rownames("sid")


cog_fac <- psych::omega(cog_mat_cu, nfactors = 4)


om_factors <- cog_fac[["scores"]]
refs <- cog_mat_cu |> select(adas_delayed_word_recall, vosp_cube, animal_fluency, trailmaking_a)


rename_factors_by_corr <- function(omega_factors, ref_vars) {
  global_cog <- omega_factors[, 1]
  omega_factors <- omega_factors[, -1]
  label_map <- c(
    adas_delayed_word_recall = "memory",
    trailmaking_a           = "executive",
    animal_fluency           = "language",
    vosp_cube                = "visuospatial"
  )
  cor_mat <- cor(ref_vars, omega_factors, use = "pairwise.complete.obs")
  strongest_ref <- apply(abs(cor_mat), 2, function(x) names(x)[which.max(x)])
  new_names <- label_map[strongest_ref]
  colnames(omega_factors) <- new_names
  omega_factors <- cbind(global_cog, omega_factors)
  return(omega_factors)
}



cog_mat <- ave_imp |> 
  select(
    sid,
    adas_delayed_word_recall,
    adas_immediate_word_recall_average,
    symbol_digit,
    trailmaking_a, 
    trailmaking_b,
    letter_s,
    animal_fluency,
    vosp_cube,
    bnt_15_2,
    aqt_color_form,
  ) |> inner_join(biofinder_cog |> select(-all_of(cog_vars))) |> 
  drop_na(pathology_ad) |>
  select(
    sid,
    adas_delayed_word_recall,
    adas_immediate_word_recall_average,
    symbol_digit,
    trailmaking_a, 
    trailmaking_b,
    letter_s,
    animal_fluency,
    vosp_cube,
    bnt_15_2,
    aqt_color_form,
  ) |> 
  mutate(
    trailmaking_a = -trailmaking_a,
    trailmaking_b = -trailmaking_b,
    adas_delayed_word_recall = -adas_delayed_word_recall,
    adas_immediate_word_recall_average = -adas_immediate_word_recall_average,
    aqt_color_form = -aqt_color_form
  ) |> 
  column_to_rownames("sid")


domain_scores <- psych::factor.scores(cog_mat, f = cog_fac$schmid$sl, method = "regression") 
domain_scores <- domain_scores$scores[, 1:5]
colnames(domain_scores) <- c("global_cog", "executive", "memory", "language", "visuospatial")


biofinder_cog_comp <- biofinder_cog |> 
  inner_join(domain_scores |> as_tibble(rownames = NA) |> rownames_to_column("sid")) 


cog_diff_plot <- function(b_size = 7, subt_size = rel(0.75)) {
  fit <- biofinder_cog_comp |> 
    mutate(
      # memory = resid(lm(memory ~ global_cog)),
      # executive = resid(lm(executive ~ global_cog))
    ) |> 
    lm(memory ~ scale(age) + scale(pathology_ad) + sex , data = _) 
  
  # Extract coefficients
  coef_df <- broom::tidy(fit) |> 
    filter(term %in% c("scale(age)", "scale(pathology_ad)")) |> 
    select(term, estimate, std.error) |> 
    mutate(ci_low = estimate - 1.96*std.error,
           ci_high = estimate + 1.96*std.error,
           outcome = "memory",
           term = recode(term,
                         "scale(age)" = "Age",
                         "scale(pathology_ad)" = "Pathology"))
  
  # Plot
  mem_plot <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    labs(y = "Memory", x = "") +
    theme_bw(base_size = b_size) 
  
  
  
  fit_exec <- biofinder_cog_comp |> 
    mutate(
      # memory = resid(lm(memory ~ global_cog)),
      # executive = resid(lm(executive ~ global_cog))
    ) |> 
    lm(executive ~ scale(age) + scale(pathology_ad) + sex, data = _) 
  
  # Extract coefficients
  coef_df <- broom::tidy(fit_exec) |> 
    filter(term %in% c("scale(age)", "scale(pathology_ad)")) |> 
    select(term, estimate, std.error) |> 
    mutate(ci_low = estimate - 1.96*std.error,
           ci_high = estimate + 1.96*std.error,
           outcome = "executive",
           term = recode(term,
                         "scale(age)" = "Age",
                         "scale(pathology_ad)" = "Pathology"))
  
  source_data <- bind_rows(
    mem_plot$data,
    coef_df
  )
  
  # Plot
  exec_plot <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    labs(y = "Executive", x = "") +
    theme_bw(base_size = b_size) 
  
  n = length(fit_exec$residuals)
  
  library(patchwork)
  cog_eff_plot <- wrap_plots(list(mem_plot, exec_plot), ncol = 1, guides = "collect", axes = "collect")  +
    plot_annotation(
      title = paste0("Effects on memory/executive function (n = ", n, ")"),
      subtitle = "Memory ~ scale(age) + scale(AD pathology) + sex \nExecutive ~ scale(age) + scale(AD pathology) + sex",
      theme = theme(plot.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
                    plot.title = element_text(size = rel(0.95), margin = margin(l = 8), hjust = 0),
                    plot.subtitle = element_text(family = "mono", margin = margin(l = 8, t = 3), size = subt_size, hjust = 0)
                    )) &
    theme(axis.text.x = element_text(size = rel(1.5)),
          axis.title.y = element_text(size = rel(1.15))
          )
  list(
    plot = cog_eff_plot,
    data = source_data
  )
}

cog_diff_res <- cog_diff_plot(subt_size = rel(0.75))
cog_diff <- cog_diff_res$plot



cognitive_composites <- plot_gradient_relationships(biofinder_cog_comp |> 
                                                      filter(fmri_bl), 
                                                    gradient_data = grad_df |> filter(study== "biofinder"), 
                                                    gradients = c(1, 3),
                                                    gradient_colors = gradient_cols,
                                                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                    brain_title_size = 6,
                                                    title_lmargin = 10,
                                                    axis_text_size = 5,
                                                    axis_title_size = 6,
                                                    scatter_title_vjust = 0,
                                                    rasterize = TRUE,
                                                    ggrastr_dpi = 300,
                                                    add_shade = TRUE, 
                                                    shade_alpha = 0.1,
                                                    shade_size = 0.1,
                                                    point_size = 0.05,
                                                    point_alpha = 0.2,
                                                    l_width = 0.05,
                                                    empty_row_height = -0.225,
                                                    base_size_ = 7,
                                                    brain_title_vjust = 1,
                                                    r_spin_size = 0.75,
                                                    vect = TRUE,
                                                    mod_formula = formula(paste0(" ~ age + pathology_ad + global_cog + executive + memory + sex + rsqa__MeanFD")),
                                                    covariates = c("age", "pathology_ad", "sex", "rsqa__MeanFD"),
                                                    plt_title = "Higher executive or memory function show positive gradient alignment",
                                                    brain_names = c("Global Cognition", "Executive", "Memory"),
                                                    plt_subtitle = TRUE,
                                                    group_n_title = FALSE,
                                                    rectangle = TRUE,
                                                    plot_net_legend = TRUE,
                                                    net_legend_x = 0.15,
                                                    net_legend_y = 0.01)




source("src/make_omega_plot.R")
omega_plot_res <- plot_cognition_factor_model()
omega_plot <- omega_plot_res$plot


source("src/mass_mediation_src.R")
df <- biofinder_cog_comp |> filter(fmri_bl
) |> 
  drop_na(age, pathology_ad, sex) |> 
  mutate(has_diag = factor(diagnosis %in% c("MCI", "AD")),
         memory = resid(lm(memory ~ global_cog)),
         executive = resid(lm(executive ~ global_cog))
  ) |> 
  #left_join(cog_edit) |> 
  inner_join(fc_measures_bf$affinity|> as_tibble(rownames = NA) |> 
               rownames_to_column("image_file") |> 
               janitor::clean_names()) 

parcs <- df |> select(starts_with("x7")) |> colnames()

pw_mediation_bar_res <- plot_mediation_barplot(subject_data = df |> 
                                                 mutate(age_ = age,
                                                        pathology_ad_ = pathology_ad), 
                                               treatments = c("age", "age_", "pathology_ad", "pathology_ad_"), 
                                               outcomes = c("executive", "memory", "executive", "memory"),
                                               covariates = list(c("pathology_ad", "sex", "rsqa__MeanFD"), 
                                                                 c("pathology_ad", "sex", "rsqa__MeanFD"),
                                                                 c("age", "sex", "rsqa__MeanFD"), 
                                                                 c("age", "sex", "rsqa__MeanFD")),
                                               gradient_data = grad_df |> filter(study=="biofinder"), 
                                               gradients = c(1, 3),
                                               parcels = parcs,
                                               gradient_colors = gradient_cols,
                                               plt_title = "    Parcelwise mediation",
                                               plt_title_size = rel(0.975),
                                               terms_nice = c("Age", "Age", "ADpath", "ADpath"),
                                               point_size = 0.75,
                                               plot_spacing = 0.3,
                                               empty_row_height = -0.15,
                                               padding_med_str = 0.5,
                                               l_width = 0.05,
                                               med_lab_size = 1.5,
                                               base_size_ = 7,
                                               n_lab_size = 4.5,
                                               add_shade = TRUE,
                                               shade_alpha = 0.1,
                                               shade_size = 0.1,
                                               rasterize = TRUE,
                                               rectangle = TRUE)
pw_mediation_bar <- pw_mediation_bar_res$plot

cognition_full_new <- 
  ggdraw() +
  draw_plot(omega_plot, x = 0, y = 0.505, width = 0.40, height = 0.495) +
  draw_plot(cognitive_composites$plot, x = 0.41, y = 0.505, width = 0.59, height = 0.495) +
  draw_plot(cog_diff #+ plot_annotation(theme = theme(plot.margin = margin(10, 5, 0, 5)))
            , x = 0.0, y = 0.0, width = 0.35, height = 0.495) +
  draw_plot(pw_mediation_bar, x = 0.36, y = 0, width = 0.64, height = 0.495) +
  draw_plot_label("A", x = 0.002, y = 1, size = 8, color = "#4e4e4e") +
  draw_plot_label("B", x = 0.412, y = 1, size = 8, color = "#4e4e4e") +
  draw_plot_label("C", x = 0.002, y = 0.495, size = 8, color = "#4e4e4e") +
  draw_plot_label("D", x = 0.362, y = 0.495, size = 8, color = "#4e4e4e") +
  #draw_text("Parcelwise mediation", x = 0.365, y = 0.495, hjust = 0, vjust = 1.5, size = 20) +
  #draw_text("Covariates: \nAge/ADpath \nSex \nMotion", x = 0.3675, y = 0.485, hjust = 0, vjust = 1.5, size = 12)
  draw_text("Covariates: \nAge/ADpath \nSex \nMotion", x = 0.3675, y = 0.46, hjust = 0, vjust = 1.5, size = 5)

img_width <-  180 
p_name <- "cognition_full.pdf"
ggsave(file.path(figure_path, p_name), cognition_full_new |> pad_plot(),
       width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, 
       units = "mm", dpi = 300, 
       device = "pdf", 
       bg = "transparent")


omega_plot_res$data |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelA.csv")))

cognitive_composites$tmaps |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelB.csv")))

cog_diff_res$data |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelC.csv")))

pw_mediation_bar_res$data |> select(-c(method, affinity, sim_method, threshold)) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelD.csv")))


###########################
# Age and pathology collinearity
###########################

pal <- ggsci::pal_futurama()(3)


age_path <- biofinder_df |> 
  dplyr::filter(fmri_bl, !is.na(abnorm_ab)) |> 
  dplyr::mutate(
    Ab_pos = as.logical(abnorm_ab),
    Diagnosis = dplyr::case_when(
      diagnosis %in% c("Normal", "SCD") ~ "CU",
      TRUE ~ diagnosis
    ) |> factor(levels = c("CU", "MCI", "AD"))
  ) |> 
  ggplot(aes(age, pathology_ad, color = Diagnosis)) +
  geom_point(alpha = 0.2, size = 0.7) +
  labs(x = "Age", y = "Pathology Score") +
  scale_color_manual(
    values = c(pal[3], pal[1], pal[2]),
    labels = parse(text = c(
      "CU",
      "MCI~(A*beta^'+')",
      "AD~(A*beta^'+')"
    ))
  ) +
  theme_bw(base_size = 7) +
  theme(legend.position = "top")

p_name <- "age_path_dist.pdf"
ggsave(file.path(figure_path, p_name), age_path |> pad_plot(),
       width = 88, height = 60, units = "mm", 
       dpi = 300, device = "pdf", bg = "white")

biofinder_df |> 
  dplyr::filter(fmri_bl, !is.na(abnorm_ab)) |> 
  dplyr::mutate(
    Ab_pos = as.logical(abnorm_ab),
    Diagnosis = dplyr::case_when(
      diagnosis %in% c("Normal", "SCD") ~ "CU",
      TRUE ~ diagnosis
    ) |> factor(levels = c("CU", "MCI", "AD"))
  ) |> 
  select(age, pathology_ad, Diagnosis) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))




########################################################
# Extended data
########################################################



########################################################
# Extended data figure with connectivity strength
########################################################


strength_supp <- figure_one(subject_data = biofinder_df,
                            subject_data_replication = adni_df |> filter(fmri_bl),
                            measures_list =list(nodal_strength = fc_measures_bf$strength), 
                            measures_list_replication = list(nodal_strength = fc_measures_adni$strength),
                            gradients_df = grad_df |> filter(study == "biofinder"),
                            gradients_df_replication = grad_df |> filter(study == "adni"),
                            selected_gradients = c(1, 3),
                            shade = TRUE,
                            shade_a = 0.1,
                            shade_s = 0.1,
                            tag_size = 10,
                            p_size = 0.05,
                            raster = TRUE,
                            ggrastr_dpi = 300,
                            empt_row_height = -0.15,
                            b_size = 7,
                            draw_size = 7,
                            plot_title_size = 0.95,
                            axes_title_size = 0.9,
                            r2_sizing2 = rel(0.7),
                            r2_sizing1 = rel(0.9),
                            ax_txt_size = 5,
                            boxed = TRUE,
                            split = FALSE,
                            brain_plot_names_f1 = c(NA, "AD Pathology"))

p_name <- "supplementary_strength.pdf"
ggsave(file.path(figure_path, p_name), strength_supp[[1]] |> pad_plot(), 
       width = img_width*scaling_factor, 
       height = img_width*0.8*scaling_factor,
       units = "mm", dpi = 300, device = "pdf")


source_data <- between_supp$tmaps$fig1$bf |> mutate(panel = "A") |> 
  bind_rows(between_supp$tmaps$fig1$adni |> mutate(panel = "B")) |> 
  bind_rows(between_supp$tmaps$fig2$health |> mutate(panel = "C")) |> 
  bind_rows(between_supp$tmaps$fig2$clin |> mutate(panel = "D")) |> 
  select(-n) |> 
  left_join(grad_df |> filter(study %in% c("adni", "biofinder")) |> 
              select(region, study, starts_with("gradient")) |>
              pivot_wider(values_from = starts_with("gradient"), names_from = "study", names_sep = "_"))

write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))



################
# Main overlaid
################

source("src/util_vis.R")
img_width = 120 
scaling_factor <-  1
net_main_res <- overlaid_main_results(subject_data = biofinder_df |> filter(fmri_bl), 
                                      fc_matrix = fc_measures_bf$affinity,
                                      parc_width = 0.05,
                                      net_width = 0.3,
                                      leg_size = 0.2,
                                      leg_text_size = 0.9) 


p_name <- "main_res_with_overlaid_nets.pdf"
ggsave(file.path(figure_path, p_name), net_main_res$plot |> pad_plot(),
       width = img_width*scaling_factor, height = img_width*0.7*scaling_factor, units = "mm", dpi = 300, device = "pdf")

net_main_res$source_data |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))


###########################
# Gradient coverage
###########################

source("src/plot_gradient_coverage.R")

gradient_cov <- plot_gradient_coverage(grad_df = grad_df |> filter(study == "biofinder"),
                                       nq_filter_terms = NULL,#test$nq_filter_terms,
                                       rasterize = TRUE,
                                       wordcloud_max_size = 4.9,
                                       wordcloud_grid_margin = 0.95,
                                       network_legend_text_scale = 1,
                                       network_legend_key_lines = 1,
                                       base_size = 7)

p_name <- "gradient_coverage.pdf"
ggsave(file.path(figure_path, p_name), 
       gradient_cov$plot |> pad_plot(), 
       width=180, 
       height=180, 
       units = "mm",
       bg="white")

gradient_cov$data$g1_wordcloud |> mutate(panel = "SA wordcloud") |> 
  bind_rows(gradient_cov$data$g3_wordcloud |> mutate(panel = "RE wordcloud")) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_wordclouds.csv")))

gradient_cov$data$network_separation |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_boxplots.csv")))



###################################
# Broken down by pathology and APOE
###################################


tau_ab_cu_adj <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, !is.na(pathology_ad)) |> 
                                               mutate(has_diagnosis = diagnosis %in% c("MCI", "AD")), 
                                             gradient_data = grad_df |> filter(study== "biofinder"), 
                                             gradients = c(1, 3),
                                             gradient_colors = gradient_cols,
                                             list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                             vect = TRUE,
                                             empty_row_height = -0.15,
                                             brain_title_size = 6,
                                             axis_text_size = 5,
                                             axis_title_size = 6,
                                             scatter_title_vjust = 0,
                                             base_size_ = 7,
                                             rasterize = TRUE,
                                             ggrastr_dpi = 300,
                                             add_shade = TRUE, 
                                             shade_alpha = 0.1,
                                             shade_size = 0.1,
                                             r_spin_size = 0.6,
                                             point_size = 0.05,
                                             point_alpha = 0.15,
                                             mod_formula = formula(paste0("~ age + tau_pathology + abnorm_ab + has_diagnosis + sex + rsqa__MeanFD")),
                                             covariates = c("sex", "rsqa__MeanFD"),
                                             plt_title = "Full sample",
                                             plt_subtitle = TRUE,
                                             subtit_lookup = c(abnorm_ab = "Ab_pos", rsqa__MeanFD = "motion", has_diagnosis = "MCI_or_AD"),
                                             brain_names = c("Age", "Tau Pathology", "A*beta*\"+\"", "MCI/AD"),
                                             brain_names_parse = c(FALSE, FALSE, TRUE, FALSE),
                                             rectangle = TRUE,
                                             plot_net_legend = FALSE,
                                             net_legend_x = 0.1,
                                             net_legend_y = 0.01)


tau_ab_cu <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, !is.na(pathology_ad), diagnosis %in% c("Normal", "SCD")) |> 
                                           mutate(diag = diagnosis %in% c("Normal", "SCD")), 
                                         gradient_data = grad_df |> filter(study== "biofinder"), 
                                         gradients = c(1, 3),
                                         gradient_colors = gradient_cols,
                                         list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                         vect = TRUE,
                                         empty_row_height = -0.15,
                                         brain_title_size = 6,
                                         axis_text_size = 5,
                                         axis_title_size = 6,
                                         scatter_title_vjust = 0,
                                         base_size_ = 7,
                                         rasterize = TRUE,
                                         ggrastr_dpi = 300,
                                         add_shade = TRUE, 
                                         shade_alpha = 0.1,
                                         shade_size = 0.1,
                                         r_spin_size = 0.60,
                                         point_size = 0.05,
                                         point_alpha = 0.15,
                                         mod_formula = formula(paste0("~ age + tau_pathology + abnorm_ab + sex + rsqa__MeanFD")),
                                         covariates = c("sex", "rsqa__MeanFD"),
                                         plt_title = "Cognitively Unimpaired",
                                         plt_title_size = rel(0.9),
                                         plt_subtitle_size = rel(0.8),
                                         plt_subtitle = TRUE,
                                         subtit_lookup = c(abnorm_ab = "Ab_pos", rsqa__MeanFD = "motion"),
                                         brain_names = c("Age", "Tau Pathology", "A*beta*\"+\""),
                                         brain_names_parse = c(FALSE, FALSE, TRUE),
                                         rectangle = TRUE,
                                         plot_net_legend = FALSE,
                                         net_legend_x = 0.2,
                                         net_legend_y = 0.01)





##################################
# APOE
#################################


apoe <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl#, abnorm_ab == 0, pathology_ad<0.15
                                                            ) |> 
                                      mutate(has_diag = diagnosis %in% c("MCI", "AD")), 
                                    gradient_data = grad_df |> filter(study== "biofinder"), 
                                    gradients = c(1, 3),
                                    gradient_colors = gradient_cols,
                                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                    empty_row_height = -0.15,
                                    brain_title_size = 6,
                                    axis_text_size = 5,
                                    axis_title_size = 6,
                                    scatter_title_vjust = 0,
                                    base_size_ = 7,
                                    vect = TRUE,
                                    rasterize = TRUE,
                                    ggrastr_dpi = 300,
                                    add_shade = TRUE, 
                                    shade_alpha = 0.1,
                                    shade_size = 0.1,
                                    r_spin_size = 0.75,
                                    point_size = 0.05,
                                    point_alpha = 0.15,
                                    plot_net_legend = TRUE,
                                    net_legend_x = 0.2,
                                    net_legend_y = 0.01,
                                    mod_formula = formula(paste0("~ age + pathology_ad + apoe4 + sex + rsqa__MeanFD")),
                                    covariates = c("sex", "rsqa__MeanFD"),
                                    plt_title = "Full sample",
                                    title_lmargin = 10,
                                    brain_names = c("Age", "AD Pathology", "epsilon*4~carrier"),
                                    brain_names_parse = c(FALSE, FALSE, TRUE),
                                    plt_subtitle = TRUE,
                                    subtit_lookup = c(apoe4 = "e4 carrier", rsqa__MeanFD = "motion", pathology_ad = "AD pathology"),
                                    rectangle = TRUE)


apoe_abneg <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, 
                                                                  diagnosis %in% c("Normal", "SCD")), 
                                          gradient_data = grad_df |> filter(study== "biofinder"), 
                                          gradients = c(1, 3),
                                          gradient_colors = gradient_cols,
                                          list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                          empty_row_height = -0.15,
                                          brain_title_size = 6,
                                          axis_text_size = 5,
                                          axis_title_size = 6,
                                          scatter_title_vjust = 0,
                                          base_size_ = 7,
                                          vect = TRUE,
                                          rasterize = TRUE,
                                          ggrastr_dpi = 300,
                                          add_shade = TRUE, 
                                          shade_alpha = 0.1,
                                          shade_size = 0.1,
                                          r_spin_size = 0.75,
                                          point_size = 0.05,
                                          point_alpha = 0.15,
                                          plot_net_legend = FALSE,
                                          net_legend_x = 0.2,
                                          net_legend_y = 0.01,
                                          mod_formula = formula(paste0("~ age + tau_pathology + apoe4 + sex + rsqa__MeanFD")),
                                          covariates = c("sex", "rsqa__MeanFD"),
                                          plt_title = '"Cognitively Unimpaired (A"*beta*"-)"',
                                          plt_title_parse = TRUE,
                                          brain_names = c("Age", "Tau Pathology", "epsilon*4~carrier"),
                                          brain_names_parse = c(FALSE, FALSE, TRUE),
                                          plt_subtitle = TRUE,
                                          subtit_lookup = c(apoe4 = "e4 carrier", rsqa__MeanFD = "motion", pathology_ad = "AD pathology"),
                                          rectangle = TRUE)

l_marg = 10
apoe_p_and_tauab <- ggdraw() +
  draw_plot(tau_ab_cu_adj$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=l_marg)),
                                                               plot.subtitle = element_text(margin = margin(l=l_marg)))),
            x = 0, y = 0.505, width = 0.555, height = 0.495) +
  draw_plot(tau_ab_cu$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=l_marg)),
                                                           plot.subtitle = element_text(margin = margin(l=l_marg))))
            , x = 0.565, y = 0.505, width = 1 - 0.565, height = 0.495) +
  draw_plot(apoe$plot
            , x = 0, y = 0, width = 0.495, height = 0.495) +
  draw_plot(apoe_abneg$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=l_marg)),
                                                           plot.subtitle = element_text(margin = margin(l=l_marg))))
            , x = 0.505, y = 0, width = 0.495, height = 0.495) +
  draw_plot_label("A", size = 9) +
  draw_plot_label("B", x = 0.570, size = 9) +
  draw_plot_label("C", x = 0.0, y = 0.495, size = 9) +
  draw_plot_label("D", x = 0.505, y = 0.495, size = 9) 


scaling_factor <- 1
img_width = 180
p_name <- "apoe_plot_v1_w_tau_ab.pdf"
ggsave(file.path(figure_path, p_name), apoe_p_and_tauab |> pad_plot(),
       width = img_width*scaling_factor, 
       height = img_width*0.8*scaling_factor,
       units = "mm", dpi = 300, device = "pdf", bg = "white")


tau_ab_cu_adj$tmaps |> mutate(panel = "A") |> 
  bind_rows(tau_ab_cu$tmaps |> mutate(panel = "B")) |> 
  bind_rows(apoe$tmaps |> mutate(panel = "C")) |> 
  bind_rows(apoe_abneg$tmaps |> mutate(panel = "D")) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))



###########################
# Nodal Tau
###########################

tau <- readRDS("data/bf_src_data/regional_tau/tau_pet_schaefer_1000_r.rds")

nodal_tau <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, !is.na(pathology_ad)), 
                                         t_mat = tau,
                                         gradient_data = grad_df |> filter(study== "biofinder"), 
                                         gradients = c(1, 3),
                                         gradient_colors = gradient_cols,
                                         list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                         empty_row_height = -0.2,
                                         brain_title_size = 6,
                                         axis_text_size = 5,
                                         axis_title_size = 6,
                                         scatter_title_vjust = 0,
                                         base_size_ = 7,
                                         vect = TRUE,
                                         rasterize = TRUE,
                                         ggrastr_dpi = 300,
                                         add_shade = TRUE, 
                                         shade_alpha = 0.1,
                                         shade_size = 0.1,
                                         r_spin_size = 0.75,
                                         point_size = 0.05,
                                         point_alpha = 0.3,
                                         plot_net_legend = TRUE,
                                         net_legend_x = 0.2,
                                         net_legend_y = 0.01,
                                         mod_formula = formula(paste0("~ age + sex + rsqa__MeanFD")),
                                         covariates = c("sex", "rsqa__MeanFD"),
                                         plt_title = "Parcelwise tau on parcelwise FCS regression",
                                         plt_subtitle = TRUE,
                                         rectangle = TRUE)


scaling_factor <- 1
img_width = 88
p_name <- "nodal_tau.pdf"
ggsave(file.path(figure_path, p_name), nodal_tau$plot |> pad_plot(),
       width = img_width*scaling_factor, height = img_width*0.9*scaling_factor, units = "mm", dpi = 300, device = "pdf", bg = "white")


source_data <- nodal_tau$tmaps |> 
  select(-n, -model_formula) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient")))
write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))



###########################
# Atrophy + vascular
###########################

new_data <- read_delim("stuff_for_revisions/jonathan__20250926_104238.csv")

x <- new_data |> 
  mutate(mri_date = as.Date(str_extract(csv_ct_composites__index, "(?<=__).*"), format = "%Y%m%d")) |> 
  filter(!is.na(mri_date)) |> 
  mutate(
    mri_bl = mri_date == min(mri_date),
    has_longitudinal = n()>1, .by = "sid"
  ) |> 
  select(sid, mri_bl, visit = Visit, visit_date, mri_date,
         contains("microbleeds"),
         contains("lacunes"),
         contains("vol"),
         samseg_wmhs_WMH_total_mm3,
         (contains("ct") & contains("lr"))) |> 
  #select(!contains("adsign")) |> 
  mutate(ct_mean_lr = rowMeans(across(contains("ct") & contains("lr") & contains("surfwt")), na.rm = TRUE),
         microbleeds = factor(as.logical(Cerebral_microbleeds_dich)),
         hippo_vol = -rowMeans(across(contains("Hippocampus") & contains("vols"))), 
         lacunes = factor(as.logical(Lacunes_dich)),
         wmhs_tot = samseg_wmhs_WMH_total_mm3,
         ct_tot = (ct_lateraltemporal_surfwt_lr + ct_medialtemporal_surfwt_lr + 
                     ct_lateralparietal_surfwt_lr + ct_medialparietal_surfwt_lr 
                   + ct_frontal_surfwt_lr + ct_occipital_surfwt_lr)/6) |> 
  filter(mri_bl)


biofinder_new <- biofinder_df |> filter(fmri_bl) |> 
  inner_join(x, by = c("sid"))

cerebro_vascular <- plot_gradient_relationships(biofinder_new %>% filter(fmri_bl) |>
                                                  mutate(
                                                    apoe4 = factor(apoe4),
                                                    has_diag = diagnosis %in% c("MCI", "AD"),
                                                    CT_adsig = -ct_adsign_lr,
                                                    ct_tot = -ct_tot,
                                                    a_syn = factor(a_syn)),
                                                gradient_data = grad_df |> filter(study== "biofinder"),
                                                gradients = c(1, 3),
                                                gradient_colors = gradient_cols,
                                                list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                empty_row_height = -0.15,
                                                brain_title_size = 6,
                                                axis_text_size = 5,
                                                axis_title_size = 6,
                                                scatter_title_vjust = 0,
                                                base_size_ = 7,
                                                rasterize = TRUE,
                                                ggrastr_dpi = 300,
                                                add_shade = TRUE, 
                                                shade_alpha = 0.1,
                                                shade_size = 0.1,
                                                r_spin_size = 0.75,
                                                point_size = 0.05,
                                                point_alpha = 0.2,
                                                plot_net_legend = TRUE,
                                                net_legend_x = 0.08,
                                                net_legend_y = 0.01,
                                                brain_title_vjust = -3,
                                                brain_subtitle_vjust = -4,
                                                l_width = 0.025,
                                                vect = TRUE,
                                                mod_formula = formula(paste0(" ~ age + pathology_ad + ct_tot + wmhs_tot + lacunes + microbleeds + a_syn + sex + rsqa__MeanFD")),
                                                covariates = c("sex", "rsqa__MeanFD"),
                                                plt_title = "Effects of age and AD pathology on FCS are not explained by other age-related factors",
                                                brain_names = c("Age", "AD Pathology", "-CT total", "WMHS tot", "Lacunes+", "Microbleeds+", 'alpha*"Syn+"'),
                                                brain_names_parse = c(rep(FALSE, 6), TRUE),
                                                #brain_title_vjust = -5,
                                                #brain_subtitle_vjust = -7.5,
                                                plt_subtitle = TRUE,
                                                subtit_lookup = c(pathology_ad = "AD pathology", ct_tot = "-ct_tot", a_syn = "a_syn_pos", rsqa__MeanFD = "motion"),
                                                group_n_title = TRUE,
                                                rectangle = TRUE)


img_width <-  180 
p_name <- "atrophy_and_cerebrovasc.pdf"
ggsave(file.path(figure_path, p_name), cerebro_vascular$plot |> pad_plot(),
       width = img_width*scaling_factor, height = img_width*0.45*scaling_factor, 
       units = "mm", dpi = 300, #device = "png",
       bg = "white")

cerebro_vascular$tmaps |>
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))


##################################
# Causality plot
##################################


orig_analysis <- nodal_regression_fits(
  biofinder_df %>% filter(fmri_bl), 
  fc_measures_bf$affinity,
  roi_names = rois,
  vectorised = TRUE,
  model_formula = formula(paste0(" ~ age + pathology_ad + sex  + rsqa__MeanFD"))
)

orig_analysis <- get_nodal_ests(orig_analysis)

g1 <- grad_df |> filter(study == "biofinder") |> pull(gradient1)
g3 <- grad_df |> filter(study == "biofinder") |> pull(gradient3)
subj_corr_g1 <- apply(fc_measures_bf$affinity[biofinder_df %>% filter(fmri_bl) |> pull(image_file), ], 1 , function(x) cor(x, g1))
subj_corr_g3 <- apply(fc_measures_bf$affinity[biofinder_df %>% filter(fmri_bl) |> pull(image_file), ], 1 , function(x) cor(x, g3))

tibble(`SA corr` = subj_corr_g1, `RE corr` = subj_corr_g3) |> 
  pivot_longer(everything(), names_to = "Gradient", values_to = "Correlation") |> 
  ggplot(aes(x = Correlation)) +
  geom_histogram(color = "black", linewidth = 0.1) +
  theme_bw(base_size = 7) +
  labs(title = "Individual correlation with gradients" , y = "") +
  facet_wrap(~Gradient) +
  theme(plot.title = element_text(size = rel(0.85)),
        plot.background = element_blank()) -> subject_corr

sexed <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                     gradient_data = grad_df |> filter(study== "biofinder"), 
                                     gradients = c(1, 3),
                                     gradient_colors = gradient_cols,
                                     list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                     brain_title_size = 6,
                                     axis_text_size = 5,
                                     axis_title_size = 6,
                                     scatter_title_vjust = 0,
                                     rasterize = TRUE,
                                     ggrastr_dpi = 300,
                                     add_shade = TRUE, 
                                     shade_alpha = 0.1,
                                     shade_size = 0.1,
                                     point_size = 0.05,
                                     point_alpha = 0.2,
                                     l_width = 0.05,
                                     empty_row_height = -0.15,
                                     base_size_ = 7,
                                     brain_title_vjust = -3,
                                     r_spin_size = 0.775,
                                     vect = TRUE,
                                     plt_subtitle = TRUE,
                                     group_n_title = FALSE,
                                     rectangle = TRUE,
                                     plot_net_legend = TRUE,
                                     net_legend_x = 0.15,
                                     net_legend_y = 0.01,
                                     mod_formula = formula(" ~ age + pathology_ad + education + sex + rsqa__MeanFD"),
                                     covariates = c("rsqa__MeanFD"),
                                     plt_title = "Full Sample",
                                     title_lmargin = 10,
                                     brain_names = c("Age", "AD Pathology", "Education", "Sex")
)


mean_aff <- fc_measures_bf$affinity[biofinder_df %>% filter(fmri_bl) |> pull(image_file), ] |> colMeans()
t_path <- orig_analysis |> filter(term == "pathology_ad") |> pull(statistic)
t_age <- orig_analysis |> filter(term == "age") |> pull(statistic)

scale_fac <- 7/20

cor_path <- ppcor::pcor(data.frame(SA_axis = g1, tmap_path = t_path, Affinity = mean_aff))$estimate |> 
  as.data.frame() |> 
  rownames_to_column("Var1") |> 
  pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") |> 
  mutate(
    Var1 = str_replace_all(Var1, "_", " "),
    Var2 = str_replace_all(Var2, "_", " "),
    Var1 = factor(Var1, levels = unique(Var1)),
    Var2 = factor(Var2, levels = rev(unique(Var1)))  # reverse order
  ) |> 
  ggplot(aes(x = Var1, y = Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = muted("blue"), high = muted("red"), midpoint = 0,   limits = c(-1, 1)) +
  scale_x_discrete(position = "top") +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 6*scale_fac, color = "darkgray") +
  coord_fixed() +
  labs(title = "Partial correlation pathology",
       x = NULL, y = NULL, fill = "Correlation") +
  theme_minimal(base_size = 7) +
  theme(
    legend.position = "",
    axis.text = element_text(size = rel(0.65)),
    axis.text.x = element_text(angle = 45, vjust = -1, hjust = 0),
    panel.grid = element_blank(),
    plot.title = element_text(size = rel(0.85)))

cor_age <- ppcor::pcor(data.frame(RE_axis = g3, tmap_age = t_age, Affinity = mean_aff))$estimate |> 
  as.data.frame() |> 
  rownames_to_column("Var1") |> 
  pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") |> 
  mutate(
    Var1 = str_replace_all(Var1, "_", " "),
    Var2 = str_replace_all(Var2, "_", " "),
    Var1 = factor(Var1, levels = unique(Var1)),
    Var2 = factor(Var2, levels = rev(unique(Var1)))  # reverse order
  ) |> 
  ggplot(aes(x = Var1, y = Var2, fill = correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = muted("blue"), high = muted("red"), midpoint = 0,   limits = c(-1, 1)) +
  geom_text(aes(label = sprintf("%.2f", correlation)), size = 6*scale_fac, color = "darkgray") +
  scale_x_discrete(position = "top") +
  coord_fixed() +
  labs(title = "Partial correlation age",
       x = NULL, y = NULL, fill = "Correlation") +
  theme_minimal(base_size = 7) +
  theme(
    legend.position = "",
    axis.text = element_text(size = rel(0.65)),
    axis.text.x = element_text(angle = 45, vjust = -1, hjust = 0),
    panel.grid = element_blank(),
    plot.title = element_text(size = rel(0.85))
  )


part_cor <- (cor_path + cor_age) / subject_corr
caus_plot <-
  ggdraw() +
  draw_plot(sexed$plot, x = 0, y = 0, width = 0.6, height = 1) +
  draw_plot_label("A", size = 7) +
  draw_plot(part_cor, x = 0.575, y = 0, width = 0.425, height = 1) +
  draw_plot_label("B", x = 0.6, size = 7) +
  draw_plot_label("C", x = 0.6, y = 0.55, size = 7) 

img_width <- 180
scaling_factor <- 1
p_name <- "caus_plot.pdf"
ggsave(file.path(figure_path, p_name), caus_plot |> pad_plot(),
       width = img_width*scaling_factor, height = img_width*0.45*scaling_factor, 
       units = "mm", dpi = 300, device = "pdf", bg = "white")

sexed$tmaps |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelA.csv")))

ppcor::pcor(data.frame(SA_axis = g1, tmap_path = t_path, Affinity = mean_aff))$estimate |> 
  as.data.frame() |> 
  rownames_to_column("Var1") |> 
  pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") |> 
  mutate(panel = "Partial correlation pathology") |> 
  bind_rows(ppcor::pcor(data.frame(RE_axis = g3, tmap_age = t_age, Affinity = mean_aff))$estimate |> 
              as.data.frame() |> 
              rownames_to_column("Var1") |> 
              pivot_longer(-Var1, names_to = "Var2", values_to = "correlation") |> 
              mutate(panel = "Partial correlation age")) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelB.csv")))

tibble(`SA corr` = subj_corr_g1, `RE corr` = subj_corr_g3) |> 
  pivot_longer(everything(), names_to = "Gradient", values_to = "Correlation") |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelC.csv")))




 ######  ##     ## ########  ########  ##       ######## ##     ## ######## ##    ## ########    ###    ########  ##    ## 
##    ## ##     ## ##     ## ##     ## ##       ##       ###   ### ##       ###   ##    ##      ## ##   ##     ##  ##  ##  
##       ##     ## ##     ## ##     ## ##       ##       #### #### ##       ####  ##    ##     ##   ##  ##     ##   ####   
 ######  ##     ## ########  ########  ##       ######   ## ### ## ######   ## ## ##    ##    ##     ## ########     ##    
      ## ##     ## ##        ##        ##       ##       ##     ## ##       ##  ####    ##    ######### ##   ##      ##    
##    ## ##     ## ##        ##        ##       ##       ##     ## ##       ##   ###    ##    ##     ## ##    ##     ##    
 ######   #######  ##        ##        ######## ######## ##     ## ######## ##    ##    ##    ##     ## ##     ##    ##    


########################################################
# Supplementary figure with diagnoses
########################################################
source("src/util_vis.R")

supplementary_figure_path <- "paper/suppfig_original"

bf_dx <-  plot_gradient_relationships(biofinder_df %>% 
                                        filter(fmri_bl) %>%
                                        mutate(motion = rsqa__MeanFD
                                        ) %>% 
                                        mutate(diagnosis = ifelse(diagnosis %in% c("Normal", "SCD"), "CN/SCD", diagnosis)) %>% 
                                        mutate(diagnosis=factor(diagnosis, levels = c("CN/SCD", "MCI", "AD"))), 
                                      gradient_data = grad_df %>% filter(study=="biofinder"), 
                                      gradients = c(1, 2, 3),
                                      empty_row_height = -0.1,
                                      brain_title_size = rel(1.05),
                                      axis_text_size = 5,
                                      axis_title_size = 6,
                                      scatter_title_vjust = 0,
                                      base_size_ = 7,
                                      vect = TRUE,
                                      rasterize = TRUE,
                                      ggrastr_dpi = 300,
                                      add_shade = TRUE, 
                                      shade_alpha = 0.1,
                                      shade_size = 0.1,
                                      r_spin_size = 0.9,
                                      point_size = 0.05,
                                      point_alpha = 0.3,
                                      rectangle = TRUE,
                                      plot_net_legend = TRUE,
                                      net_legend_x = 0.12,
                                      net_legend_y = 0.009,
                                      gradient_colors = gradient_cols,
                                      list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                      mod_formula = formula(paste0(" ~ age + diagnosis + sex + rsqa__MeanFD")),
                                      covariates = c("sex", "rsqa__MeanFD"),
                                      layout_construction = "horizontal",
                                      plt_title = "BioFINDER",
                                      plt_subtitle = TRUE,
                                      brain_names = c("Age", "MCI vs. CN", "AD vs. CN")
                                      )


img_width <- 140
p_name <- "supplementary_diagnosis.pdf"
ggsave(file.path(supplementary_figure_path, p_name), bf_dx$plot |> pad_plot(), width = img_width, 
       bg = "white", height = img_width*0.9, units = "mm", dpi = 300)

source_data <- bf_dx$tmaps |> 
  select(-n, -model_formula) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient")))
write_csv(source_data, file.path(source_figure_path, paste0("supplementary_diagnosis_estimates", ".csv")))

########################################################
# Supplementary figure with gradient 2
########################################################
source("src/util_vis.R")
fig_one <- figure_one(subject_data = biofinder_df,
                      subject_data_replication = adni_df |> filter(fmri_bl),
                      measures_list =list(nodal_affinity = fc_measures_bf$affinity), 
                      measures_list_replication = list(nodal_affinity = fc_measures_adni$affinity),
                      gradients_df = grad_df |> filter(study == "biofinder"),
                      gradients_df_replication = grad_df |> filter(study == "adni"),
                      selected_gradients = c(2), 
                      shade = TRUE,
                      shade_a = 0.1,
                      shade_s = 0.1,
                      tag_size = 10,
                      p_size = 0.2,
                      raster = TRUE,
                      ggrastr_dpi = 300,
                      empt_row_height = -0.15,
                      b_size = 7,
                      draw_size = 7,
                      plot_title_size = 0.9,
                      axes_title_size = 0.9,
                      ax_txt_size = rel(0.8),
                      r2_sizing1 = 0.85,
                      brain_title_size = rel(1.05),
                      boxed = TRUE,
                      split = TRUE,
                      fig = 1,
                      brain_plot_names_f1 = c(NA, "AD Pathology"))

img_width = 180 
scaling_factor <-  1

p_name <- "supplementary_gradient2.pdf"
ggsave(file.path(supplementary_figure_path, p_name), fig_one[[1]] |> pad_plot(),
       width = img_width*scaling_factor, height = img_width*0.35*scaling_factor, 
       units = "mm", dpi = 300,  bg = "white")

source_data <- fig_one$tmaps$bf |> mutate(study = "biofinder") |> 
  bind_rows(fig_one$tmaps$adni |> mutate(study = "adni")) |> 
  select(-n, -model_formula) |> 
  left_join(grad_df |> filter(study %in% c("adni", "biofinder")) |> 
              select(region, study, starts_with("gradient")))
write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))



########################################################
# Supplementary figure within network affinity
########################################################


within_supp <- figure_one(subject_data = biofinder_df,
                          subject_data_replication = adni_df |> filter(fmri_bl),
                          measures_list =list(within_network_affinity = fc_measures_bf$affinity_within), 
                          measures_list_replication = list(within_network_affinity = fc_measures_adni$affinity_within),
                          gradients_df = grad_df |> filter(study == "biofinder"),
                          gradients_df_replication = grad_df |> filter(study == "adni"),
                          selected_gradients = c(1, 3),
                          shade = TRUE,
                          shade_a = 0.1,
                          shade_s = 0.1,
                          tag_size = 10,
                          p_size = 0.05,
                          raster = TRUE,
                          ggrastr_dpi = 300,
                          empt_row_height = -0.15,
                          b_size = 7,
                          draw_size = 7,
                          plot_title_size = 0.95,
                          axes_title_size = 0.9,
                          r2_sizing2 = rel(0.7),
                          r2_sizing1 = rel(0.9),
                          ax_txt_size = 5,
                          boxed = TRUE,
                          split = FALSE,
                          brain_plot_names_f1 = c(NA, "AD Pathology"))
img_width = 180 
scaling_factor <-  1
p_name <- "within_supplementary.pdf"
ggsave(file.path(supplementary_figure_path, p_name), within_supp[[1]] |> pad_plot(), 
       width = img_width*scaling_factor, 
       height = img_width*0.8*scaling_factor,
       units = "mm", dpi = 300, device = "pdf")


source_data <- within_supp$tmaps$fig1$bf |> mutate(panel = "A") |> 
  bind_rows(within_supp$tmaps$fig1$adni |> mutate(panel = "B")) |> 
  bind_rows(within_supp$tmaps$fig2$health |> mutate(panel = "C")) |> 
  bind_rows(within_supp$tmaps$fig2$clin |> mutate(panel = "D")) |> 
  select(-n) |> 
  left_join(grad_df |> filter(study %in% c("adni", "biofinder")) |> 
              select(region, study, starts_with("gradient")) |>
              pivot_wider(values_from = starts_with("gradient"), names_from = "study", names_sep = "_"))
  
write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))



########################################################
# Supplementary figure between network affinity
########################################################

between_supp <- figure_one(subject_data = biofinder_df,
                           subject_data_replication = adni_df |> filter(fmri_bl),
                           measures_list =list(between_network_affinity = fc_measures_bf$affinity_between), 
                           measures_list_replication = list(between_network_affinity = fc_measures_adni$affinity_between),
                           gradients_df = grad_df |> filter(study == "biofinder"),
                           gradients_df_replication = grad_df |> filter(study == "adni"),
                           selected_gradients = c(1, 3),
                           shade = TRUE,
                           shade_a = 0.1,
                           shade_s = 0.1,
                           tag_size = 10,
                           p_size = 0.05,
                           raster = TRUE,
                           ggrastr_dpi = 300,
                           empt_row_height = -0.15,
                           b_size = 7,
                           draw_size = 7,
                           plot_title_size = 0.95,
                           axes_title_size = 0.9,
                           r2_sizing2 = rel(0.7),
                           r2_sizing1 = rel(0.9),
                           ax_txt_size = 5,
                           boxed = TRUE,
                           split = FALSE,
                           brain_plot_names_f1 = c(NA, "AD Pathology"))

img_width = 180
p_name <- "between_supplementary.pdf"
ggsave(file.path(supplementary_figure_path, p_name), between_supp[[1]] |> pad_plot(), 
       width = img_width*scaling_factor, 
       height = img_width*0.8*scaling_factor,
       units = "mm", dpi = 300, device = "pdf")


source_data <- between_supp$tmaps$fig1$bf |> mutate(panel = "A") |> 
  bind_rows(between_supp$tmaps$fig1$adni |> mutate(panel = "B")) |> 
  bind_rows(between_supp$tmaps$fig2$health |> mutate(panel = "C")) |> 
  bind_rows(between_supp$tmaps$fig2$clin |> mutate(panel = "D")) |> 
  select(-n) |> 
  left_join(grad_df |> filter(study %in% c("adni", "biofinder")) |> 
              select(region, study, starts_with("gradient")) |>
              pivot_wider(values_from = starts_with("gradient"), names_from = "study", names_sep = "_"))

write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))


########################################################
# Supplementary figure with corr + no threshold affinity
########################################################

aff_no_thresh <- figure_one(subject_data = biofinder_df,
                            subject_data_replication = adni_df |> filter(fmri_bl),
                            measures_list =list(nodal_affinity = bf_aff_nothres$affinity), 
                            measures_list_replication = list(nodal_affinity = adni_affinity_no_thresh$affinity),
                            gradients_df = grad_df |> filter(study == "biofinder"),
                            gradients_df_replication = grad_df |> filter(study == "adni"),
                            selected_gradients = c(1, 3),
                            shade = TRUE,
                            shade_a = 0.1,
                            shade_s = 0.1,
                            tag_size = 10,
                            p_size = 0.05,
                            raster = TRUE,
                            ggrastr_dpi = 300,
                            empt_row_height = -0.15,
                            b_size = 7,
                            draw_size = 7,
                            plot_title_size = 0.95,
                            axes_title_size = 0.9,
                            r2_sizing2 = rel(0.7),
                            r2_sizing1 = rel(0.9),
                            ax_txt_size = 5,
                            boxed = TRUE,
                            split = FALSE,
                            brain_plot_names_f1 = c(NA, "AD Pathology"))

img_width = 180
p_name <- "supplementary_affinity_no_thresh_correlation.pdf"
ggsave(file.path(supplementary_figure_path, p_name), aff_no_thresh[[1]] |> pad_plot(), 
       width = img_width*scaling_factor, 
       height = img_width*0.8*scaling_factor,
       units = "mm", dpi = 300, device = "pdf")


source_data <- aff_no_thresh$tmaps$fig1$bf |> mutate(panel = "A") |> 
  bind_rows(aff_no_thresh$tmaps$fig1$adni |> mutate(panel = "B")) |> 
  bind_rows(aff_no_thresh$tmaps$fig2$health |> mutate(panel = "C")) |> 
  bind_rows(aff_no_thresh$tmaps$fig2$clin |> mutate(panel = "D")) |> 
  select(-n) |> 
  left_join(grad_df |> filter(study %in% c("adni", "biofinder")) |> 
              select(region, study, starts_with("gradient")) |>
              pivot_wider(values_from = starts_with("gradient"), names_from = "study", names_sep = "_"))

write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))






######################################################################
# Supplementary figure for analysis in clinical group with interaction
######################################################################

clinical_cog_int <-  plot_gradient_relationships(biofinder_df |> filter(fmri_bl, diagnosis=="MCI" | diagnosis=="AD", !is.na(mPACC_v1)) |> 
                                                   mutate(`-mPACC_v1` = -mPACC_v1), 
                                                 gradient_data = grad_df |> filter(study=="biofinder"), 
                                                 gradients = c(1, 3),
                                                 vect = TRUE,
                                                 gradient_colors = gradient_cols,
                                                 list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                 mod_formula = formula(paste0(" ~ age + scale(pathology_ad)* scale(-mPACC_v1)  + sex + rsqa__MeanFD")),
                                                 covariates = c("sex", "rsqa__MeanFD"),
                                                 empty_row_height = -0.15,
                                                 brain_title_size = 6,
                                                 axis_text_size = 5,
                                                 axis_title_size = 6,
                                                 scatter_title_vjust = 0,
                                                 base_size_ = 7,
                                                 rasterize = TRUE,
                                                 ggrastr_dpi = 300,
                                                 add_shade = TRUE, 
                                                 shade_alpha = 0.1,
                                                 shade_size = 0.1,
                                                 r_spin_size = 0.75,
                                                 point_size = 0.05,
                                                 point_alpha = 0.3,
                                                 rectangle = TRUE,
                                                 plt_title = c("Diagnosed MCI/AD (Ab+)"),
                                                 plt_subtitle = TRUE,
                                                 plot_net_legend = TRUE,
                                                 brain_names = c("Age", "AD Pathology", "-mPACC", "AD Pathology × -mPACC"),
                                                 net_legend_x = 0.1,
                                                 net_legend_y = 0.01)

img_width = 140 
p_name <- "supplementary_clin_int.pdf"
ggsave(file.path(supplementary_figure_path, p_name), clinical_cog_int$plot |> pad_plot(), width = img_width*scaling_factor, 
       height = img_width*0.6*scaling_factor, units = "mm", dpi = 300, #device = "pdf",
       bg = "white")

source_data <- clinical_cog_int$tmaps |> 
  select(-n, -model_formula) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient")))
write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))


######################################################################
# Supplementary figure for analysis in clinical group without mPACC
######################################################################

clinical_wo_cog <-  plot_gradient_relationships(biofinder_df |> filter(fmri_bl, diagnosis=="MCI" | diagnosis=="AD", !is.na(mPACC_v1)) |> 
                                                   mutate(`-mPACC_v1` = -mPACC_v1), 
                                                gradient_data = grad_df |> filter(study=="biofinder"), 
                                                gradients = c(1, 3),
                                                vect = TRUE,
                                                gradient_colors = gradient_cols,
                                                list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                mod_formula = formula(paste0(" ~ age + pathology_ad + sex + rsqa__MeanFD")),
                                                covariates = c("sex", "rsqa__MeanFD"),
                                                empty_row_height = -0.15,
                                                brain_title_size = 6,
                                                axis_text_size = 5,
                                                axis_title_size = 6,
                                                scatter_title_vjust = 0,
                                                base_size_ = 7,
                                                rasterize = TRUE,
                                                ggrastr_dpi = 300,
                                                add_shade = TRUE, 
                                                shade_alpha = 0.1,
                                                shade_size = 0.1,
                                                r_spin_size = 0.75,
                                                point_size = 0.05,
                                                point_alpha = 0.3,
                                                rectangle = TRUE,
                                                plt_title = c("Diagnosed MCI/AD (Ab+)"),
                                                plt_subtitle = TRUE,
                                                plot_net_legend = TRUE,
                                                brain_names = c("Age", "AD Pathology"),
                                                net_legend_x = 0.2,
                                                net_legend_y = 0.01)



img_width = 90 
p_name <- "supplementary_clin_without_cognition.pdf"
ggsave(file.path(supplementary_figure_path, p_name), clinical_wo_cog$plot |> pad_plot(), 
       width = img_width*scaling_factor, 
       height = img_width*0.85*scaling_factor, 
       units = "mm", dpi = 300, #device = "pdf", 
       bg = "white")

source_data <- clinical_wo_cog$tmaps |> 
  select(-n, -model_formula) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient")))
write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))


##############################################################################
# Supplementary figure for cognitively unimpaired without interaction
##############################################################################

health_cog_no_interact <-  plot_gradient_relationships(biofinder_df |> filter(fmri_bl, diagnosis=="Normal" | diagnosis=="SCD", abnorm_ab==0, !apoe4) |> 
                                                         mutate(`-mPACC_v1` = -mPACC_v1),
                                                       gradient_data = grad_df |> filter(study == "biofinder"), 
                                                       gradients = c(1, 3),
                                                       gradient_colors = gradient_cols,
                                                       list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                       mod_formula = formula(paste0(" ~ age + `-mPACC_v1` + pathology_ad + sex + rsqa__MeanFD")),
                                                       empty_row_height = -0.15,
                                                       brain_title_size = 6,
                                                       axis_text_size = 5,
                                                       axis_title_size = 6,
                                                       scatter_title_vjust = 0,
                                                       base_size_ = 7,
                                                       rasterize = TRUE,
                                                       ggrastr_dpi = 300,
                                                       add_shade = TRUE, 
                                                       shade_alpha = 0.1,
                                                       shade_size = 0.1,
                                                       r_spin_size = 0.65,
                                                       point_size = 0.05,
                                                       point_alpha = 0.3,
                                                       rectangle = FALSE,
                                                       plt_title = c("Cognitively unimpaired w/o APOE e4 (Ab-)"),
                                                       plt_subtitle = FALSE,
                                                       plot_net_legend = FALSE,
                                                       brain_names = c("Age", "-mPACC", "AD Pathology"),
                                                       net_legend_x = 0.2,
                                                       net_legend_y = 0.01)

n_health <- health_cog_no_interact$n
health_l_marg <- 0
p_health_cog_no_interact <- health_cog_no_interact$plot 

p_health_cog_no_interact <- p_health_cog_no_interact + 
  plot_annotation(title = bquote(
    "Cognitively unimpaired w/o APOE " * epsilon * 4 ~
      "(A" * beta^"-" * ") (N=" * .(n_health) * ")"
  ), 
  subtitle = expression(FC[parcel] ~ "~" ~ age + "-mPACC" + pathology + sex + motion),
  theme = theme(plot.subtitle = element_text(size = rel(0.9),
                                             family = "mono",
                                             face = "italic",
                                             hjust = 0, 
                                             vjust = -0.05, 
                                             margin = margin(l = health_l_marg, unit = "npc")), 
                plot.title = element_text(size = rel(1), 
                                          hjust =0, 
                                          margin = margin(l = health_l_marg, unit = "npc")))) 


plt_idx <- 3:5
#plt_idx <- plt_idx-1
p_health_cog_no_interact[[plt_idx[1]]] <- p_health_cog_no_interact[[plt_idx[1]]] + labs(title = "Age")
p_health_cog_no_interact[[plt_idx[2]]] <- p_health_cog_no_interact[[plt_idx[2]]] + labs(title = "-mPACC", subtitle = "(Inverted cognition)") + theme(plot.subtitle = element_text(hjust = 0.5))
p_health_cog_no_interact[[plt_idx[3]]] <- p_health_cog_no_interact[[plt_idx[3]]] + labs(title = "AD pathology")

img_width = 90
p_name <- "supplementary_health_no_interaction.pdf"
ggsave(file.path(supplementary_figure_path, p_name), p_health_cog_no_interact |> pad_plot(), 
       width = img_width*scaling_factor, 
       height = img_width*0.8*scaling_factor, 
       units = "mm", dpi = 300, #device = "pdf", 
       bg = "white")

source_data <- health_cog_no_interact$tmaps |> 
  select(-n, -model_formula) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient")))
write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))


##############################################################################
# Supplementary figure for cognitively unimpaired without pathology adjustment
##############################################################################


health_cog <-  plot_gradient_relationships(biofinder_df |> filter(fmri_bl, diagnosis=="Normal" | diagnosis=="SCD", abnorm_ab==0, !apoe4),
                                           gradient_data = grad_df |> filter(study == "biofinder"), 
                                           gradients = c(1, 3),
                                           gradient_colors = gradient_cols,
                                           list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                           vect = TRUE,
                                           mod_formula = formula(paste0(" ~ scale(age) * scale(-mPACC_v1) + sex + rsqa__MeanFD")),
                                           covariates = c("sex", "rsqa__MeanFD"),
                                           empty_row_height = -0.15,
                                           brain_title_size = 6,
                                           axis_text_size = 5,
                                           axis_title_size = 6,
                                           scatter_title_vjust = 0,
                                           base_size_ = 7,
                                           rasterize = TRUE,
                                           ggrastr_dpi = 300,
                                           add_shade = TRUE, 
                                           shade_alpha = 0.1,
                                           shade_size = 0.1,
                                           r_spin_size = 0.75,
                                           point_size = 0.05,
                                           point_alpha = 0.3,
                                           rectangle = FALSE,
                                           #plt_title = c("Cognitively unimpaired w/o APOE e4 (Ab-)"),
                                           plt_subtitle = FALSE,
                                           plot_net_legend = FALSE,
                                           brain_names = c("Age", "-mPACC", "AD Pathology"),
                                           net_legend_x = 0.2,
                                           net_legend_y = 0.01)

n_health <- health_cog$n
health_l_marg <- 0
p_health_cog <- health_cog$plot 

p_health_cog <- p_health_cog + 
  plot_annotation(title = bquote(
    "Cognitively unimpaired w/o APOE " * epsilon * 4 ~
      "(A" * beta^"-" * ") (N=" * .(n_health) * ")"
  ), 
  subtitle = expression(FC[parcel] ~ "~" ~ age * "×" * "-mPACC" + sex + motion),
  theme = theme(plot.subtitle = element_text(size = rel(1),
                                             family = "mono",
                                             face = "italic",
                                             hjust = 0, 
                                             vjust = -0.05, 
                                             margin = margin(l = health_l_marg, unit = "npc")), 
                plot.title = element_text(size = rel(1), 
                                          hjust =0, 
                                          margin = margin(l = health_l_marg, unit = "npc")))) 


plt_idx <- 3:5
#plt_idx <- plt_idx-1
p_health_cog[[plt_idx[1]]] <- p_health_cog[[plt_idx[1]]] + labs(title = "Age")
p_health_cog[[plt_idx[2]]] <- p_health_cog[[plt_idx[2]]] + labs(title = "-mPACC", subtitle = "(Inverted cognition)") + theme(plot.subtitle = element_text(hjust = 0.5))
p_health_cog[[plt_idx[3]]] <- p_health_cog[[plt_idx[3]]] + labs(title = "-mPACC×Age")


tau_in_health <- biofinder_df |> filter(fmri_bl, diagnosis %in% c("Normal", "SCD"), abnorm_ab==0, !apoe4) |> 
  pivot_longer(starts_with("cho"), values_to = "Tau PET SUVR", names_to = "region") |> 
  mutate(region = case_when(
    region == "cho12" ~ "Braak12",
    region == "cho34" ~ "Braak34",
    region == "cho56" ~ "Braak56"
  )) |> 
  ggplot(aes(`Tau PET SUVR`, fill = region)) +
  geom_density( alpha = 0.5, linewidth = 0.2) +
  labs(fill = "", y = "Density") +
  theme_bw(base_size = 7) +
  theme(legend.position = "bottom",
        plot.margin = unit(c(5, 5, 5, 5), "mm")) +
  ggsci::scale_fill_nejm()

library(cowplot)
health_no_pat <- ggdraw() +
  draw_plot(p_health_cog, x = 0, y = 0, height = 0.95, width = 0.6) +
  draw_plot_label("A", size = 7) +
  draw_plot(tau_in_health, x = 0.6, y = 0, height = 0.85, width = 0.4) +
  draw_plot_label("B", x = 0.6, size = 7)

img_width = 180
p_name <- "supplementary_health_no_pat_adj.pdf"
ggsave(file.path(supplementary_figure_path, p_name), health_no_pat |> pad_plot(), 
       width = img_width*scaling_factor, 
       height = img_width*0.45*scaling_factor, units = "mm", 
       dpi = 300, #device = "pdf",
       bg = "white")


source_data <- health_cog$tmaps |> 
  select(-n, -model_formula) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient")))
write_csv(source_data, file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelA.csv")))

biofinder_df |> filter(fmri_bl, diagnosis %in% c("Normal", "SCD"), abnorm_ab==0, !apoe4) |> 
  pivot_longer(starts_with("cho"), values_to = "Tau PET SUVR", names_to = "region") |> 
  mutate(region = case_when(
    region == "cho12" ~ "Braak12",
    region == "cho34" ~ "Braak34",
    region == "cho56" ~ "Braak56"
  )) |> select(diagnosis, `Tau PET SUVR`, region) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelB.csv")))



#################
# Pathology_score
#################

scale_factor <- 1
old <- theme_set(theme_bw(base_size = 7))
theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
             plot.background = element_rect(fill = "transparent", colour = NA),
             legend.background = element_rect(fill = "transparent", colour = NA),
             legend.box.background = element_rect(fill = "transparent", colour = NA))

legend_labs <-  c("Ab42/40", "Braak12", "Braak34", "Braak56")
biof_path_data <- biofinder_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |> 
  select(pathology_ad, ab_ratio, 
         starts_with("braak"), starts_with("cho")) |> 
  mutate(across(where(is.numeric) & !pathology_ad, scale)) 
ab_mean <- attr(biof_path_data$ab_ratio, c("scaled:center"))
ab_sd <- attr(biof_path_data$ab_ratio, c("scaled:scale"))
ab_pos <- (0.08 - ab_mean)/ab_sd

biof_path_plot <- biof_path_data |> 
  pivot_longer(-pathology_ad, names_to = "scaled_pat_measures", values_to = "value") %>% 
  ggplot(aes(pathology_ad, value, color = scaled_pat_measures)) +
  geom_smooth() +
  guides(color = guide_legend(
    nrow=1, byrow=TRUE,
    title.position="left", 
    label.position = "bottom",
    title.hjust = 0)) +
  geom_hline(yintercept = ab_pos, linetype = 6) +
  geom_text(inherit.aes = FALSE, 
            data = data.frame(x = 0.85, y = ab_pos),
            aes(x = x, y = y),
            label = expression("A" * beta^"+"), 
            parse = TRUE,
            nudge_y = 0.5, size = 5) +
  labs(x = "Pathology score", y = "Scaled value")+
  ggsci::scale_color_nejm(name = "Pathology", labels = legend_labs) +
  theme(legend.position = "bottom",
        legend.margin = margin(-4, 20, 0, 0),
        legend.box.margin = margin(-4, 20, 0, 0),
        #ggside.panel.scale = 0.2,
        legend.text = element_text(size = rel(0.6)),
        legend.title = element_text(size = rel(0.8)))  


dens <- biofinder_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |> 
  mutate(Diagnosis = ifelse(diagnosis == "Normal", "CN", diagnosis) |> factor(levels = c("CN", "SCD", "MCI", "AD")) ) |> 
  ggplot(aes(x = pathology_ad, fill = Diagnosis)) +
  geom_density(alpha = 0.5) +
  labs(x = "Pathology score", y = "Density") +
  ggsci::scale_fill_futurama() +
  guides(fill = guide_legend(theme = theme(
    legend.direction = "horizontal",
    legend.title.position = "left",
    legend.text.position = "bottom",
    # legend.text = element_text(hjust = 0.5, vjust = 1, #angle = 90
    #                            )
  ))) +
  theme(legend.position = "bottom", 
        legend.margin = margin(-4, 0, 0, 0),
        legend.box.margin = margin(-4, 0, 0, 0))

path_plot <- biof_path_plot + dens + plot_layout(axis_titles = "collect") 

# path_plot <- path_plot &
#   theme(text = element_text(size = 15))

img_width <- 180
p_name <- "path_plot_bf.pdf"
scale_factor = 1
ggsave(file.path(figure_path, p_name), path_plot |> pad_plot(),
       width = img_width*scale_factor, height = img_width/2.4*scale_factor, device = "pdf", units = "mm", bg = "transparent")

biofinder_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |> 
  mutate(Diagnosis = ifelse(diagnosis == "Normal", "CN", diagnosis) |> factor(levels = c("CN", "SCD", "MCI", "AD")) ) |> 
  select(Diagnosis, pathology_ad) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelB.csv")))


biof_path_data |> 
  pivot_longer(-pathology_ad, names_to = "scaled_pat_measures", values_to = "value")  |> 
  mutate(value = as.vector(value)) |> 
  select(scaled_pat_measures, value, pathology_ad) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_panelA.csv")))


###############
# Interaction 
###############

int_plot <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                        gradient_data = grad_df |> filter(study== "biofinder"), 
                                        gradients = c(1, 3),
                                        gradient_colors = gradient_cols,
                                        list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                        brain_title_size = 6,
                                        axis_text_size = 5,
                                        axis_title_size = 6,
                                        scatter_title_vjust = 0,
                                        rasterize = TRUE,
                                        ggrastr_dpi = 300,
                                        add_shade = TRUE, 
                                        shade_alpha = 0.1,
                                        shade_size = 0.1,
                                        point_size = 0.05,
                                        point_alpha = 0.2,
                                        l_width = 0.05,
                                        empty_row_height = -0.15,
                                        base_size_ = 7,
                                        brain_title_vjust = -3,
                                        r_spin_size = 0.775,
                                        vect = TRUE,
                                        plt_subtitle = TRUE,
                                        group_n_title = FALSE,
                                        rectangle = TRUE,
                                        plot_net_legend = TRUE,
                                        net_legend_x = 0.15,
                                        net_legend_y = 0.01,
                                        mod_formula = formula(" ~ scale(age) * pathology_ad + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        brain_names = c("Age (scaled)", "AD Pathology", "Age × ADpath"),
                                        plt_title = "Full sample")



####################
# Cognitive subtypes
####################



domains <- c("executive", "memory", "language", "visuospatial")
biofinder_cog_comp_res <- biofinder_cog_comp


ctrl_stats <- biofinder_cog_comp_res |>
  filter(diagnosis %in% c("Normal")) |>
  summarise(
    across(
      all_of(domains),
      list(
        mean = \(x) mean(x, na.rm = TRUE),
        sd   = \(x) sd(x, na.rm = TRUE)
      )
    )
  )

df_cog <- biofinder_cog_comp_res |>
  mutate(
    across(
      all_of(domains),
      ~ (.x - ctrl_stats[[paste0(cur_column(), "_mean")]]) /
        ctrl_stats[[paste0(cur_column(), "_sd")]]
    )
  )

impair_cutoff <- -1

df_cog <- df_cog |>
  mutate(
    across(
      all_of(domains),
      ~ .x < impair_cutoff,
      .names = "imp_{.col}"
    )
  ) |>
  rowwise() |>
  mutate(
    n_impairments = sum(c_across(starts_with("imp_"))),
    min_domain = which.min(c_across(all_of(domains))),
    cog_subtype = case_when(
      imp_memory ~ "amnestic",
      !imp_memory & n_impairments > 0 ~ "non-amnestic",
      TRUE ~ NA_character_
    )
  ) |>
  ungroup() |>
  mutate(cog_subtype = ifelse(diagnosis %in% c("MCI", "AD"), cog_subtype, NA),
         min_domain = ifelse(diagnosis %in% c("MCI", "AD"), min_domain, NA))


source("src/util.R")


plot_gam_gradient_effects <- function(
    gam_nonamnestic,
    gam_amnestic,
    gam_load,
    gam_eoad,
    grad_df,
    study = "biofinder",
    gradient = "gradient1",
    n_nonamnestic = NULL,
    n_amnestic = NULL,
    n_load = NULL,
    n_eoad = NULL
) {
  
  extract_gam_cors <- function(gam_obj, sample_name, n = NULL) {
    
    point_effects <- do.call(rbind, gam_obj$pat_derivs)
    
    pathology <- point_effects[nrow(point_effects), ]
    point_effects <- point_effects[-nrow(point_effects), , drop = FALSE]
    
    grad_vec <- grad_df |>
      dplyr::filter(study == !!study) |>
      dplyr::pull(.data[[gradient]])
    
    cors <- apply(point_effects, 2, function(x) cor(x, grad_vec, use = "pairwise.complete.obs"))
    
    label <- if (!is.null(n)) {
      paste0(sample_name, " (n = ", n, ")")
    } else {
      sample_name
    }
    
    tibble::tibble(
      pathology = as.numeric(pathology),
      gradient_correlation = as.numeric(cors),
      sample = label,
      sample_raw = sample_name,
      n = n,
      gradient = gradient
    )
  }
  
  source_data_age <- dplyr::bind_rows(
    extract_gam_cors(gam_load, "MCI/AD >65", n_load),
    extract_gam_cors(gam_eoad, "MCI/AD <65", n_eoad)
  ) |>
    dplyr::mutate(
      panel = "Full control sample, subsets of patients",
      sample = factor(
        sample,
        levels = unique(sample)
      )
    )
  
  source_data_cog <- dplyr::bind_rows(
    extract_gam_cors(gam_amnestic, "Amnestic", n_amnestic),
    extract_gam_cors(gam_nonamnestic, "Non-amnestic", n_nonamnestic)
  ) |>
    dplyr::mutate(
      panel = "Cognitive subtype",
      sample = factor(
        sample,
        levels = unique(sample)
      )
    )
  
  source_data <- dplyr::bind_rows(
    source_data_age,
    source_data_cog
  )
  
  age_plot <- source_data_age |>
    ggplot2::ggplot(ggplot2::aes(x = pathology, y = gradient_correlation)) +
    ggplot2::geom_line(show.legend = FALSE, linewidth = 0.5) +
    ggplot2::facet_wrap(~sample) +
    ggplot2::theme_bw(base_size = 7) +
    ggplot2::labs(
      y = "Gradient Correlation (r)",
      x = "AD Pathology"
    ) +
    ggplot2::ggtitle("Full control sample, subsets of patients") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = ggplot2::rel(1))
    )
  
  cog_plot <- source_data_cog |>
    ggplot2::ggplot(ggplot2::aes(x = pathology, y = gradient_correlation)) +
    ggplot2::geom_line(show.legend = FALSE, linewidth = 0.5) +
    ggplot2::facet_wrap(~sample) +
    ggplot2::theme_bw(base_size = 7) +
    ggplot2::labs(
      y = "Gradient Correlation (r)",
      x = "AD Pathology"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = ggplot2::rel(1))
    )
  
  plot <- patchwork::wrap_plots(
    list(age_plot, cog_plot),
    ncol = 1,
    guides = "collect",
    axis_titles = "collect",
    axes = "collect"
  ) +
    patchwork::plot_annotation(
      theme = ggplot2::theme(
        plot.background = ggplot2::element_rect(
          color = "black",
          linewidth = 0.5
        )
      )
    )
  
  list(
    plot = plot,
    source_data = source_data
  )
}

gam_preds_wo_amnestic <- gam_pred_nodes(df_cog |> filter(is.na(cog_subtype) | cog_subtype == "non-amnestic" ),
                                        fc_matrix = fc_measures_bf$affinity,
                                        roi_names = rois,
                                        id_var = "image_file",
                                        print_ticker = TRUE,
                                        model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))

gam_preds_only_amnestic <- gam_pred_nodes(df_cog |>  filter(is.na(cog_subtype) | cog_subtype == "amnestic" ),
                                          fc_matrix = fc_measures_bf$affinity,
                                          roi_names = rois,
                                          id_var = "image_file",
                                          print_ticker = TRUE,
                                          model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))

eoad_gam <- gam_pred_nodes(biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age>=65)),
                           fc_matrix = fc_measures_bf$affinity,
                           roi_names = rois, 
                           id_var = "image_file",
                           print_ticker = TRUE,
                           model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))


load_gam  <- gam_pred_nodes(biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age<=65)),
                            fc_matrix = fc_measures_bf$affinity,
                            roi_names = rois, 
                            id_var = "image_file",
                            print_ticker = TRUE,
                            model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))

# write_rds(list(gam_preds_wo_amnestic = gam_preds_wo_amnestic, 
#                gam_preds_only_amnestic = gam_preds_only_amnestic,
#                eoad_gam = eoad_gam,
#                load_gam = load_gam), "data/modelled_data/gams_supplementary.rds")


n_nonamn <- df_cog |> filter(is.na(cog_subtype) | cog_subtype == "non-amnestic" ) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
n_amn <- df_cog |>  filter(is.na(cog_subtype) | cog_subtype == "amnestic" )|> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()

n_eoad <- biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age>=65)) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
n_load <- biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age<=65)) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()


plot_gam_gradient_effects_res <- plot_gam_gradient_effects(
  gam_preds_wo_amnestic,
  gam_preds_only_amnestic,
  load_gam,
  eoad_gam,
  grad_df = grad_df,
  gradient = "gradient1",
  n_nonamnestic = n_nonamn$n,
  n_amnestic = n_amn$n,
  n_load = n_load$n,
  n_eoad = n_eoad$n
)

gam_plots <- plot_gam_gradient_effects_res$plot

age_plots <- ggdraw()+ 
  draw_plot(int_plot$plot, x = 0, y = 0, width = 0.60,  height = 1) + 
  draw_plot(gam_plots, x = 0.61, y = 0.0, width = 0.39, height = 1) 

scaling_factor <- 1
img_width <- 180
p_name <- "age_plot.pdf"
ggsave(file.path(supplementary_figure_path, p_name), age_plots |> pad_plot(),
       width = img_width*scaling_factor, height = img_width*0.5*scaling_factor, 
       units = "mm", dpi = 300, device = "pdf", bg = "white")

plot_gam_gradient_effects_res$source_data |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_LeftPanel.csv")))

int_plot$tmaps |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), "_RightPanel.csv")))


###########################
# Gradients from AD
###########################


con_all <- readRDS(file.path(atlas_dir, "average_connectome_ALL.rds"))

marg_gradients <- read_csv("data/gradients/margulies_gradients.csv")

grad_list_all <- get_gradients(connectome_ests = list(everyone = con_all),
                               reference_gradients = marg_gradients,
                               n_gradients = c(1,2,3),
                               threshold = 0.0,
                               similarity_method = "cosine",
                               on_affinity = FALSE,
                               method = "pca",
                               visualize = FALSE)

grad_all <- grad_list_all$gradients |> filter(study == "everyone") 

all_grad <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl, diagnosis %in% c("MCI", "AD")) |> 
                                          mutate(`-mPACC` = -mPACC_v1), 
                                        gradient_data = grad_all, 
                                        gradients = c(1, 3),
                                        gradient_colors = gradient_cols,
                                        list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                        brain_title_size = 6,
                                        axis_text_size = 5,
                                        axis_title_size = 5,
                                        scatter_title_vjust = 0,
                                        rasterize = TRUE,
                                        ggrastr_dpi = 300,
                                        add_shade = TRUE, 
                                        shade_alpha = 0.1,
                                        shade_size = 0.1,
                                        point_size = 0.05,
                                        point_alpha = 0.2,
                                        l_width = 0.05,
                                        empty_row_height = -0.25,
                                        base_size_ = 7,
                                        brain_title_vjust = -1,
                                        r_spin_size = 0.775,
                                        vect = TRUE,
                                        plt_subtitle = TRUE,
                                        group_n_title = FALSE,
                                        rectangle = TRUE,
                                        plot_net_legend = TRUE,
                                        net_legend_x = 0.2,
                                        net_legend_y = 0.015,
                                        mod_formula = formula(" ~ age + pathology_ad + `-mPACC` + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        plt_title = "Diagnosed MCI/AD, gradients derived from full sample ",
                                        brain_names = c("Age", "AD Pathology", "-mPACC"))

img_width = 88
p_name <- "all_grad_plot.pdf"

ggsave(file.path(supplementary_figure_path, p_name), all_grad$plot |> pad_plot(),
       width = img_width*scaling_factor,
       height = img_width*0.8*scaling_factor, 
       units = "mm", dpi = 300, device = "pdf", bg = "white")

all_grad$tmaps |> 
  left_join(grad_all |> filter(study %in% c("everyone")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))


###############################
# Across cohort pathology score and switched gradients
###############################

adni_df_patscore <- adni_df |> filter(fmri_bl) |> 
  drop_na(ab_ratio,
          braak1,
          braak34,
          braak56)

adni_mat <- adni_df_patscore |> filter(fmri_bl) |> 
  select(
    ab_ratio,
    braak1,
    braak34,
    braak56
  ) |> 
  rename(cho12 = braak1,
         cho34 = braak34,
         cho56 = braak56) |> as.matrix()

biofinder_df_patscore <- biofinder_df |>  filter(fmri_bl) |> 
  drop_na(
    ab_ratio,
    cho12,
    cho34,
    cho56
  )

bf_mat <- biofinder_df_patscore |> 
  select(ab_ratio,
         cho12,
         cho34,
         cho56) |> as.matrix()


withr::with_seed(12345, {
  space <- reduce_dimensionality(rbind(adni_mat, bf_mat), "euclidean")
  trajectory <- infer_trajectory(space)
})
trajectory <- reverse_trajectory(trajectory)

traject <- trajectory$time

adni_pat <- traject[1:nrow(adni_mat)]

bf_pat <- traject[(nrow(adni_mat)+1):(nrow(adni_mat) + nrow(bf_mat))]

biofinder_df_patscore$path_all <- bf_pat

adni_df_patscore$path_all <- adni_pat



bf_adni_grad <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                            gradient_data = grad_df |> filter(study== "adni"), 
                                            gradients = c(1, 3),
                                            gradient_colors = gradient_cols,
                                            list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                            brain_title_size = 6,
                                            axis_text_size = 5,
                                            axis_title_size = 5,
                                            scatter_title_vjust = 0,
                                            rasterize = TRUE,
                                            ggrastr_dpi = 300,
                                            add_shade = TRUE, 
                                            shade_alpha = 0.1,
                                            shade_size = 0.1,
                                            point_size = 0.05,
                                            point_alpha = 0.2,
                                            l_width = 0.05,
                                            empty_row_height = -0.25,
                                            base_size_ = 7,
                                            brain_title_vjust = -1,
                                            r_spin_size = 0.75,
                                            vect = TRUE,
                                            plt_subtitle = TRUE,
                                            group_n_title = FALSE,
                                            rectangle = TRUE,
                                            mod_formula = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"),
                                            covariates = c("sex", "rsqa__MeanFD"),
                                            plt_title = "Gradients derived from ADNI in BioFINDER sample",
                                            brain_names = c("Age", "AD Pathology"))

adni_bf_grad <- plot_gradient_relationships(adni_df %>% filter(fmri_bl), 
                                            gradient_data = grad_df |> filter(study== "biofinder"), 
                                            gradients = c(1, 3),
                                            gradient_colors = gradient_cols,
                                            list_of_parcel_data = list(nodal_affinity = fc_measures_adni$affinity),
                                            brain_title_size = 6,
                                            axis_text_size = 5,
                                            axis_title_size = 5,
                                            scatter_title_vjust = 0,
                                            rasterize = TRUE,
                                            ggrastr_dpi = 300,
                                            add_shade = TRUE, 
                                            shade_alpha = 0.1,
                                            shade_size = 0.1,
                                            point_size = 0.05,
                                            point_alpha = 0.2,
                                            l_width = 0.05,
                                            empty_row_height = -0.25,
                                            base_size_ = 7,
                                            brain_title_vjust = -1,
                                            r_spin_size = 0.75,
                                            vect = TRUE,
                                            plt_subtitle = TRUE,
                                            group_n_title = FALSE,
                                            rectangle = TRUE,
                                            id_var = "file_func",
                                            mod_formula = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD"),
                                            covariates = c("sex", "rsqa__MeanFD"),
                                            plt_title = "Gradients derived from BioFINDER in ADNI sample",
                                            brain_names = c("Age", "AD Pathology"))



bf_path_all <- plot_gradient_relationships(biofinder_df_patscore %>% filter(fmri_bl), 
                                           gradient_data = grad_df |> filter(study== "biofinder"), 
                                           gradients = c(1, 3),
                                           gradient_colors = gradient_cols,
                                           list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                           brain_title_size = 6,
                                           axis_text_size = 5,
                                           axis_title_size = 5,
                                           scatter_title_vjust = 0,
                                           rasterize = TRUE,
                                           ggrastr_dpi = 300,
                                           add_shade = TRUE, 
                                           shade_alpha = 0.1,
                                           shade_size = 0.1,
                                           point_size = 0.05,
                                           point_alpha = 0.2,
                                           l_width = 0.05,
                                           empty_row_height = -0.25,
                                           base_size_ = 7,
                                           brain_title_vjust = -1,
                                           r_spin_size = 0.75,
                                           vect = TRUE,
                                           plt_subtitle = TRUE,
                                           group_n_title = FALSE,
                                           rectangle = TRUE,
                                           mod_formula = formula(" ~ age + path_all + sex + rsqa__MeanFD"),
                                           covariates = c("sex", "rsqa__MeanFD"),
                                           plt_title = "BioFINDER, pathology score derived from both cohorts",
                                           brain_names = c("Age", "AD Pathology (all)"))

adni_path_all <- plot_gradient_relationships(adni_df_patscore %>% filter(fmri_bl), 
                                             gradient_data = grad_df |> filter(study== "adni"), 
                                             gradients = c(1, 3),
                                             gradient_colors = gradient_cols,
                                             list_of_parcel_data = list(nodal_affinity = fc_measures_adni$affinity),
                                             brain_title_size = 6,
                                             axis_text_size = 5,
                                             axis_title_size = 5,
                                             scatter_title_vjust = 0,
                                             rasterize = TRUE,
                                             ggrastr_dpi = 300,
                                             add_shade = TRUE, 
                                             shade_alpha = 0.1,
                                             shade_size = 0.1,
                                             point_size = 0.05,
                                             point_alpha = 0.2,
                                             l_width = 0.05,
                                             empty_row_height = -0.25,
                                             base_size_ = 7,
                                             brain_title_vjust = -1,
                                             r_spin_size = 0.75,
                                             vect = TRUE,
                                             plt_subtitle = TRUE,
                                             group_n_title = FALSE,
                                             rectangle = TRUE,
                                             id_var = "file_func",
                                             mod_formula = formula(" ~ age + path_all + sex + rsqa__MeanFD"),
                                             covariates = c("sex", "rsqa__MeanFD"),
                                             plt_title = "ADNI, pathology score derived from both cohorts",
                                             brain_names = c("Age", "AD Pathology (all)"))


switched_grads_all_path <- ggdraw() +
  draw_plot(bf_adni_grad$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=10)),
                                                              plot.subtitle = element_text(margin = margin(l=10)))),
            x = 0, y = 0.505, width = 0.495, height = 0.495) +
  draw_plot(adni_bf_grad$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=10)),
                                                              plot.subtitle = element_text(margin = margin(l=10))))
            , x = 0.505, y = 0.505, width = 0.495, height = 0.495) +
  draw_plot(bf_path_all$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=10)),
                                                            plot.subtitle = element_text(margin = margin(l=10))))
            , x = 0, y = 0, width = 0.495, height = 0.495) +
  draw_plot(adni_path_all$plot+ plot_annotation(theme = theme(plot.title = element_text(margin = margin(l=10)),
                                                              plot.subtitle = element_text(margin = margin(l=10))))
            , x = 0.505, y = 0, width = 0.495, height = 0.495) +
  draw_plot_label("A", size = 9) +
  draw_plot_label("B", x = 0.505, size = 9) +
  draw_plot_label("C", x = 0.0, y = 0.495, size = 9) +
  draw_plot_label("D", x = 0.505, y = 0.495, size = 9) 

scaling_factor <- 1
img_width <- 180
p_name <- "switched_grad_all_path.pdf"
ggsave(file.path(supplementary_figure_path, p_name), switched_grads_all_path |> pad_plot(),
       width = img_width*scaling_factor, 
       height = img_width*0.85*scaling_factor, 
       units = "mm", dpi = 300, device = "pdf", bg = "white")


bf_adni_grad$tmaps |> mutate(panel = "A") |> 
  bind_rows(adni_bf_grad$tmaps |> mutate(panel = "B")) |> 
  bind_rows(bf_path_all$tmaps |> mutate(panel = "C")) |> 
  bind_rows(adni_path_all$tmaps |> mutate(panel = "D")) |> 
  left_join(grad_df |> filter(study %in% c("biofinder")) |> 
              select(region, study, starts_with("gradient"))) |> 
  write_csv(file.path(source_figure_path, paste0(tools::file_path_sans_ext(p_name), ".csv")))



########    ###    ########  ##       ########  ######  
   ##      ## ##   ##     ## ##       ##       ##    ## 
   ##     ##   ##  ##     ## ##       ##       ##       
   ##    ##     ## ########  ##       ######    ######  
   ##    ######### ##     ## ##       ##             ## 
   ##    ##     ## ##     ## ##       ##       ##    ## 
   ##    ##     ## ########  ######## ########  ######  

library(finalfit)


descriptives <- biofinder_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |> 
  select(age, pathology_ad, sex, rsqa__MeanFD, diagnosis, abnorm_ab, education, cho12, cho34, cho56, ab_ratio) |> 
  mutate(abnorm_ab = ifelse(abnorm_ab == 1, "Positive", "Negative")) |> 
  rename(MeanFD = rsqa__MeanFD,
         Education = education,
         Aβ = abnorm_ab,
         Braak12 = cho12,
         Braak34 = cho34,
         Braak56 = cho56) |> 
  rename(Age = age) |> 
  rename(Pathology = pathology_ad) |> 
  mutate(sex = ifelse(sex == 1, "Female", "Male"),
         diagnosis = ifelse(diagnosis %in% c("AD", "MCI", "SCD", "Normal"), diagnosis, "Other"),
         diagnosis = ifelse(diagnosis %in% c("SCD", "Normal"), "CN", diagnosis), 
         diagnosis = ifelse(diagnosis %in% c("MCI"), "MCI_ab+", diagnosis), 
         diagnosis = factor(diagnosis, levels= c("CN", "MCI_ab+", "AD"))) |> 
  rename(Diagnosis = diagnosis) |> 
  rename(Sex = sex) 

explanatory = c("Age",  "Pathology", "Aβ", "ab_ratio", "Braak12", "Braak34", "Braak56", "MeanFD", "Education", "Sex")
dependent = 'Diagnosis'
descriptives |>
  summary_factorlist(dependent, explanatory, cont = "mean",
                     p=FALSE, add_dependent_label=FALSE,
                     dependent_label_prefix = "",
                     add_col_totals = TRUE,
                     digits = c(2, 2, 3, 1, 0)
  ) |> 
  mutate(Cohort = "BioFINDER",
         label = case_when(
           grepl("Braak", label) ~ paste0(label, " (SUVR)"),
           label == "MeanFD" ~ "Mean FD (mm)",
           label == "Pathology" ~ "Pathology Score",
           label == "Education" ~ "Education (Years)",
           label == "ab_ratio" ~ "Aβ42/40",
           TRUE ~ label
         )
         ) |> 
  relocate(Cohort)-> t1_bf



descriptives_adni <- adni_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |>
  select(age, pathology_ad, sex, diagnosis, rsqa__MeanFD, abnorm_ab, education_yrs, braak1, braak34, braak56, ab_ratio) |> 
  mutate(abnorm_ab = ifelse(abnorm_ab == 1, "Positive", "Negative")) |> 
  mutate(sex = ifelse(sex == 0, "Male", "Female")) |> 
  mutate(Diagnosis = factor(ifelse(diagnosis=="MCI",  "MCI_ab+", 
                                   ifelse(diagnosis=="Dementia",  "AD", diagnosis)), 
                            levels= c("CN", "MCI_ab+", "AD"))) |> 
  rename(Age = age,
         MeanFD = rsqa__MeanFD,
         Pathology = pathology_ad,
         Education = education_yrs,
         Aβ = abnorm_ab,
         Braak12 = braak1,
         Braak34 = braak34,
         Braak56 = braak56,
         Sex = sex)

dependent = 'Diagnosis'
descriptives_adni |>
  summary_factorlist(dependent, explanatory, cont = "mean",
                     p=FALSE, add_dependent_label=FALSE,
                     dependent_label_prefix = "",
                     add_col_totals = TRUE,
                     digits = c(2, 2, 3, 1, 0)
  ) |> 
  mutate(Cohort = "ADNI", 
         label = case_when(
           grepl("Braak", label) ~ paste0(label, " (SUVR)"),
           label == "MeanFD" ~ "Mean FD (mm)",
           label == "Pathology" ~ "Pathology Score",
           label == "Education" ~ "Education (Years)",
           label == "ab_ratio" ~ "Aβ42/40",
           TRUE ~ label
         )
         ) |> 
  relocate(Cohort)-> t1_adni

t1 <- rbind(t1_bf, t1_adni)

write_rds(t1, file.path(table_path, "CS_tbl1.rds"))
write_csv(t1, file.path(table_path, "CS_tbl1.csv"))

###########################################
# Longitudinal
###########################################

###########
# BioFINDER
###########

descriptives <- biofinder_df |> 
  drop_na(age, tau_pathology, sex, # mPACC_v1, 
          rsqa__MeanFD) |> 
  mutate(abnorm_ab = ifelse(abnorm_ab == 1, "Positive", "Negative")) |> 
  group_by(sid) |> 
  arrange(sid, fmri_date) |> 
  mutate(
    Diagnosis_BL = first(diagnosis),
    education_bl = first(education),
    education_follow = last(education),
    aβ_bl = first(abnorm_ab),
    Age_bl = first(age),
    Age_follow = last(age),
    path_bl = first(tau_pathology),
    path_follow = last(tau_pathology),
    Braak12_bl = first(cho12),
    Braak34_bl = first(cho34),
    Braak56_bl = first(cho56),
    Braak12_follow = last(cho12),
    Braak34_follow = last(cho34),
    Braak56_follow = last(cho56),
    
    interval = (difftime(last(fmri_date), first(fmri_date), units = "days")/365) |> as.numeric(),

    pathΔ = last(tau_pathology) - path_bl,
    time = (difftime(fmri_date, first(fmri_date), units = "days")/365) |> as.numeric(),
    n_visits = n()
  ) |>  
  fill(diagnosis, .direction = "downup") |> 
  ungroup() |>   
  filter(n()>1, .by = "sid") |> 
  mutate(sex = ifelse(sex == 1, "Female", "Male"),
         Diagnosis_BL = ifelse(Diagnosis_BL %in% c("AD", "MCI", "SCD", "Normal"), Diagnosis_BL, "Other"),
         Diagnosis_BL = ifelse(Diagnosis_BL %in% c("SCD", "Normal"), "CN", Diagnosis_BL), 
         Diagnosis_BL = ifelse(Diagnosis_BL %in% c("MCI"), "MCI_ab+", Diagnosis_BL), 
         Diagnosis_BL = factor(Diagnosis_BL, levels= c("CN", "MCI_ab+", "AD"))) |> 
  rename(Diagnosis_BL = Diagnosis_BL) |> 
  rename(Sex = sex) |> 
  filter(fmri_date == max(fmri_date), .by = "sid") 



longi_tbl1_bf <- descriptives |> 
  summary_factorlist("Diagnosis_BL", c(
    "n_visits",
    "interval",
    "pathΔ",
    "Sex", 
    "education",
    "Age_bl", 
    "path_bl", 
    "aβ_bl",
    "Braak12_bl",
    "Braak34_bl",
    "Braak56_bl",
    "Age_follow",
    "path_follow",
    "Braak12_follow",
    "Braak34_follow",
    "Braak56_follow"), cont = "mean",
                     p=FALSE, add_dependent_label=FALSE,
                     dependent_label_prefix = "",
                     add_col_totals = TRUE,
                     digits = c(2, 2, 3, 1, 0)
  ) |> 
  mutate(Cohort = "BioFINDER",
         label = ifelse(label == "", NA, label)) |> 
  fill(label, .direction = "down") |> 
  mutate(TimePoint = case_when(
    grepl("follow", label) ~ "Last follow up",
    grepl("bl", label) ~ "Baseline",
    TRUE ~ "Constant"
  ),
  TimePoint = as_factor(TimePoint)) |> 
  mutate(label = str_remove_all(label, "_follow"),
         label = str_remove_all(label, "_bl"),
         label = str_replace_all(label, "pathΔ", "ΔPathology score"),
         label = str_replace_all(label, "path", "Pathology score"),
         label = str_replace_all(label, "interval", "Follow up (Years)"),
         label = str_replace_all(label, "n_visits", "Number of visits"),
         label = str_replace_all(label, "education", "Education (Years)"),
         label = ifelse(grepl("Braak", label), paste0(label, " (SUVR)"), label)) |> 
  relocate(Cohort, TimePoint) 


longitudinal_table <- rbind(longi_tbl1_bf)


write_rds(longitudinal_table, file.path(table_path, "LT_tbl1.rds"))
write_csv(longitudinal_table, file.path(table_path, "LT_tbl1.csv"))


##############################
# Gradient sensitivity analyses
##############################

source("src/util.R")
if (from_start) {
  supp_table_grads <- write_gradient_supp_table(params) 
  write_rds(supp_table_grads, file.path(table_path, "supp_table.rds"))
}

########  ######## ##    ## ########  ######## ########     ########     ###    ########     ###    ######## ########  
##     ## ##       ###   ## ##     ## ##       ##     ##    ##     ##   ## ##   ##     ##   ## ##   ##       ##     ## 
##     ## ##       ####  ## ##     ## ##       ##     ##    ##     ##  ##   ##  ##     ##  ##   ##  ##       ##     ## 
########  ######   ## ## ## ##     ## ######   ########     ########  ##     ## ########  ##     ## ######   ########  
##   ##   ##       ##  #### ##     ## ##       ##   ##      ##        ######### ##        ######### ##       ##   ##   
##    ##  ##       ##   ### ##     ## ##       ##    ##     ##        ##     ## ##        ##     ## ##       ##    ##  
##     ## ######## ##    ## ########  ######## ##     ##    ##        ##     ## ##        ##     ## ######## ##     ##

source("src/render_paper.R")
