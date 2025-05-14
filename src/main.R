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


figure_path <- "paper/figures"
dir.create(figure_path, showWarnings = FALSE)
table_path <- "paper/tables"
dir.create(table_path, showWarnings = FALSE)

from_start <- FALSE
create_brain_permutations <- FALSE
extract_timeseries <- FALSE
real_data <- TRUE


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
    image_file = csv_rsqa__index,
    abnorm_ab = Abnormal_CSF_Ab42_Ab40_Ratio,
    sex = gender_baseline_variable,
    education = education_level_years_baseline_variable,
    diagnosis = diagnosis_baseline_variable,
    ab_ratio = CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022,
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
         image_file,
         has_longitudinal,
         age,
         education,
         diagnosis, 
         mPACC_v1, 
         apoe4,
         ab_ratio,
         abnorm_ab,
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
  write_rds(average_connectome, file.path(atlas_dir, "average_connectome_normalyoung.rds"))
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
  
  # This is max millimiter (of journals) in inches
  img_width = 180 / 25.4
  
  # Scale image to balance figure elements and fontsize, but need "shrink" image after
  # so the fontsize chosen should be large enough to shrink it by the scaling factor
  scaling_factor <-  3
  magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")
  p_name <- "gradient_param_comparison_bf.png"
  
  print("Writing supplementary figure 9")
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
adni_df_ <- adni_df |> select(ID, file_func, DX, age, sex, amyloid_status,
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
  mutate(EXAMDATE_func = as.Date(EXAMDATE_func),
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
}

rsqa <- lapply(adni_fd, function(x) c(rsqa__MeanFD = x |> pull(1) |> mean(), rsqa__MaxFD = x |> pull(1) |> max()))
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
  inner_join(rsqa_fd) 

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

source("src/methods_figure.R")

########################
# Main figures
########################

# This script contains an array of functions that creates the figures in the paper
source("src/util_vis.R")

# These parameters set the final figure width, and the scaling to use to 
# get elements and font sizes aligned
img_width <-  180 / 25.4
scaling_factor <-  3
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")

################
# FIGURE 1
################


fig_one <- figure_one(subject_data = biofinder_df,
                      subject_data_replication = adni_df |> filter(fmri_bl),
                      measures_list =list(nodal_affinity = fc_measures_bf$affinity), 
                      measures_list_replication = list(nodal_affinity = fc_measures_adni$affinity),
                      gradients_df = grad_df |> filter(study == "biofinder"),
                      gradients_df_replication = grad_df |> filter(study == "adni"),
                      selected_gradients = c(1, 3), 
                      empt_row_height = -0.1,
                      b_size = 21,
                      draw_size = 21,
                      plot_title_size = 0.9,
                      axes_title_size = 0.9,
                      r2_sizing1 = 6,
                      boxed = FALSE,
                      split = TRUE)


p_name <- "figure1.png"
ggsave(file.path(figure_path, p_name), fig_one[[1]],
       width = img_width*scaling_factor, height = img_width*0.5*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


###########################
# FIGURE 2 (non-linear)
###########################

source("src/util.R")
source("src/util_vis.R")

gam_preds_affinity <- gam_pred_nodes(biofinder_df |> filter(fmri_bl),
                                     fc_matrix = fc_measures_bf$affinity,
                                     roi_names = rois, 
                                     id_var = "image_file",
                                     print_ticker = TRUE,
                                     model_formula = formula("FC ~ s(age) + s(pathology_ad) + sex + rsqa__MeanFD"))

nonlin_p <- plot_gams_v1(gam_predictions = gam_preds_affinity, 
                         grad_df = grad_df |> filter(study == "biofinder"),
                         biofinder_data = biofinder_df |> filter(fmri_bl), scale_fac = 3)

p_name <- "figure2.png"
img_width <- 90/25.4
scale_factor = 3
magick_geom_scaling <- paste0(100*0.66, "%x", 100*0.66, "%")
ggsave(file.path(figure_path, p_name), nonlin_p, 
       width = img_width*scale_factor, 
       height = img_width*scale_factor*1.3, 
       units = "in", dpi = 300,
       bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)



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

# Very case specific function for creating figure 3
longitudindal_figs <- longitudinal_and_window_analysis(long_df = long_bf_, 
                                                       b_size = 18,
                                                       run_windowing = FALSE
                                                         #ifelse(from_start, TRUE, FALSE) 
                                                       )

img_width = 180 / 25.4
scaling_factor <- 2
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")

p_name <- "figure3.png"
ggsave(file.path(figure_path, p_name), longitudindal_figs[["main_fig"]],
       width = img_width*scaling_factor, height = img_width*1*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


p_name <- "supplementary_windowing.png"
ggsave(file.path(figure_path, p_name), longitudindal_figs[["supp_fig"]],
       width = img_width*scaling_factor, height = img_width*0.5*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


######################
# FIGURE 4 (cognition)
######################

img_width = 180 / 25.4
scaling_factor <-  3
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")

fig4 <- figure_one(subject_data = biofinder_df,
                   subject_data_replication = adni_df |> filter(fmri_bl),
                   measures_list =list(nodal_affinity = fc_measures_bf$affinity), 
                   measures_list_replication = list(nodal_affinity = fc_measures_adni$affinity),
                   gradients_df = grad_df |> filter(study == "biofinder"),
                   gradients_df_replication = grad_df |> filter(study == "adni"),
                   selected_gradients = c(1, 3), 
                   tag_size = 18,
                   draw_size = 18,
                   empt_row_height = -0.1,
                   b_size = 18,
                   plot_title_size = 0.85,
                   axes_title_size = 0.85,
                   r2_sizing2 = 4.4,
                   boxed = FALSE,
                   split = TRUE)

p_name <- "figure4.png"
ggsave(file.path(figure_path, p_name), fig4[[2]], 
       width = img_width*scaling_factor, height = img_width*0.4*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


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

bf_dx <-  plot_gradient_relationships(biofinder_df %>% 
                                        filter(fmri_bl) %>%
                                        mutate(motion = rsqa__MeanFD
                                        ) %>% 
                                        mutate(diagnosis = ifelse(diagnosis %in% c("Normal", "SCD"), "CN/SCD", diagnosis)) %>% 
                                        mutate(diagnosis=factor(diagnosis, levels = c("CN/SCD", "MCI", "AD"))), 
                                      gradient_data = grad_df %>% filter(study=="biofinder"), 
                                      gradients = c(1, 2, 3),
                                      vect = TRUE,
                                      gray_out = FALSE,
                                      r2_size = rel(4.5),
                                      gradient_colors = gradient_cols,
                                      list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                      mod_formula = formula(paste0(" ~ age + diagnosis + sex + rsqa__MeanFD")),
                                      covariates = c("sex", "rsqa__MeanFD"),
                                      filter_criteria = quo(),
                                      show_networks = FALSE,
                                      tag_prefix = "",
                                      tag_sep = "",
                                      layout_construction = "horizontal",
                                      include_gradient_plots = TRUE,
                                      right_term_side = FALSE,
                                      plt_title = "",
                                      cache_runs = FALSE)

p_dx <- bf_dx$plot + 
  plot_annotation(title = paste0("BioFINDER", " (N=", bf_dx$n, ")"), subtitle = expression(FC[parcel] ~ "~" ~ age + diagnosis + sex + motion),
                  theme = theme(plot.subtitle = element_text(hjust = 0, vjust = -0.05,
                                                             family = "mono",
                                                             margin = margin(l = 0, unit = "npc")), 
                                plot.background = element_rect(colour = "black"),
                                plot.title.position = "plot",
                                plot.title = element_text( hjust =0, margin = margin(l = 0, unit = "npc"))
                  )) &
  theme(text = element_text(size = 15),
        plot.tag = element_blank())         


p_dx[[4]] <- p_dx[[4]] + ggtitle("Age")
p_dx[[5]] <- p_dx[[5]] + ggtitle("MCI vs. CN")
p_dx[[6]] <- p_dx[[6]] + ggtitle("AD vs. CN")

img_width <- 180/25.4/2

p_name <- "supplementary_diagnosis.png"
ggsave(file.path(figure_path, p_name), p_dx, width = img_width*3, 
       bg = "white", height = img_width*0.9*3, units = "in", dpi = 300, device = "png")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, "33%x33%")
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


########################################################
# Supplementary figure with gradient 2
########################################################

fig_one <- figure_one(subject_data = biofinder_df,
                      subject_data_replication = adni_df |> filter(fmri_bl),
                      measures_list =list(nodal_affinity = fc_measures_bf$affinity), 
                      measures_list_replication = list(nodal_affinity = fc_measures_adni$affinity),
                      gradients_df = grad_df |> filter(study == "biofinder"),
                      gradients_df_replication = grad_df |> filter(study == "adni"),
                      selected_gradients = c(2), 
                      empt_row_height = -0.1,
                      b_size = 21,
                      draw_size = 21,
                      plot_title_size = 0.9,
                      axes_title_size = 0.9,
                      r2_sizing1 = 6,
                      boxed = FALSE,
                      split = TRUE)


p_name <- "supplementary_gradient2.png"
ggsave(file.path(figure_path, p_name), fig_one[[1]],
       width = img_width*scaling_factor, height = img_width*0.35*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)

########################################################
# Supplementary figure within network affinity
########################################################

img_width = 180 / 25.4
scaling_factor <-  3
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")


within_supp <- figure_one(subject_data = biofinder_df,
                          subject_data_replication = adni_df |> filter(fmri_bl),
                          measures_list =list(within_network_affinity = fc_measures_bf$affinity_within), 
                          measures_list_replication = list(within_network_affinity = fc_measures_adni$affinity_within),
                          gradients_df = grad_df |> filter(study == "biofinder"),
                          gradients_df_replication = grad_df |> filter(study == "adni"),
                          selected_gradients = c(1, 3),
                          tag_size = 18,
                          draw_size = 18,
                          empt_row_height = -0.1,
                          b_size = 18,
                          plot_title_size = 0.90,
                          axes_title_size = 0.85,
                          r2_sizing2 = 4.5,
                          r2_sizing1 = 6,
                          boxed = TRUE,
                          split = FALSE)

p_name <- "within_supplementary.png"
ggsave(file.path(figure_path, p_name), within_supp, width = img_width*scaling_factor, height = img_width*0.8*scaling_factor,
       units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


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
                           tag_size = 18,
                           draw_size = 18,
                           empt_row_height = -0.1,
                           b_size = 18,
                           plot_title_size = 0.90,
                           axes_title_size = 0.85,
                           r2_sizing2 = 4.5,
                           r2_sizing1 = 6,
                           boxed = TRUE,
                           split = FALSE)

p_name <- "between_supplementary.png"
ggsave(file.path(figure_path, p_name), between_supp, width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)



########################################################
# Supplementary figure with connectivity strength
########################################################


strength_supp <- figure_one(subject_data = biofinder_df,
                            subject_data_replication = adni_df |> filter(fmri_bl),
                            measures_list =list(nodal_strength = fc_measures_bf$strength), 
                            measures_list_replication = list(nodal_strength = fc_measures_adni$strength),
                            gradients_df = grad_df |> filter(study == "biofinder"),
                            gradients_df_replication = grad_df |> filter(study == "adni"),
                            selected_gradients = c(1, 3),
                            tag_size = 18,
                            draw_size = 18,
                            empt_row_height = -0.1,
                            b_size = 18,
                            plot_title_size = 0.90,
                            axes_title_size = 0.85,
                            r2_sizing2 = 4.5,
                            r2_sizing1 = 6,
                            boxed = TRUE,
                            split = FALSE)

p_name <- "supplementary_strength.png"
ggsave(file.path(figure_path, p_name), strength_supp, width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


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
                            tag_size = 18,
                            draw_size = 18,
                            empt_row_height = -0.1,
                            b_size = 18,
                            plot_title_size = 0.90,
                            axes_title_size = 0.85,
                            r2_sizing2 = 4.5,
                            r2_sizing1 = 6,
                            boxed = TRUE,
                            split = FALSE)

p_name <- "supplementary_affinity_no_thresh_correlation.png"
ggsave(file.path(figure_path, p_name), aff_no_thresh, width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


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
                                                 logistic_fit = FALSE,
                                                 covariates = c("sex", "rsqa__MeanFD"),
                                                 filter_criteria = quo(),
                                                 r2_size = 5,
                                                 show_networks = FALSE,
                                                 tag_prefix = "",
                                                 tag_sep = "",
                                                 layout_construction = "horizontal",
                                                 include_gradient_plots = TRUE,
                                                 right_term_side = FALSE,
                                                 plt_title = "",
                                                 cache_runs = FALSE)

n_clin <- clinical_cog_int$n
clin_l_marg = 0
p_clinical_cog <-
  clinical_cog_int$plot + plot_annotation(title = paste0("Diagnosed MCI/AD (Ab+)", " (N=", n_clin, ")"), subtitle = expression(FC[parcel] ~ "~" ~ age + pathology ~ "× -mPACC" + sex + motion),
                                          theme = theme(plot.subtitle = element_text(hjust = 0, vjust = -0.05,
                                                                                     family = "mono",
                                                                                     margin = margin(l = clin_l_marg, unit = "npc")), 
                                                        plot.title.position = "plot",
                                                        plot.title = element_text( hjust =0, margin = margin(l = clin_l_marg, unit = "npc"))
                                          )
  ) &
  theme(plot.tag = element_blank(),
        text = element_text(size = 16)
  )

plt_idx <- 4:6
p_clinical_cog[[plt_idx[1]]] <- p_clinical_cog[[plt_idx[1]]] + labs(title = "AD Pathology")
p_clinical_cog[[plt_idx[2]]] <- p_clinical_cog[[plt_idx[2]]] + labs(title = "-mPACC", subtitle = "(Inverted cognition)") + theme(plot.subtitle = element_text(hjust = 0.5))
p_clinical_cog[[plt_idx[3]]] <- p_clinical_cog[[plt_idx[3]]] + labs(title = "AD Pathology × -mPACC") + theme(plot.subtitle = element_text(hjust = 0.5))


img_width = 140 / 25.4
p_name <- "supplementary_clin_int.png"
ggsave(file.path(figure_path, p_name), p_clinical_cog, width = img_width*scaling_factor, 
       height = img_width*0.6*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)

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
                                                 logistic_fit = FALSE,
                                                 covariates = c("sex", "rsqa__MeanFD"),
                                                 filter_criteria = quo(),
                                                 r2_size = 6,
                                                 show_networks = FALSE,
                                                 tag_prefix = "",
                                                 tag_sep = "",
                                                 layout_construction = "horizontal",
                                                 include_gradient_plots = TRUE,
                                                 right_term_side = FALSE,
                                                 plt_title = "",
                                                 cache_runs = FALSE)

n_clin <- clinical_wo_cog$n
clin_l_marg = 0
p_clinical_wo_cog <-
  clinical_wo_cog$plot + plot_annotation(title = paste0("Diagnosed MCI/AD (Ab+)", " (N=", n_clin, ")"), 
                                          subtitle = expression(FC[parcel] ~ "~" ~ age + pathology + sex + motion),
                                          theme = theme(plot.subtitle = element_text(hjust = 0, vjust = -0.05,
                                                                                     family = "mono",
                                                                                     margin = margin(l = clin_l_marg, unit = "npc")), 
                                                        plot.title.position = "plot",
                                                        plot.title = element_text( hjust =0, margin = margin(l = clin_l_marg, unit = "npc"))
                                          )
  ) &
  theme(plot.tag = element_blank(),
        text = element_text(size = 14)
  )

plt_idx <- 3:4
p_clinical_wo_cog[[plt_idx[1]]] <- p_clinical_wo_cog[[plt_idx[1]]] + labs(title = "Age")
p_clinical_wo_cog[[plt_idx[2]]] <- p_clinical_wo_cog[[plt_idx[2]]] + labs(title = "AD Pathology")

img_width = 90 / 25.4
p_name <- "supplementary_clin_without_cognition.png"
ggsave(file.path(figure_path, p_name), p_clinical_wo_cog, width = img_width*scaling_factor, height = img_width*0.85*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


##############################################################################
# Supplementary figure for cognitively unimpaired without interaction
##############################################################################

health_cog_no_interact <-  plot_gradient_relationships(biofinder_df |> filter(fmri_bl, diagnosis=="Normal" | diagnosis=="SCD", abnorm_ab==0, !apoe4) |> 
                                                         mutate(`-mPACC_v1` = -mPACC_v1),
                                                       gradient_data = grad_df |> filter(study == "biofinder"), 
                                                       gradients = c(1, 3),
                                                       gradient_colors = gradient_cols,
                                                       r2_size = rel(7),
                                                       base_size_ = 18,
                                                       list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                       vect = TRUE,
                                                       mod_formula = formula(paste0(" ~ age + `-mPACC_v1` + pathology_ad + sex + rsqa__MeanFD")),
                                                       logistic_fit = FALSE,
                                                       covariates = c("sex", "rsqa__MeanFD"),
                                                       filter_criteria = quo(),
                                                       show_networks = FALSE,
                                                       tag_prefix = "",
                                                       tag_sep = "",
                                                       layout_construction = "horizontal",
                                                       plot_spacing = 0.2,
                                                       include_gradient_plots = TRUE,
                                                       right_term_side = FALSE,
                                                       plt_title = "",
                                                       cache_runs = FALSE)

n_health <- health_cog_no_interact$n
health_l_marg <- 0
p_health_cog_no_interact <- health_cog_no_interact$plot &
  theme(plot.tag = element_blank(),
        title  = element_text(size = rel(0.8))) 

p_health_cog_no_interact <- p_health_cog_no_interact + 
  plot_annotation(title = paste0("Cognitively unimpaired, no APOE e4, Ab-", 
                                 " (N=", n_health, ")"), 
                  subtitle = expression(FC[parcel] ~ "~" ~ age + "-mPACC" + pathology + sex + motion),
                  theme = theme(plot.subtitle = element_text(size = rel(1),
                                                             family = "mono",
                                                             face = "italic",
                                                             hjust = 0, 
                                                             vjust = -0.05, 
                                                             margin = margin(l = health_l_marg, unit = "npc")), 
                                plot.title = element_text(size = rel(1.2), 
                                                          hjust =0, 
                                                          margin = margin(l = health_l_marg, unit = "npc")))) 


plt_idx <- 3:5
#plt_idx <- plt_idx-1
p_health_cog_no_interact[[plt_idx[1]]] <- p_health_cog_no_interact[[plt_idx[1]]] + labs(title = "Age")
p_health_cog_no_interact[[plt_idx[2]]] <- p_health_cog_no_interact[[plt_idx[2]]] + labs(title = "-mPACC", subtitle = "(Inverted cognition)") + theme(plot.subtitle = element_text(hjust = 0.5))
p_health_cog_no_interact[[plt_idx[3]]] <- p_health_cog_no_interact[[plt_idx[3]]] + labs(title = "AD pathology")

img_width = 150 / 25.4
p_name <- "supplementary_health_no_interaction.png"
ggsave(file.path(figure_path, p_name), p_health_cog_no_interact, width = img_width*scaling_factor, height = img_width*0.6*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)


##############################################################################
# Supplementary figure for cognitively unimpaired without pathology adjustment
##############################################################################


health_cog <-  plot_gradient_relationships(biofinder_df |> filter(fmri_bl, diagnosis=="Normal" | diagnosis=="SCD", abnorm_ab==0, !apoe4),
                                           gradient_data = grad_df |> filter(study == "biofinder"), 
                                           gradients = c(1, 3),
                                           gradient_colors = gradient_cols,
                                           r2_size = rel(6),
                                           base_size_ = 18,
                                           list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                           vect = TRUE,
                                           mod_formula = formula(paste0(" ~ scale(age) * scale(-mPACC_v1) + sex + rsqa__MeanFD")),
                                           logistic_fit = FALSE,
                                           covariates = c("sex", "rsqa__MeanFD"),
                                           filter_criteria = quo(),
                                           show_networks = FALSE,
                                           tag_prefix = "",
                                           tag_sep = "",
                                           layout_construction = "horizontal",
                                           plot_spacing = 0.2,
                                           include_gradient_plots = TRUE,
                                           right_term_side = FALSE,
                                           plt_title = "",
                                           cache_runs = FALSE)

n_health <- health_cog$n
health_l_marg <- 0
p_health_cog <- health_cog$plot &
  theme(plot.tag = element_blank(),
        title  = element_text(size = rel(0.8))) 

p_health_cog <- p_health_cog + 
  plot_annotation(title = paste0("Cognitively unimpaired, no APOE e4, Ab-", 
                                 " (N=", n_health, ")"), 
                  subtitle = expression(FC[parcel] ~ "~" ~ age * "×" * "-mPACC" + sex + motion),
                  theme = theme(plot.subtitle = element_text(size = rel(1.2),
                                                             family = "mono",
                                                             face = "italic",
                                                             hjust = 0, 
                                                             vjust = -0.05, 
                                                             margin = margin(l = health_l_marg, unit = "npc")), 
                                plot.title = element_text(size = rel(1.2), 
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
  geom_density( alpha = 0.5) +
  labs(fill = "", y = "Density") +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom",
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
  ggsci::scale_fill_nejm()

library(cowplot)
health_no_pat <- ggdraw() +
  draw_plot(p_health_cog, x = 0, y = 0, height = 0.95, width = 0.7) +
  draw_plot_label("A", size = 21) +
  draw_plot(tau_in_health, x = 0.7, y = 0, height = 0.85, width = 0.3) +
  draw_plot_label("B", x = 0.7, size = 21)

img_width = 180 / 25.4
p_name <- "supplementary_health_no_pat_adj.png"
ggsave(file.path(figure_path, p_name), health_no_pat, width = img_width*scaling_factor, height = img_width*0.45*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)





#################
# Pathology_score
#################

scale_factor <- 4
old <- theme_set(theme_bw(base_size = 5*scale_factor))
theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
             plot.background = element_rect(fill = "transparent", colour = NA),
             legend.background = element_rect(fill = "transparent", colour = NA),
             legend.box.background = element_rect(fill = "transparent", colour = NA))

legend_labs <-  c("Ab42/40", "Braak12", "Braak34", "Braak56")
biof_path_plot <- biofinder_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |> 
  select(pathology_ad, ab_ratio, 
         starts_with("braak"), starts_with("cho")) |> 
  mutate(across(where(is.numeric) & !pathology_ad, scale)) |> 
  pivot_longer(-pathology_ad, names_to = "scaled_pat_measures", values_to = "value") |> 
  ggplot(aes(pathology_ad, value, color = scaled_pat_measures)) +
  geom_smooth() +
  guides(color = guide_legend(theme = theme(
    legend.direction = "horizontal",
    legend.title.position = "left",
    legend.text.position = "bottom",
    # legend.text = element_text(hjust = 0.5, vjust = 1, #angle = 90
    #                            )
  ))) +
  labs(x = "Pathology score", y = "Scaled value")+
  ggsci::scale_color_nejm(name = "Pathology", labels = legend_labs) +
  theme(legend.position = "bottom")  

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
  theme(legend.position = "bottom")

path_plot <- biof_path_plot + dens + plot_layout(axis_titles = "collect") 

# path_plot <- path_plot &
#   theme(text = element_text(size = 15))

img_width <- 80/25.4
ggsave(file.path(figure_path, "path_plot_bf.png"), path_plot, width = img_width*scale_factor, height = img_width/2.4*scale_factor, bg = "transparent")

################
# Main overlaid
################

source("src/util_vis.R")
img_width = 120 / 25.4
scaling_factor <-  3
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")

net_main_res <- overlaid_main_results(biofinder_df |> filter(fmri_bl), fc_measures_bf$affinity) 


p_name <- "main_res_with_overlaid_nets.png"
ggsave(file.path(figure_path, p_name), net_main_res,
       width = img_width*scaling_factor, height = img_width*0.7*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")

img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)



########    ###    ########  ##       ########  ######  
   ##      ## ##   ##     ## ##       ##       ##    ## 
   ##     ##   ##  ##     ## ##       ##       ##       
   ##    ##     ## ########  ##       ######    ######  
   ##    ######### ##     ## ##       ##             ## 
   ##    ##     ## ##     ## ##       ##       ##    ## 
   ##    ##     ## ########  ######## ########  ######  

library(finalfit)


descriptives <- biofinder_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |> 
  select(age, pathology_ad, sex, rsqa__MeanFD, diagnosis, abnorm_ab, education, cho12, cho34, cho56) |> 
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

explanatory = c("Age",  "Pathology", "Aβ", "Braak12", "Braak34", "Braak56", "MeanFD", "Education", "Sex")
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
           TRUE ~ label
         )
         ) |> 
  relocate(Cohort)-> t1_bf



descriptives_adni <- adni_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |>
  select(age, pathology_ad, sex, diagnosis, rsqa__MeanFD, abnorm_ab, education_yrs, braak1, braak34, braak56) |> 
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
           TRUE ~ label
         )
         ) |> 
  relocate(Cohort)-> t1_adni

t1 <- rbind(t1_bf, t1_adni)

write_rds(t1, file.path(table_path, "CS_tbl1.rds"))

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

quarto::quarto_render("paper/fc_changes_paper.qmd")
