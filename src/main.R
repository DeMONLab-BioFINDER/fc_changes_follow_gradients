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

library(tidyverse)
library(broom)
library(SCORPIUS)
library(conflicted)
library(magick)
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
extract_timeseries <- FALSE

ts_dir <- "data/bf_src_data/timeseries"
connectome_dir <- "data/bf_src_data/connectomes"
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


if (from_start) {
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

##################
# Images
##################
if (extract_timeseries){
  
  require(reticulate)
  nilearn <- import("nilearn", convert = FALSE)
  np <- import("numpy", convert = FALSE)
  
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
  img_dir <- "data/bf_src_data/images"
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



##################
# Subject data
##################

sub_data <- read_csv('data/bf_src_data/data_bf.csv', show_col_types = FALSE)
cog_data <- read_csv("data/bf_src_data/cognitive_BF2.csv") |> 
  dplyr::select(uid, sid, Visit, mmse_score, mPACC_v1, mPACC_v2, mPACC_v3)


biofinder_df_ <- sub_data
temp <- biofinder_df_ |> inner_join(sub_data |> select(uid, csv_rsqa__index, rsqa__fd_max, rsqa__MeanFD))  |> 
  select(sid, csv_rsqa__index, rsqa__fd_max, rsqa__MeanFD) |> drop_na() 

biofinder_df_ <- biofinder_df_ |> 
  inner_join(temp)
rm(temp)

biofinder_df__ <- biofinder_df_ |> 
  left_join(cog_data) |>
  mutate(
    fmri_date = as.Date(str_split_i(csv_rsqa__index, "_", 5), format = "%Y%m%d"),
    image_file = csv_rsqa__index,
    tau_file = csv_btnic_sr_mr_fs__index,
    abnorm_ab = Abnormal_CSF_Ab42_Ab40_Ratio,
    sex = gender_baseline_variable,
    education = education_level_years_baseline_variable,
    diagnosis = diagnosis_baseline_variable,
    ab_ratio = CSF_Ab42_Ab40_ratio_imputed_Elecsys_2020_2022,
    apoe4 = grepl("4", as.character(apoe_genotype_baseline_variable)),
    ptau217 = CSF_ptau217_pgml_Lilly_2019,
    ptau181 = CSF_ptau181_pgml_Lilly_2019,
    alpha_syn = CSF_SynucleinAlpha_RTQuIC_2022,
    cho12 = tnic_cho_com_I_II,
    cho34 = tnic_cho_com_III_IV,
    cho56 = tnic_cho_com_V_VI,
    ento_tau_mean = (tnic_sr_mr_fs_ctx_lh_entorhinal + tnic_sr_mr_fs_ctx_rh_entorhinal)/2,
    hippo_tau_mean = (tnic_sr_mr_fs_Left_Hippocampus + tnic_sr_mr_fs_Right_Hippocampus)/2,
    inftemp_tau_mean = (tnic_sr_mr_fs_ctx_rh_inferiortemporal + tnic_sr_mr_fs_ctx_lh_inferiortemporal)/2,
    precuneus_tau_mean = (tnic_sr_mr_fs_ctx_lh_precuneus + tnic_sr_mr_fs_ctx_rh_precuneus)/2,
    postcing_tau_mean = (tnic_sr_mr_fs_ctx_lh_posteriorcingulate + tnic_sr_mr_fs_ctx_rh_posteriorcingulate)/2,
    motion_filter = (rsqa__fd_max < 3 & rsqa__MeanFD < 0.3)
  ) |> 
  mutate(fmri_bl = fmri_date == min(fmri_date),
         has_longitudinal = n()>1, .by = "sid") |> 
  select(sid, 
         uid,
         Visit, 
         fmri_bl,
         fmri_date,
         image_file,
         tau_file,
         has_longitudinal,
         age,
         education,
         diagnosis, 
         mPACC_v1, 
         mPACC_v2, 
         mPACC_v3,
         mmse_score,
         apoe4,
         ab_ratio,
         abnorm_ab,
         ptau217,
         ptau181,
         alpha_syn,
         sex, 
         rsqa__MeanFD,
         rsqa__fd_max,
         motion_filter,
         cho12,
         cho34,
         cho56,
         ento_tau_mean,
         hippo_tau_mean,
         inftemp_tau_mean,
         precuneus_tau_mean,
         postcing_tau_mean
  ) |> 
  arrange(sid, Visit) |> 
  group_by(sid) |> 
  filter(
    (any(abnorm_ab == 1, na.rm = TRUE) | !all(diagnosis == "MCI" & abnorm_ab == 0, na.rm = TRUE)) &
      any(diagnosis %in% c("AD", "MCI", "SCD", "Normal") | is.na(diagnosis), na.rm = TRUE)
  ) |> 
  filter(any(diagnosis != "MSA", na.rm = TRUE)) |> 
  filter(any(diagnosis != "PD", na.rm = TRUE)) |> 
  filter(any(diagnosis != "PPA_NOS", na.rm = TRUE)) |> 
  ungroup() 


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
  dplyr::select(uid,
                ab_ratio,
                cho12,
                cho34,
                cho56
  ) |>  drop_na()

withr::with_seed(123132, {
  space <- reduce_dimensionality(as.matrix(patvars |> select(-uid)), "euclidean")
  trajectory <- infer_trajectory(space)
})
traject <- trajectory$time
#biofinder_df__ |> left_join(patvars |> mutate(pathology_ad = traject) |> select(uid, pathology_ad)) |> ggplot(aes(pathology_ad, fill = diagnosis)) +geom_density(alpha = 0.5)
biofinder_df__ <- biofinder_df__ |> left_join(patvars |> mutate(pathology_ad = traject) |> select(uid, pathology_ad)) 

patvars_tau <- biofinder_df__ |> 
  dplyr::select(uid,
                cho12,
                cho34,
                cho56
  ) |>  drop_na()
withr::with_seed(123, {
  space <- reduce_dimensionality(as.matrix(patvars_tau |> select(-uid)), "euclidean")
  trajectory <- infer_trajectory(space)
})
traject <- trajectory$time
#biofinder_df__ |> left_join(patvars_tau |> mutate(tau_path = traject) |> select(uid, tau_path)) |> ggplot(aes(tau_path, fill = diagnosis)) +geom_density(alpha = 0.5)
biofinder_df <- biofinder_df__ |> left_join(patvars_tau |> mutate(tau_pathology = traject) |> select(uid, tau_pathology)) |> 
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
rm(biofinder_df__, space, trajectory, patvars_tau, patvars, cog_data, biofinder_df_)


########  ########  ########          ######     ###    ##        ######  
##     ## ##     ## ##               ##    ##   ## ##   ##       ##    ## 
##     ## ##     ## ##               ##        ##   ##  ##       ##       
########  ########  ######   ####### ##       ##     ## ##       ##       
##        ##   ##   ##               ##       ######### ##       ##       
##        ##    ##  ##               ##    ## ##     ## ##       ##    ## 
##        ##     ## ########          ######  ##     ## ########  ######  

#########################################################
## Extract connectomes
#########################################################


if (from_start) {
  
  ts_files <- list.files(ts_dir, pattern = "\\.rds$", full.names = TRUE)
  ts_list <- lapply(ts_files, readRDS)
  names(ts_list) <- tools::file_path_sans_ext(basename(ts_files))
  
  # Take only the timeseries from subjects of interest
  timeseries <- ts_list[biofinder_df$image_file]
  
  dir.create(connectome_dir, showWarnings = FALSE)
  already_proc_conn <- list.files(connectome_dir) |> str_remove_all(".rds")
  proc_conn <- names(timeseries)[!(names(timeseries) %in% already_proc_conn)]
  pb = txtProgressBar(min = 0, 
                      max = ifelse(length(proc_conn) != 0, length(proc_conn), 1),
                      initial = 0)
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
        p  },
      warning = function(cond) {
        print(cond)
      })
    i = i + 1
    setTxtProgressBar(pb,i)
  }
  close(pb)
}

success_vec <- list.files(connectome_dir) |> str_remove_all(".rds")
biofinder_df <- biofinder_df |> filter(image_file %in% success_vec)
biof_motion_unfilt <- biof_motion_unfilt |> filter(image_file %in% success_vec)


dir.create(clean_dir, showWarnings = FALSE)
write_rds(biofinder_df, file.path(clean_dir, "biofinder_df.rds"))
write_rds(biof_motion_unfilt, file.path(clean_dir, "biofinder_motion_unfiltered.rds"))

biofinder_df <- biofinder_df |> filter(motion_filter)

if (from_start) {
  
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
  
  # images have been accessed from https://neurovault.org/collections/1598/
  nifti_images <- list.files(gradient_dir, pattern = "\\.nii.gz$")
  
  
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



########  ######## ########  ##       ####  ######     ###    ######## ####  #######  ##    ## 
##     ## ##       ##     ## ##        ##  ##    ##   ## ##      ##     ##  ##     ## ###   ## 
##     ## ##       ##     ## ##        ##  ##        ##   ##     ##     ##  ##     ## ####  ## 
########  ######   ########  ##        ##  ##       ##     ##    ##     ##  ##     ## ## ## ## 
##   ##   ##       ##        ##        ##  ##       #########    ##     ##  ##     ## ##  #### 
##    ##  ##       ##        ##        ##  ##    ## ##     ##    ##     ##  ##     ## ##   ### 
##     ## ######## ##        ######## ####  ######  ##     ##    ##    ####  #######  ##    ## 


df_full <- readxl::read_xlsx("data/adni_src_data/ADNI_selected_Jake_full.xlsx")
df_selected <- readxl::read_xlsx("data/adni_src_data/ADNI_selected_Jake.xlsx")
adni_df_ <- df_full |> select(any_of(colnames(df_selected)), age, sex, APOE4_alleles, CSF_AlphaSyn_seeding,
                               contains("braak"),
                               contains("entorhinal"), 
                               contains("inferiortemporal"), contains("precuneus")
) 

rm(df_full, df_selected)

id_ses <- c()
for (i in adni_df_$file_func) {
  id_ses <- c(id_ses, str_split_1(str_split_i(i, "/", 10), "_")[1:2] |> paste0(collapse = "_") )
}

adni_df__ <- adni_df_ |> 
  mutate(id_ses = id_ses,
         EXAMDATE_func = as.Date(EXAMDATE_func)
  ) |> 
  relocate(id_ses) |> 
  mutate(fmri_bl = EXAMDATE_func == min(EXAMDATE_func), .by = ID) |> 
  mutate(
    apoe4 = APOE4_alleles > 0,
    diagnosis = DX,
    fmri_date = EXAMDATE_func,
    asyn_pos = CSF_AlphaSyn_seeding == "positive",
    abnorm_ab = ifelse(amyloid_status == "Ab.neg", 0, 1),
    braak1 = tau.SUVR.DK.braak1,
    braak34 = tau.SUVR.DK.braak34,
    braak56 = tau.SUVR.DK.braak56,
    prec_amy = centiloid.amyloid.SUVR.DK.lprecuneus + centiloid.amyloid.SUVR.DK.rprecuneus,
    ento_amy = centiloid.amyloid.SUVR.DK.lentorhinal + centiloid.amyloid.SUVR.DK.rentorhinal,
    inftemp_amy = centiloid.amyloid.SUVR.DK.linferiortemporal + centiloid.amyloid.SUVR.DK.rinferiortemporal,
    ento_tau_mean = (tau.SUVR.DK.lentorhinal + tau.SUVR.DK.rentorhinal)/2,
    inftemp_tau_mean = (tau.SUVR.DK.linferiortemporal + tau.SUVR.DK.rinferiortemporal)/2,
    precuneus_tau_mean = (tau.SUVR.DK.lprecuneus + tau.SUVR.DK.rprecuneus)/2,
    ab_ratio = ABETA42/ABETA40,
    has_ab_csf = !is.na(ab_ratio),
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
  space <- reduce_dimensionality(as.matrix(patvars |> select(-id_ses)), "euclidean")
  trajectory <- infer_trajectory(space)
})
traject <- trajectory$time

adni_df___ <- adni_df__ |> left_join(patvars |> mutate(pathology_ad = traject) |> select(id_ses, pathology_ad))

patvars_tau <- adni_df__ |> 
  dplyr::select(id_ses,
                braak1,
                braak34,
                braak56
  )  |>  drop_na()

withr::with_seed(12345, {
  space <- reduce_dimensionality(as.matrix(patvars_tau |> select(-id_ses)), "euclidean")
  trajectory <- infer_trajectory(space)
})
traject <- trajectory$time
adni_df____ <- adni_df___ |> left_join(patvars_tau |> mutate(tau_pathology = traject) |> select(id_ses, tau_pathology)) 


rm(adni_df_, adni_df__, adni_df___)

############
# Timeseries
############

fd_files <- list.files("data/adni_src_data/timeseries", pattern = "displacement\\.xlsx")

adni_fd <- list()
for (file in fd_files) {
  sub_and_ses <- paste0(str_split_i(file, "_", 1), "_", str_split_i(file, "_", 2))
  adni_fd[sub_and_ses] <- read_xlsx(paste0("data/adni_src_data/timeseries/", file))
}

rsqa <- lapply(adni_fd, function(x) c(rsqa__MeanFD = mean(x), rsqa__MaxFD = max(x)))
rsqa_fd <- do.call(rbind, rsqa) |> as_tibble(rownames = NA) |> rownames_to_column("id_ses")


if (from_start) {
  ts_dir <- "data/adni_src_data/timeseries"
  timeseries_files <- list.files(ts_dir, pattern = "timeseries\\.xlsx", full.names = TRUE)
  
  sub_and_ses <- lapply(timeseries_files, function(file) {
    file |>  
      str_remove(paste0(ts_dir, "/")) |> 
      str_remove("_task.*")
    }
  ) |> unlist()
  
  adni_timeseries <- pblapply(timeseries_files, function(file) {
    read_xlsx(file, progress = FALSE) |> 
      column_to_rownames("ROI") |> 
      t() |> 
      as.matrix()
  })
  names(adni_timeseries) <- sub_and_ses
  
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
    fd_filt <- set_false_window(fd<0.5)
    scrubbed_time_series[[file]] <- adni_timeseries[[file]][fd_filt, ]
    
  }
  frame_length <- sapply(scrubbed_time_series, nrow)
  timeseries <- scrubbed_time_series[frame_length>100]
}


##############
# Connectomes
##############
if (from_start) {
  connectome_dir <- "data/adni_src_data/connectomes"
  dir.create(connectome_dir, showWarnings = FALSE)
  already_proc_conn <- list.files(connectome_dir) |> str_remove_all(".rds")
  proc_conn <- names(timeseries)[!(names(timeseries) %in% already_proc_conn)]
  
  pb = txtProgressBar(min = 0, max = length(proc_conn), style = 3)
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


success_vec <- list.files("data/adni_src_data/connectomes") |> str_remove_all(".rds")
adni_df <- adni_df____ |> filter(id_ses %in% success_vec) |> 
  inner_join(rsqa_fd) 

adni_df_unfilt <- adni_df |> mutate(motion_filter = (rsqa__MeanFD<0.3 & rsqa__MaxFD<3))
adni_df <- adni_df_unfilt |>  filter(motion_filter)


write_rds(adni_df, file.path(clean_dir, "adni_df.rds"))
write_rds(adni_df_unfilt, file.path(clean_dir,"adni_motion_unfilt.rds"))

if (from_start) {
  con_cube_adni <- array(dim = c(length(rois), 
                                 length(rois), 
                                 length(success_vec)), 
                         dimnames = list(rois, rois, success_vec))
  pb = txtProgressBar(min = 0, max = length(success_vec), initial = 0, style = 3)
  i = 0
  for (img_file in success_vec){
    tryCatch(
      {
        con_cube_adni[, , img_file] <- read_rds(file.path(connectome_dir, paste0(img_file, ".rds")))
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
  adni_fc <- fc_strength(con_cube_adni, 
                         roi_names = rois, 
                         threshold = 0, 
                         replace = 0)
  
  adni_affinity <- get_affinity(con_cube_adni, 
                                similarity_method = "cosine",
                                roi_names = rois, 
                                threshold = 0.75, 
                                atlas = yeo_msk)
  
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
  yh_filt <- adni_df |> filter(fmri_bl, abnorm_ab==0, !apoe4, DX == "CN") |> pull(id_ses)
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

# Source the methods_figure.R scrip to generate all plots for the methods figure.
# However, the figures have then been manually put together in inkscape so I will leave this
# out of the main pipeline 

# source("src/methods_figure.R")

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

###### I AM HERE


fig_one <- figure_one(subject_data = biofinder_df,
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
                          r2_sizing1 = 5.75,
                          boxed = FALSE,
                          split = TRUE)


p_name <- "figure1.png"
ggsave(file.path(figure_path, p_name), fig_one[[1]],
       width = img_width*scaling_factor, height = img_width*0.45*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
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
                         biofinder_data = biofinder_df, scale_fac = 3)

p_name <- "figure2.png"
img_width <- 90/25.4
scale_factor = 3
magick_geom_scaling <- paste0(100/scale_factor, "%x", 100/scale_factor, "%")
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
                                                       run_windowing = ifelse(from_start, TRUE, FALSE) )

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

fig2 <- figure_one(subject_data = biofinder_df,
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
ggsave(file.path(figure_path, p_name), fig2[[2]], 
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
# Supplementary figure within network affinity
########################################################

img_width = 180 / 25.4
scaling_factor <-  3
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")


within_supp <- figure_one(subject_data = biofinder_df,
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
                                                 r2_size = 3.8,
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
        #text = element_text(size = text_size)
  )

plt_idx <- 4:6
p_clinical_cog[[plt_idx[1]]] <- p_clinical_cog[[plt_idx[1]]] + labs(title = "AD Pathology")
p_clinical_cog[[plt_idx[2]]] <- p_clinical_cog[[plt_idx[2]]] + labs(title = "-mPACC", subtitle = "(Inverted cognition)") + theme(plot.subtitle = element_text(hjust = 0.5))
p_clinical_cog[[plt_idx[3]]] <- p_clinical_cog[[plt_idx[3]]] + labs(title = "AD Pathology × -mPACC") + theme(plot.subtitle = element_text(hjust = 0.5))


img_width = 90 / 25.4
p_name <- "supplementary_clin_int.png"
ggsave(file.path(figure_path, p_name), p_clinical_cog, width = img_width*scaling_factor, height = img_width*0.6*scaling_factor, units = "in", dpi = 300, device = "png", bg = "white")
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
  select(age, pathology_ad, sex, rsqa__MeanFD, diagnosis) |> 
  rename(MeanFD = rsqa__MeanFD) |> 
  rename(Age = age) |> 
  rename(Pathology = pathology_ad) |> 
  mutate(sex = ifelse(sex == 1, "Female", "Male"),
         diagnosis = ifelse(diagnosis %in% c("AD", "MCI", "SCD", "Normal"), diagnosis, "Other"),
         diagnosis = ifelse(diagnosis %in% c("SCD", "Normal"), "CN", diagnosis), 
         diagnosis = ifelse(diagnosis %in% c("MCI"), "MCI_ab+", diagnosis), 
         diagnosis = factor(diagnosis, levels= c("CN", "MCI_ab+", "AD"))) |> 
  rename(Diagnosis = diagnosis) |> 
  rename(Sex = sex) 

explanatory = c("Age", "Pathology", "MeanFD", "Sex")
dependent = 'Diagnosis'
descriptives |>
  summary_factorlist(dependent, explanatory, cont = "median",
                     p=FALSE, add_dependent_label=FALSE,
                     dependent_label_prefix = "",
                     add_col_totals = TRUE,
                     digits = c(2, 2, 3, 1, 0)
  ) |> 
  mutate(Cohort = "BioFINDER") |> 
  relocate(Cohort)-> t1_bf



descriptives_adni <- adni_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |>
  select(age, pathology_ad, sex, diagnosis, rsqa__MeanFD) |> 
  mutate(sex = ifelse(sex == 0, "Male", "Female")) |> 
  mutate(Diagnosis = factor(ifelse(diagnosis=="MCI",  "MCI_ab+", 
                                   ifelse(diagnosis=="Dementia",  "AD", diagnosis)), 
                            levels= c("CN", "MCI_ab+", "AD"))) |> 
  rename(MeanFD = rsqa__MeanFD) |> 
  rename(Age = age) |> 
  rename(Pathology = pathology_ad) |> 
  rename(Sex = sex) 

explanatory = c("Age", "Pathology", "Sex", "MeanFD")
dependent = 'Diagnosis'
descriptives_adni |>
  summary_factorlist(dependent, explanatory, cont = "median",
                     p=FALSE, add_dependent_label=FALSE,
                     dependent_label_prefix = "",
                     add_col_totals = TRUE,
                     digits = c(2, 2, 3, 1, 0)
  ) |> 
  mutate(Cohort = "ADNI") |> 
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
  group_by(sid) |> 
  arrange(sid, fmri_date) |> 
  mutate(
    Diagnosis_BL = first(diagnosis),
    ab_bl = first(abnorm_ab),
    age_bl = min(age),
    path_bl = first(tau_pathology),
    # path_bl_grp = cut(path_bl, quantile(path_bl, probs = seq(0, 1, length.out = 4))),
    mPACC_bl = first(mPACC_v1),
    #ageΔ  = age - age_bl,
    pathΔ = last(tau_pathology) - path_bl,
    time = (difftime(fmri_date, first(fmri_date), units = "days")/365) |> as.numeric(),
    time_diff = as.numeric(difftime(fmri_date, lag(fmri_date), units = "days")) / 365,
    mean_time = mean(time_diff, na.rm = TRUE),
    n_visits = n()
  ) |>  
  fill(diagnosis, .direction = "downup") |> 
  ungroup() |>   
  filter(n()>1, .by = "sid") |> 
  rename(Age = age) |> 
  mutate(sex = ifelse(sex == 1, "Female", "Male"),
         Diagnosis_BL = ifelse(Diagnosis_BL %in% c("AD", "MCI", "SCD", "Normal"), Diagnosis_BL, "Other"),
         Diagnosis_BL = ifelse(Diagnosis_BL %in% c("SCD", "Normal"), "CN", Diagnosis_BL), 
         Diagnosis_BL = ifelse(Diagnosis_BL %in% c("MCI"), "MCI_ab+", Diagnosis_BL), 
         Diagnosis_BL = factor(Diagnosis_BL, levels= c("CN", "MCI_ab+", "AD"))) |> 
  rename(Diagnosis_BL = Diagnosis_BL) |> 
  rename(Sex = sex) |> 
  filter(fmri_date == max(fmri_date), .by = "sid") 

longi_tbl1_bf <- descriptives |> 
  summary_factorlist("Diagnosis_BL", c("Age", "mean_time", "pathΔ", "Sex"), cont = "median",
                     p=FALSE, add_dependent_label=FALSE,
                     dependent_label_prefix = "",
                     add_col_totals = TRUE,
                     digits = c(2, 2, 3, 1, 0)
  ) |> 
  mutate(Cohort = "BioFINDER") |> 
  relocate(Cohort)


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
