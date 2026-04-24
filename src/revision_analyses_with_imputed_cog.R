###########################
# Cognition
###########################

library(psych)
library(mice)



###############
# Linear stuff
###############
# new_cog <- read_delim("stuff_for_revisions/jonathan__20250926_104238.csv")
# 
# cogs <- new_cog |> select(aqt_form:`letter_s_45-60`)

#cog <- read_delim("stuff_for_revisions/cognitive_tests-2025-08-29.csv") 
#cog <- read_delim("stuff_for_revisions/cognitive_tests-2025-08-29.csv") 
nic <- read_delim("stuff_for_revisions/nicola__20240326_143046.csv") 
# nic <- read_delim("stuff_for_revisions/jonathan__20250926_104238.csv")
# cog <- nic |> select(1:5, aqt_form:`letter_s_45-60`,
#                      # adas_delayed_word_recall, adas_immediate_word_recall_average,
#                      # letter_s, trailmaking_a, trailmaking_b, symbol_digit,
#                      # animal_fluency, vosp_cube, bnt_15_2,
#                      cognitive_test_date
# )

cog <- nic |> select(1:5, mPACC_v1:fcsrt_immediate,
                     # adas_delayed_word_recall, adas_immediate_word_recall_average,
                     # letter_s, trailmaking_a, trailmaking_b, symbol_digit,
                     # animal_fluency, vosp_cube, bnt_15_2,
                     cognitive_test_date
)

# biofinder_cog |> 
#   select(
#     sid,
#     diagnosis,
#     adas_delayed_word_recall,
#     adas_immediate_word_recall_average,
#     symbol_digit,
#     trailmaking_a, 
#     trailmaking_b,
#     letter_s,
#     animal_fluency,
#     vosp_cube,
#     bnt_15_2,
#     aqt_color_form
#   )
# 
# biofinder_cog %>%
#   drop_na(pathology_ad) |> 
#   select(
#     diagnosis,
#     adas_delayed_word_recall,
#     adas_immediate_word_recall_average,
#     symbol_digit,
#     trailmaking_a, 
#     trailmaking_b,
#     letter_s,
#     animal_fluency,
#     vosp_cube,
#     bnt_15_2,
#     aqt_color_form
#   ) %>%
#   mutate(has_missing = if_any(-diagnosis, is.na)) %>%
#   summarise(
#     n_rows = n(),
#     n_with_missing = sum(has_missing),
#     .by = diagnosis
#   ) |>  drop_na() |> arrange(n_with_missing)



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
# cat_vars  <- c("sex", "diagnosis")
# num_vars  <- c("age", "education", "pathology_ad", cog_vars)

library(mice)

# ini  <- mice(imp_df, maxit = 0, print = FALSE)
# meth <- ini$method
# pred <- ini$predictorMatrix
# 
# # ID
# meth["sid"] <- ""
# pred["sid", ] <- 0
# pred[, "sid"] <- 0
# 
# # Not imputing these (rows = 0), but allowing them as predictors (columns untouched)
# no_impute <- c("age", "sex", "education", "diagnosis", "pathology_ad", "mPACC_v1")
# pred[no_impute, ] <- 0
# meth[no_impute] <- ""
# 
# # Impute cog vars using selected predictors
# pred[cog_vars, c("age","sex","education","diagnosis","pathology_ad", "mPACC_v1")] <- 1
# pred[cog_vars, cog_vars] <- 1
# diag(pred) <- 0
# # Methods for cog vars
# meth[cog_vars] <- "pmm"
# 
# 
# 

# --- initialise -------------------------------------------------------
ini  <- mice(imp_df, maxit = 0, print = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

# --- ID: never imputed, never used -----------------------------------
meth[id_vars] <- ""
pred[id_vars, ] <- 0
pred[, id_vars] <- 0

# --- variables that should NOT be imputed (rows = 0) ------------------
no_impute_vars <- c(
  "sex",
  "diagnosis",
  "age",
  "education",
  "mPACC_v1",
  "pathology_ad"
)

pred[no_impute_vars, ] <- 0

# --- predictors for cognitive variables -------------------------------
cog_predictors <- c(
  "age",
  "sex",
  "education",
  "pathology_ad",
  "mPACC_v1"
)

pred[cog_vars, cog_predictors] <- 1
pred[cog_vars, cog_vars]       <- 1

# --- no self-prediction ----------------------------------------------
diag(pred) <- 0


# # 
# ini <- mice(imp_df, maxit = 0, print = FALSE)
# 
# meth <- ini$method
# pred <- ini$predictorMatrix
# 
# meth[id_vars] <- ""
# pred[id_vars, ] <- 0
# pred[, id_vars] <- 0
# 
# pred[c("sex", "diagnosis"), ] <- 0
# pred["age", ]        <- 0
# pred["education", ]     <- 0
# ##### Maybe remove
# pred["mPACC_v1", ]     <- 0
# pred["pathology_ad", ]     <- 0
# 
# pred[cog_vars, c(
#   "age",
#   "sex",
#   "education",
#   "pathology_ad",
#   "mPACC_v1"
# )] <- 1
# 
# pred[cog_vars, cog_vars] <- 1
# diag(pred) <- 0


# # --- init -------------------------------------------------------------
# ini  <- mice(imp_df, maxit = 0, print = FALSE)
# meth <- ini$method
# pred <- ini$predictorMatrix
# 
# # --- helpers ----------------------------------------------------------
# set_rows0 <- function(m, vars) { m[vars, ] <- 0; m }
# set_cols0 <- function(m, vars) { m[, vars] <- 0; m }
# 
# # --- configuration ----------------------------------------------------
# id_vars        <- "sid"
# no_impute_vars <- c("age", "sex", "education", "diagnosis", "pathology_ad", "mPACC_v1")
# cog_pred_vars  <- c("age", "sex", "education", "pathology_ad", "mPACC_v1")
# 
# # --- ID: never imputed, never used -----------------------------------
# meth[id_vars] <- ""
# pred <- pred |> set_rows0(id_vars) |> set_cols0(id_vars)
# 
# # --- variables we don’t impute (but may still be used as predictors) --
# pred <- pred |> set_rows0(no_impute_vars)
# 
# # --- impute cognitive tests using selected predictors + other tests ---
# pred[cog_vars, cog_pred_vars] <- 1
# pred[cog_vars, cog_vars]      <- 1
# 
# # --- final sanity: no self-prediction --------------------------------
# diag(pred) <- 0


# md.pattern(imp_df)
# 
# imp_df |> summarise(
#   across(everything(),
#          ~ sum(is.na(.x))/nrow(imp_df),
#          .names = "{.col}"
#   )
# )

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

#domain_scores_cu <- rename_factors_by_corr(om_factors, refs)


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
# 
# 
# mu <- colMeans(cog_mat_cu, na.rm = TRUE)
# sdv <- apply(cog_mat_cu, 2, sd, na.rm = TRUE)
# 
# test_z <- sweep(cog_mat_ad, 2, mu, "-")
# test_z <- sweep(test_z, 2, sdv, "/")

# mu  <- colMeans(cog_mat_cu, na.rm=TRUE)
# sdv <- apply(cog_mat_cu, 2, sd, na.rm=TRUE)
# Z_train <- scale(cog_mat_cu, center = mu, scale = sdv)
# W <- domain_scores_cu$weights
# Z_all <- scale(cog_mat_ad, center = mu, scale = sdv)
# domain_scores <- Z_all %*% W
# domain_scores <- domain_scores[, 1:5]
#colnames(domain_scores) <- c("global_cog", domains)

# domain_scores_cu <- factor.scores(cog_mat_cu, f = cog_fac$schmid$sl, method = "regression")
# domain_scores_ad <- factor.scores(cog_mat_ad, f = cog_fac$schmid$sl, method = "regression") 
# 
# domain_scores <- factor.scores(cog_mat_ad, f = cog_fac$schmid$sl, method = "tenBerge") 


# domain_scores_ad$scores[, 3] |>  hist()
# domain_scores_cu$scores[, 3] |>  hist()

# Check so that we get correct factor scores:
# scores_train_recomputed <- factor.scores(cog_mat_cu, f = cog_fac$schmid$sl, method = "regression")
# cor(scores_train_recomputed$scores[, 1:5], om_factors)
# Diagonal should be 1

# domain_scores <- domain_scores_ad$scores[, 1:5]
# colnames(domain_scores) <- c("global_cog", domains)
# domain_scores_cu <- rename_factors_by_corr(domain_scores_cu$scores[, 1:5], refs)
# domain_scores <- rbind(domain_scores_cu, domain_scores_ad$scores[, 1:5])


biofinder_cog_comp <- biofinder_cog |> 
  inner_join(domain_scores |> as_tibble(rownames = NA) |> rownames_to_column("sid")) 


cog_diff_plot <- function(b_size = 16, subt_size = rel(0.75)) {
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
           term = recode(term,
                         "scale(age)" = "Age",
                         "scale(pathology_ad)" = "Pathology"))
  
  # Plot
  mem_plot <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
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
           term = recode(term,
                         "scale(age)" = "Age",
                         "scale(pathology_ad)" = "Pathology"))
  
  # Plot
  exec_plot <- ggplot(coef_df, aes(x = term, y = estimate)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.1) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(y = "Executive", x = "") +
    theme_bw(base_size = b_size) 
  
  n = length(fit_exec$residuals)
  
  library(patchwork)
  cog_eff_plot <- wrap_plots(list(mem_plot, exec_plot), ncol = 1, guides = "collect", axes = "collect")  +
    plot_annotation(
      title = paste0("Effects on memory/executive function (n = ", n, ")"),
      subtitle = "Memory ~ scale(age) + scale(AD pathology) + sex \nExecutive ~ scale(age) + scale(AD pathology) + sex",
      theme = theme(plot.background = element_rect(color = "black", fill = NA))) &
    theme(axis.text.x = element_text(size = rel(1.5)),
          axis.title.y = element_text(size = rel(1.15)),
          plot.subtitle = element_text(family = "mono", size = subt_size))
  cog_eff_plot
}

cog_diff <- cog_diff_plot(subt_size = rel(0.9))


biofinder_cog_comp |>
  ggplot(aes(memory, fill = diagnosis)) +
  geom_density(alpha = 0.3)

biofinder_cog_comp |> drop_na(age) |>
  ggplot(aes(memory, fill = cut(age, 4))) +
  geom_density(alpha = 0.3)

biofinder_cog_comp |>
  ggplot(aes(global_cog, memory, color = diagnosis)) +
  geom_point(alpha = 0.8)

biofinder_cog_comp |>
  ggplot(aes(pathology_ad, memory)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm")

biofinder_cog_comp |> drop_na(age) |>
  ggplot(aes(executive, fill = cut(age, 4))) +
  geom_density(alpha = 0.3)

biofinder_cog_comp |>
  ggplot(aes(age, executive)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm")

biofinder_cog_comp |>
  ggplot(aes(pathology_ad, memory)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm")

biofinder_cog_comp |>
  ggplot(aes(global_cog, executive, color = diagnosis)) +
  geom_point(alpha = 0.8)


cognitive_composites <- plot_gradient_relationships(biofinder_cog_comp |> 
                                                      filter(fmri_bl), 
                                                    gradient_data = grad_df |> filter(study== "biofinder"), 
                                                    gradients = c(1, 3),
                                                    gradient_colors = gradient_cols,
                                                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                                    #atlas_geometry = readRDS(file.path(atlas_dir, "Schaefer2018_1000Parcels_geometry.rds")),
                                                    add_shade = TRUE, 
                                                    shade_alpha = 0.0075,
                                                    shade_size = 0.001,
                                                    l_width = 0.05,
                                                    point_size = 0.75,
                                                    empty_row_height = -0.225,
                                                    base_size_ = 16,
                                                    brain_title_vjust = -0.5,
                                                    r_spin_size = 1,
                                                    vect = TRUE,
                                                    mod_formula = formula(paste0(" ~ age + pathology_ad + global_cog + executive + memory + sex + rsqa__MeanFD")),
                                                    covariates = c("age", "pathology_ad", "sex", "rsqa__MeanFD"),
                                                    plt_title = "Higher executive or memory function show positive gradient alignment",
                                                    brain_names = c("Global Cognition", "Executive", "Memory"),
                                                    plt_subtitle = TRUE,
                                                    group_n_title = FALSE,
                                                    rectangle = TRUE)

cognitive_composites$plot <- cognitive_composites$plot + plot_annotation(theme = theme(plot.subtitle = element_text(vjust = -6)))
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

# pw_mediation <- plot_mediation_gradients(subject_data = df, 
#                                       treatments = c("age", "pathology_ad"), 
#                                       outcomes = c("global_cog"),
#                                       covariates = list(c("pathology_ad", "sex", "rsqa__MeanFD"), 
#                                                         c("age", "sex", "rsqa__MeanFD")),
#                                       gradient_data = grad_df |> filter(study=="biofinder"), 
#                                       gradients = c(1, 3),
#                                       parcels = parcs,
#                                       gradient_colors = gradient_cols,
#                                       plt_title = "Parcelwise Mediation Analysis",
#                                       terms_nice = c("age", "AD path"),
#                                       plot_spacing = 0.3,
#                                       empty_row_height = -0.1,
#                                       base_size_ = 16,
#                                       add_shade = TRUE,
#                                       shade_size = 0.01, 
#                                       shade_alpha = 0.01,
#                                       rectangle = TRUE)

# pw_mediation <- plot_mediation_gradients(subject_data = df |> 
#                                            mutate(age_ = age,
#                                                   pathology_ad_ = pathology_ad), 
#                                          treatments = c("age", "age_", "pathology_ad", "pathology_ad_"), 
#                                          outcomes = c("executive", "memory", "executive", "memory"),
#                                          covariates = list(c("pathology_ad", "sex", "rsqa__MeanFD"), 
#                                                            c("pathology_ad", "sex", "rsqa__MeanFD"),
#                                                            c("age", "sex", "rsqa__MeanFD"), 
#                                                            c("age", "sex", "rsqa__MeanFD")),
#                                          gradient_data = grad_df |> filter(study=="biofinder"), 
#                                          gradients = c(1, 3),
#                                          parcels = parcs,
#                                          gradient_colors = gradient_cols,
#                                          plt_title = waiver(),
#                                          terms_nice = c("Age", "Age", "ADpath", "ADpath"),
#                                          point_size = 0.75,
#                                          plot_spacing = 0.3,
#                                          empty_row_height = -0.15,
#                                          padding_med_str = 0.5,
#                                          l_width = 0.05,
#                                          med_lab_size = 5,
#                                          base_size_ = 16,
#                                          n_lab_size = 11,
#                                          add_shade = TRUE,
#                                          shade_alpha = 0.0075,
#                                          shade_size = 0.001,
#                                          rectangle = TRUE)

pw_mediation_bar <- plot_mediation_barplot(subject_data = df |> 
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
                                           plt_title = "Parcelwise mediation",
                                           terms_nice = c("Age", "Age", "ADpath", "ADpath"),
                                           point_size = 0.75,
                                           plot_spacing = 0.3,
                                           empty_row_height = -0.15,
                                           padding_med_str = 0.5,
                                           l_width = 0.05,
                                           med_lab_size = 5,
                                           base_size_ = 16,
                                           n_lab_size = 11,
                                           add_shade = TRUE,
                                           shade_alpha = 0.0075,
                                           shade_size = 0.001,
                                           rectangle = TRUE)

library(psych)
library(ggraph)
library(tidygraph)

sl <- cog_fac$schmid$sl
sl <- sl[, 1:5]

latent_factors <- colnames(renamed_omega) #c("g", "Memory", "Visuosp", "Language", "Executive")

latent_factors <- str_to_title(str_replace(latent_factors, "_", " "))
domains <- latent_factors[-1]
colnames(sl) <- latent_factors

sl <- as.data.frame(sl) |>
  rownames_to_column("test") |>
  mutate(
    max_domain = domains[ max.col(across(all_of(domains)), ties.method = "first") ],
    max_loading = pmax(!!!rlang::syms(domains))
  ) |>
  filter(max_loading > 0.2) |>
  mutate(
    max_domain = factor(max_domain, levels = latent_factors)
  ) |>
  arrange(max_domain, desc(max_loading)) |> 
  select(-max_domain, -max_loading) |> 
  column_to_rownames("test")


#sl[abs(sl)<0.2] <- NA

item_names <- c(adas_delayed_word_recall = "ADAS Delayed Recall",
                adas_immediate_word_recall_average = "ADAS Immediate Recall",
                bnt_15_2 = "Boston Naming",
                vosp_cube = "VOSP Cube",
                animal_fluency = "Animal Fluency",
                letter_s = "Letter Fluency",
                trailmaking_a = "Trail Making A",
                trailmaking_b = "Trail Making B",
                symbol_digit = "Symbol Digit",
                aqt_color_form = "AQT Color–Form")

item_names <- item_names[rownames(sl)]

sl <- sl[enframe(item_names)$name, ]
rownames(sl) <- item_names



df <- as_tibble(sl, rownames = NA) |> rownames_to_column("item") 

edges <- df %>%
  pivot_longer(cols = -item, names_to = "factor", values_to = "loading") %>%
  filter(!is.na(loading), loading > 0.2) %>%
  mutate(loading = round(loading, 1)) |> 
  rename(from = factor, to = item)

nodes <- tibble(
  name = c(latent_factors, df$item),
  type = c(
    rep("latent", length(latent_factors)),
    rep("item", nrow(df))
  ),
  x = c(
    0, rep(1, 4),   # left side (latent variables)
    rep(0.5, nrow(df))                  # middle (items)
  ),
  y = c(0.5,
        seq(0.8, 0.2, length.out = 4),  # spread latent factors vertically
        seq(1, 0, length.out = nrow(df))                 # spread items vertically
  )
)

graph <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE)


make_push_vector <- function(edges, push_df, default = 0, unit = "mm") {
  merged <- edges %>%
    left_join(push_df, by = c("to", "from")) %>%
    mutate(push = ifelse(is.na(push), default, push))
  unit_list <- lapply(merged$push, function(x) unit(x, unit))
  push_vec <- do.call(grid::unit.c, unit_list)
  return(push_vec)
}

push_vec <- make_push_vector(edges = edges, push_df = tibble(
  to = c("VOSP Cube", "ADAS Delayed Recall", "ADAS Immediate Recall"),
  from = c("Executive", "Memory", "Memory"),
  push = c(15, 5, 5)
))

omega_plot <- ggraph(graph, layout = "manual", x = x, y = y) +
  geom_edge_link(
    aes(label = loading, 
        #start_cap = circle(radius = 1.075),#label_rect(node1.name),
        end_cap = label_rect(node2.name, padding = margin(3, 3, 3, 3, "mm"))),
    angle_calc = 'along',
    arrow = arrow(angle = 15, length = unit(3, "mm"), type = "closed"),
    label_dodge = unit(2, 'mm'),
    #check_overlap = TRUE,
    label_push = push_vec,
    color = "grey30"
  ) +
  geom_node_label(aes(label = name, filter = type == "item"), 
                  label.padding = unit(2, "mm"),
                  label.r = unit(3, "mm"), size = 5) +
  geom_node_circle(aes(filter = type == "latent", r = 0.095), fill = "white") +
  geom_node_text(aes(label = name, filter = type == "latent"), size = 5) +
  theme_void(base_size = 16) +
  ggtitle(" Bifactor model of cognition (Schmid–Leiman)") +
  theme(plot.background = element_rect(color = "black", linewidth = 1))


cognition_full_new <- 
  ggdraw() +
  draw_plot(omega_plot, x = 0, y = 0.505, width = 0.40, height = 0.495) +
  draw_plot(cognitive_composites$plot + plot_annotation(theme = theme(plot.title = element_text(margin = margin(b = -8)))),
            x = 0.41, y = 0.505, width = 0.59, height = 0.495) +
  draw_plot(cog_diff #+ plot_annotation(theme = theme(plot.margin = margin(10, 5, 0, 5)))
            , x = 0.0, y = 0.0, width = 0.35, height = 0.495) +
  draw_plot(pw_mediation_bar, x = 0.36, y = 0, width = 0.64, height = 0.495) + 
  draw_plot_label("A", x = 0.382, y = 1, size = 18, color = "#4e4e4e") +
  draw_plot_label("B", x = 0.98, y = 1, size = 18, color = "#4e4e4e") +
  draw_plot_label("C", x = 0.33, y = 0.495, size = 18, color = "#4e4e4e") +
  draw_plot_label("D", x = 0.98, y = 0.495, size = 18, color = "#4e4e4e") +
  #draw_text("Parcelwise mediation", x = 0.365, y = 0.495, hjust = 0, vjust = 1.5, size = 20) +
  #draw_text("Covariates: \nAge/ADpath \nSex \nMotion", x = 0.3675, y = 0.485, hjust = 0, vjust = 1.5, size = 12)
  draw_text("Covariates: \nAge/ADpath \nSex \nMotion", x = 0.3675, y = 0.445, hjust = 0, vjust = 1.5, size = 11)

scaling_factor <- 3
img_width <-  180 / 25.4
magick_geom_scaling <- paste0(100/scaling_factor, "%x", 100/scaling_factor, "%")
figure_path <- "paper/revision_figures/imputed_cog"
p_name <- "cognition_full.png"
ggsave(file.path(figure_path, p_name), cognition_full_new,
       width = img_width*scaling_factor, height = img_width*0.8*scaling_factor, 
       units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)

############
# Subtypes
############



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

df_cog$cog_subtype |> table()
df_cog$min_domain |> table()

source("src/util.R")
source("src/util_vis.R")

gam_preds_wo_amnestic <- gam_pred_nodes(df_cog |> filter(is.na(cog_subtype) | cog_subtype == "non-amnestic" ) |> 
                                          mutate(pathology_ad == tau_pathology),
                                        fc_matrix = fc_measures_bf$affinity,
                                        roi_names = rois, 
                                        id_var = "image_file",
                                        print_ticker = TRUE,
                                        model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))

point_effects <- do.call(rbind, gam_preds_wo_amnestic$pat_derivs)
point_effects_pat <- point_effects[nrow(point_effects), ]
point_effects <- point_effects[-nrow(point_effects), ]

g_cors <- apply(point_effects, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient1)))

plot(point_effects_pat, g_cors, type = "l")


gam_preds_only_amnestic <- gam_pred_nodes(df_cog |>  filter(is.na(cog_subtype) | cog_subtype == "amnestic" )|> 
                                            mutate(pathology_ad == tau_pathology),
                                          fc_matrix = fc_measures_bf$affinity,
                                          roi_names = rois, 
                                          id_var = "image_file",
                                          print_ticker = TRUE,
                                          model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))

point_effects_amn <- do.call(rbind, gam_preds_only_amnestic$pat_derivs)
point_effects_pat_amn <- point_effects_amn[nrow(point_effects_amn), ]
point_effects_amn <- point_effects_amn[-nrow(point_effects_amn), ]

g_cors_amn <- apply(point_effects_amn, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient1)))

plot(point_effects_pat_amn, g_cors_amn, type = "l")


g3_cors <- c(apply(point_effects, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient3))),
             apply(point_effects_amn, 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient3))))


n_nonamn <- df_cog |> filter(is.na(cog_subtype) | cog_subtype == "non-amnestic" ) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
n_amn <- df_cog |>  filter(is.na(cog_subtype) | cog_subtype == "amnestic" )|> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()

cog_sub_gam <- tibble(g1 = c(g_cors_amn, g_cors), 
                      #g3 = g3_cors,
                      pathology = rep(point_effects_pat_amn, 2),
                      sample = factor(rep(c("Amnestic" |> paste0(" (n = ", n_amn, ")"), "Non-amnestic" |> paste0(" (n = ", n_nonamn, ")")
                      ), each = length(g_cors)),
                      levels = c(c("Amnestic" |> paste0(" (n = ", n_amn, ")"), "Non-amnestic" |> paste0(" (n = ", n_nonamn, ")")))
                      ))   |> 
  #pivot_longer(starts_with("g"), values_to = "g_corr", names_to = "gradient") |> 
  ggplot(aes(x = pathology, y = g1, color = sample)) +
  geom_line(show.legend = F, linewidth = 1) +
  ggsci::scale_color_nejm() +
  facet_wrap(~sample) +
  labs(y = "Gradient Correlation (r)", x = "AD Pathology") +
  ggtitle("Full control sample, subsets of patients")


###############3
# Interaction 
###############




int_plot <- plot_gradient_relationships(biofinder_df %>% filter(fmri_bl), 
                                        gradient_data = grad_df |> filter(study== "biofinder"), 
                                        gradients = c(1, 3),
                                        gradient_colors = gradient_cols,
                                        list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                                        add_shade = TRUE, 
                                        shade_alpha = 0.0075,
                                        shade_size = 0.001,
                                        empty_row_height = -0.15,
                                        base_size_ = 12,
                                        vect = TRUE,
                                        mod_formula = formula(" ~ scale(age) * pathology_ad + sex + rsqa__MeanFD"),
                                        covariates = c("sex", "rsqa__MeanFD"),
                                        brain_names = c("Age (scaled)", "AD Pathology", "Age × ADpath"),
                                        rectangle = TRUE,
                                        plt_title = "Full sample",
                                        plt_subtitle = TRUE)


# variant_plot_agexpath <- 
#   ggdraw() +
#   draw_plot(domain_scores, x = 0, y = 0.505, width = 0.35, height = 0.495) +
#   draw_plot(int_plot$plot, x = 0.36, y = 0.505, width = 0.64, height = 0.495) +
#   draw_plot(cluster_cog$plot, x = 0.0, y = 0, width = 1.0 , height = 0.495) +
#   draw_plot_label("A", x = 0.33, size = 22) +
#   draw_plot_label("B", x = 0.98, size = 22) +
#   draw_plot_label("C", x = 0.98, y = 0.49, size = 22) 


scaling_factor <- 3
p_name <- "interaction_plot.png"
ggsave(file.path(figure_path, p_name), int_plot$plot,
       width = img_width*scaling_factor, height = img_width*0.9*scaling_factor, 
       units = "in", dpi = 300, device = "png", bg = "white")
img <- magick::image_read(file.path(figure_path, p_name))
img_resized <- magick::image_resize(img, magick_geom_scaling)
magick::image_write(img_resized, file.path(figure_path, p_name), density = 300)





eoad_gam <- gam_pred_nodes(biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age>=65)),
                           fc_matrix = fc_measures_bf$affinity,
                           roi_names = rois, 
                           id_var = "image_file",
                           print_ticker = TRUE,
                           model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))


eff_eoad_gam <- apply(do.call(rbind, eoad_gam$pat_derivs)[-1001, ], 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient1)))
plot(point_effects_pat, eff_eoad_gam, type = "l")


load_gam  <- gam_pred_nodes(biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age<=65)),
                            fc_matrix = fc_measures_bf$affinity,
                            roi_names = rois, 
                            id_var = "image_file",
                            print_ticker = TRUE,
                            model_formula = formula("FC ~ s(age) + s(pathology_ad, k = 5) + sex + rsqa__MeanFD"))


eff_load_gam <- apply(do.call(rbind, load_gam$pat_derivs)[-1001, ], 2, function(x) cor(x, grad_df |> filter(study == "biofinder") |> pull(gradient1)))
plot(point_effects_pat, eff_load_gam, type = "l")

n_eoad <- biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age>=65)) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()
n_load <- biofinder_df %>% filter(fmri_bl) |> filter(!(diagnosis %in% c("MCI", "AD") & age<=65)) |> drop_na(age, pathology_ad, sex, rsqa__MeanFD) |> count()

age_gams <- tibble(g1 = c(eff_load_gam , eff_eoad_gam),
                   #g3 = g3_cors,
                   pathology = rep(point_effects_pat_amn, 2),
                   sample = factor(rep(c("MCI/AD >65" |> paste0(" (n = ", n_load, ")"), "MCI/AD <65"  |> paste0(" (n = ", n_eoad, ")")
                   ), each = length(g_cors)),
                   levels = c(c("MCI/AD >65" |> paste0(" (n = ", n_load, ")"), "MCI/AD <65"  |> paste0(" (n = ", n_eoad, ")")
                   ))
                   )) |> 
  ggplot(aes(x = pathology, y = g1, #color = sample
  )) +
  geom_line(show.legend = F, linewidth = 1) +
  ggsci::scale_color_nejm() +
  #ylim(-0.1, 0.75) +
  theme_bw(base_size = 12) +
  facet_wrap(~sample) +
  #ylim(-0.1, 0.7) +
  labs(y = "Gradient Correlation (r)", x = "AD Pathology") +
  ggtitle("Full control sample, subsets of patients") +
  theme(#plot.background = element_rect(color = "black"),
    plot.title = element_text(size = rel(1)
    ))

cog_sub_gam <- tibble(g1 = c(g_cors_amn, g_cors), 
                      #g3 = g3_cors,
                      pathology = rep(point_effects_pat, 2),
                      sample = factor(rep(c("Amnestic" |> paste0(" (n = ", n_amn, ")"), "Non-amnestic" |> paste0(" (n = ", n_nonamn, ")")
                      ), each = length(g_cors)),
                      levels = c(c("Amnestic" |> paste0(" (n = ", n_amn, ")"), "Non-amnestic" |> paste0(" (n = ", n_nonamn, ")")))
                      ))   |> 
  #pivot_longer(starts_with("g"), values_to = "g_corr", names_to = "gradient") |> 
  ggplot(aes(x = pathology, y = g1, #color = sample
  )) +
  geom_line(show.legend = F, linewidth = 1) +
  ggsci::scale_color_nejm() +
  #ylim(-0.1, 0.75) +
  facet_wrap(~sample) +
  theme_bw(base_size = 12) +
  labs(y = "Gradient Correlation (r)", x = "AD Pathology") +
  #ggtitle("Full control sample, subsets of patients") +
  theme(#plot.background = element_rect(color = "black"),
    plot.title = element_text(size = rel(1)
    ))


pal <- ggsci::pal_futurama()(3)

biofinder_df |> 
  filter(fmri_bl, !is.na(abnorm_ab)) |> 
  mutate(`Aβ+` = as.logical(abnorm_ab),
         Diagnosis = case_when(
           diagnosis %in% c("Normal", "SCD") ~ "CU",
           TRUE ~ diagnosis
         ) |> factor(levels = c("CU", "MCI", "AD"))) |> 
  ggplot(aes(pathology_ad, age, color = Diagnosis#`Aβ+`
  )) +
  geom_point(alpha=0.2, size = 0.75) +
  geom_smooth(data = biofinder_df |> 
                filter(fmri_bl, !is.na(abnorm_ab), pathology_ad>0.5), 
              color = "black",
              method = "lm") +
  geom_smooth(data = biofinder_df |> 
                filter(fmri_bl, !is.na(abnorm_ab), pathology_ad<0.5), 
              color = "black",
              method = "lm") +
  labs(x = "Pathology Score", y = "Age") +
  scale_color_manual(
    values = c(pal[3], pal[1], pal[2]),
    labels = c("CU", "MCI (Aβ+)", "AD (Aβ+)")
  )+
  theme_bw(base_size = 12) +
  theme(legend.position = "top") +
  theme(plot.background = element_rect(color = "black")) -> age_path

gam_plots <- 
  wrap_plots(list(age_gams, cog_sub_gam), ncol = 1, guides = "collect", axis_titles = "collect", axes = "collect") +
  plot_annotation(theme = theme(plot.background = element_rect(color = "black")))



age_plots <- ggdraw()+ 
  draw_plot(int_plot$plot, x = 0, y = 0, width = 0.60,  height = 1) + 
  draw_plot(gam_plots, x = 0.61, y = 0.0, width = 0.39, height = 1) 

#draw_plot(cog_sub_gam, x = 0.61, y = 0, width = 0.39, height = 0.495) 


scaling_factor <- 1
img_width <- 12
p_name <- "age_plot.png"
ggsave(file.path(figure_path, p_name), age_plots,
       width = img_width*scaling_factor, height = img_width*0.5*scaling_factor, 
       units = "in", dpi = 300, device = "png", bg = "white")
