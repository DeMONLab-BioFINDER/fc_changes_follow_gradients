fc_strength <- function(connectomes, atlas = yeo_msk, roi_names, threshold = 0, replace = 0){
  
  zero_out_mat <- function(mat, thresh, replacement) {
    library(matrixStats)
    thresholds <- colQuantiles(abs(mat), probs = thresh)
    mat[abs(mat) <= thresholds[col(mat)]] <- replacement
    return(mat)
  }
  
  sub_ids <- dimnames(connectomes)[[3]]
  
  n_subs <- dim(connectomes)[3]
  n_nodes <- dim(connectomes)[1]
  
  globals <- matrix(0, nrow = n_subs, ncol = 4, dimnames = list(sub_ids, c("strength", "within", "between", "segregation")))
  node_strength <- matrix(0, nrow = n_subs, ncol = n_nodes, dimnames = list(sub_ids, roi_names))
  node_strength_within <- matrix(0, nrow = n_subs, ncol = n_nodes, dimnames = list(sub_ids, roi_names))
  node_strength_between <- matrix(0, nrow = n_subs, ncol = n_nodes, dimnames = list(sub_ids, roi_names))
  nodal_sys_seg <- matrix(0, nrow = n_subs, ncol = n_nodes, dimnames = list(sub_ids, roi_names))
  
  msk <- atlas
  network_same_mat <- outer(msk, msk, FUN = "==")
  diag(network_same_mat) <- FALSE 
  network_diff_mat <- !network_same_mat
  diag(network_diff_mat) <- FALSE 
  
  within_lower <- network_same_mat
  within_lower[upper.tri(within_lower)] <- FALSE
  
  between_lower <- network_diff_mat
  between_lower[upper.tri(between_lower)] <- FALSE
  
  
  network_same_counts <- rowSums(network_same_mat)
  network_diff_counts <- rowSums(network_diff_mat)
  
  
  pb = txtProgressBar(min = 0, max = length(sub_ids), initial = 0) 
  for (sub in 1:dim(connectomes)[3]){
    
    conn <- connectomes[, , sub]
    if (threshold > 0) {
      conn <- zero_out_mat(conn, thresh = threshold, replacement = replace) 
    }
    conn[conn < 0] <- replace
    diag(conn) <- NA
    
    
    glob_str <- mean(conn[upper.tri(conn)], na.rm = TRUE)
    glob_within <- mean(conn[within_lower], na.rm = TRUE)
    glob_between <- mean(conn[between_lower], na.rm = TRUE)
    glob_seg <- (glob_within - glob_between) / glob_within
    globals[sub, ] <- c(glob_str, glob_within, glob_between, glob_seg)
    
    node_str <- rowMeans(conn, na.rm = TRUE)
    within_fc <- conn
    within_fc[network_diff_mat] <- NA
    node_within <- rowMeans(within_fc, na.rm = TRUE)
    
    
    between_fc <- conn
    between_fc[network_same_mat] <- NA
    node_between <- rowMeans(between_fc, na.rm = TRUE)
    
    node_seg <- (node_within - node_between) / node_within
    
    node_strength[sub, ] <- node_str
    node_strength_within[sub, ] <- node_within
    node_strength_between[sub, ] <- node_between
    nodal_sys_seg[sub, ] <- node_seg
    
    setTxtProgressBar(pb, sub)
  }
  
  close(pb)
  return(list(strength = node_strength, within = node_strength_within, 
              between = node_strength_between, segregation = nodal_sys_seg, globals = globals))
}

zero_out_mat <- function(mat, thresh) {
  library(matrixStats)
  thresholds <- colQuantiles(abs(mat), probs = thresh)
  mat[abs(mat) <= thresholds[col(mat)]] <- 0
  return(mat)
}


get_affinity <- function(con_cube, similarity_method="correlation",
                         threshold = 0, 
                         roi_names = rois,
                         atlas = yeo_msk) {
  library(proxyC)
  library(matrixStats)
  
  
  zero_out_mat <- function(mat, thresh) {
    thresholds <- colQuantiles(abs(mat), probs = thresh)
    mat[abs(mat) <= thresholds[col(mat)]] <- 0
    return(mat)
  }
  
  n_sub <- dim(con_cube)[3]
  n_nodes <- dim(con_cube)[2]
  sub_names <- dimnames(con_cube)[[3]]
  
  mat_similarity <- matrix(0, nrow = n_sub, ncol = n_nodes, dimnames = list(sub_names, roi_names))
  within_affinity <- matrix(0, nrow = n_sub, ncol = n_nodes, dimnames = list(sub_names, roi_names))
  between_affinity <- matrix(0, nrow = n_sub, ncol = n_nodes, dimnames = list(sub_names, roi_names))
  
  msk <- atlas
  network_same_mat <- outer(msk, msk, FUN = "==")
  diag(network_same_mat) <- FALSE 
  network_diff_mat <- !network_same_mat
  diag(network_diff_mat) <- FALSE 
  
  
  pb = txtProgressBar(min = 0, max = n_sub, initial = 0) 
  for(i in 1:n_sub) {
    connectome <- con_cube[, , i]
    L <- matrix(connectome[!diag(nrow(connectome))], nrow = nrow(connectome)-1) # mask it
    L <- zero_out_mat(L, threshold)
    L[L<0] <- 0
    if (similarity_method == "correlation") {
      L <- cor(L)
      diag(L) <- NA
      
    } else {
      L <- Matrix::Matrix(L, sparse = TRUE)
      L <- proxyC::simil(L, margin = 2, method = similarity_method#, use_nan = TRUE
      )
      L <- as.matrix(L)
      diag(L) <- NA
    }
    
    within_aff <- L
    within_aff[network_diff_mat] <- NA
    node_within <- rowMeans(within_aff, na.rm = TRUE)
    within_affinity[i, ] <- node_within
    
    between_aff <- L
    between_aff[network_same_mat] <- NA
    node_between <- rowMeans(between_aff, na.rm = TRUE)
    between_affinity[i, ] <- node_between
    
    
    mean_similarity_for_sub <- colMeans(L , na.rm = TRUE)
    mat_similarity[i, ] <- mean_similarity_for_sub
    setTxtProgressBar(pb,i)
  }
  close(pb)
  
  return(list(affinity = mat_similarity, affinity_within = within_affinity, affinity_between = between_affinity))
  
}




nodal_regression_fits <- function(subject_data,
                                  inter_sub_sim_matrix, 
                                  roi_names,
                                  logistic = FALSE,
                                  scale_fc = FALSE,
                                  vectorised = FALSE,
                                  id_var = "image_file",
                                  dep_var = "FC",
                                  model_formula = formula(" ~ age + pathology_ad + sex + rsqa__MeanFD")) {
  
  
  
  if (vectorised) {
    
    bio_filtered <- subject_data %>% 
      drop_na(all_of(all.vars(model_formula)))
    X <- model.matrix(model_formula, data = bio_filtered)
    Y <- as.data.frame(inter_sub_sim_matrix[bio_filtered[[id_var]], ])
    
    similarity_fits <- list()
    for (roi in roi_names) {
      fit <- lm.fit(X, Y[, roi])
      fit$formula <- model_formula
      similarity_fits[[roi]] <- fit
    }
    
    return(similarity_fits)
    
  } else {
    
    inter_sub_long <- as.data.frame(inter_sub_sim_matrix) %>% rownames_to_column(id_var) %>%
      inner_join(subject_data, by = id_var) %>%
      pivot_longer(starts_with('7Networks'), names_to = 'region', values_to = dep_var)
    similarity_fits <- list()
    ticker = 0
    pb = txtProgressBar(min = 0, max = length(roi_names), initial = 0) 
    for (roi in roi_names) {
      if (scale_fc) {
        reg_dat <- inter_sub_long %>% filter(region == roi) %>% mutate(dep_var = scale(.data[[dep_var]]))  
      } else {
        reg_dat <- inter_sub_long %>% filter(region == roi)  
      }
      
      if (!logistic) {
        fit <- lm(model_formula, data = reg_dat)
      }
      if (logistic) {
        fit <- glm(model_formula, data = reg_dat, family = binomial(link = "logit"))
      }
      similarity_fits[[roi]] <- fit

        ticker <- ticker + 1
        setTxtProgressBar(pb,ticker)

    }
    close(pb)
    return(similarity_fits)
  } 
    
}

get_nodal_ests <- function(fit_list, roi_names = rois, vectorised = TRUE, mc = FALSE){
  
 
  tictoc::tic()
  
      if(vectorised) {
        
        if (!mc) {
          fits_ests <- data.frame()
          pb = txtProgressBar(min = 0, max = length(roi_names), initial = 0, style = 3) 
          i = 0
          for (roi in roi_names) {
            z <- fit_list[[roi]]
            
            p <- z$rank
            rdf <- z$df.residual
            r <- z$residuals
            f <- z$fitted.values
            rss <- sum(r^2)
            resvar <- rss/rdf
            p1 <- 1L:p
            R <- chol2inv(z$qr$qr[1:p, 1:p, drop=FALSE])
            se <- sqrt(diag(R) * resvar)
            est <- z$coefficients[z$qr$pivot[p1]]
            tval <- est/se
            
            t_values <- tval %>% enframe("term", "statistic") %>%
              mutate(region = roi,
                     n = length(r),
                     model_formula = deparse1(z$formula))
            fits_ests <- rbind(fits_ests, t_values)
            setTxtProgressBar(pb, i)
            i = i + 1
          }
          close(pb)
          
        } else if(mc) {
          library(parallel)
          library(pbmcapply)
          chunk_size <- 20
          
          roi_chunks <- split(roi_names, ceiling(seq_along(roi_names) / chunk_size))
          
          n_cores <- detectCores()-1
          
          res_chunks <- pbmclapply(roi_chunks, FUN = function(rois_subset) {
            
            chunk_results <- lapply(rois_subset, function(roi) {
              z <- fit_list[[roi]]
              
              p <- z$rank
              rdf <- z$df.residual
              r <- z$residuals
              f <- z$fitted.values
              rss <- sum(r^2)
              resvar <- rss/rdf
              p1 <- 1L:p
              R <- chol2inv(z$qr$qr[1:p, 1:p, drop=FALSE])
              se <- sqrt(diag(R) * resvar)
              est <- z$coefficients[z$qr$pivot[p1]]
              tval <- est/se
              
              df_tvals <- tibble(
                term          = names(tval),
                statistic     = tval,
                region        = roi,
                n             = length(r),
                model_formula = deparse1(deparse1(z$formula))
              )
              
              df_tvals
            })
            
            do.call(rbind, chunk_results)
            
          }, mc.cores = n_cores)
          
          fits_ests <- do.call(rbind, res_chunks)
        }

      } else {
        library(broom)
        fits_ests <- data.frame()
        for (roi in roi_names) {
          fit <- fit_list[[roi]]
          coeffs <- tidy(fit) 
          coeffs$region <- roi
          coeffs$n <- nobs(fit)
          coeffs$model_formula <-  deparse1(formula(fit))
          fits_ests <- rbind(fits_ests, coeffs)
      }
      }
  tictoc::toc()
  return(fits_ests)
}



nodal_lmm_ests <- function(subject_data,
                           inter_sub_sim_matrix, 
                           roi_names,
                           mc = TRUE,
                           id_var = "image_file",
                           subject_id = "sid",
                           dep_var = "FC",
                           model_formula = formula("FC ~ age + tau_pathology + sex + rsqa__MeanFD + (1 | sid)")){
  
  
  if(mc) {
    
    library(parallel)
    library(pbmcapply)
    library(lme4)
    
    chunk_size <- 20
    
    roi_chunks <- split(roi_names, ceiling(seq_along(roi_names) / chunk_size))
    
    n_cores <- detectCores()-1
    
    fixedfx_chunks <- pbmclapply(roi_chunks, FUN = function(chunk) {
      
      chunk_list <- lapply(chunk, function(roi) {
        
        roi_fc <- inter_sub_sim_matrix[, roi] %>% 
          enframe(id_var, "FC")
        
        reg_df <- subject_data %>% 
          inner_join(roi_fc, by = id_var)
        
        fit <- lmer(model_formula, data = reg_df)
        
        tvals <- coef(summary(fit))[, "t value"]
        
        x <- tibble(
          term          = rownames(coef(summary(fit))),
          statistic     = tvals,
          region        = roi,
          n             = ngrps(fit),
          model_formula = deparse1(formula(fit))
        )
        
      })
      
      
      do.call(rbind, chunk_list)
      
    }, mc.cores = n_cores)
    
    fixedfx <- do.call(rbind, fixedfx_chunks)
    return(fixedfx)
    
  } else {
    library(lme4)
    longitudinal_subject_data <- subject_data 
    fixedfx <- data.frame()
    pb = txtProgressBar(min = 0, max = length(rois), initial = 0)
    i=0
    for(roi in rois)  {
      roi_fc <- inter_sub_sim_matrix[, roi] %>% enframe(id_var, "FC")
      
      reg_df <- longitudinal_subject_data %>% inner_join(roi_fc, by = id_var)
      
      fit <- lmer(model_formula, data = reg_df)
      
      fixedfx <- rbind(fixedfx, summary(fit)[["coefficients"]][, "t value"] %>% 
                         enframe("term", "statistic") %>%
                         mutate(region = roi, n = ngrps(fit), model_formula = deparse1(formula(fit))))
      i = i +1
      setTxtProgressBar(pb,i)
    }
    close(pb)
    return(fixedfx)
  }

}

seq_range <- function(x, n) {
  seq(min(x), max(x), length.out = n)
}


gam_pred_nodes <- function(subject_data,
                           fc_matrix, 
                           roi_names,
                           print_ticker = TRUE, 
                           id_var = "image_file",
                           dep_var = "FC",
                           model_formula = formula("FC ~ s(age) + s(pathology_ad) + sex")) {
  library(mgcv)
  library(gratia)
  fc_matrix_long <- as.data.frame(fc_matrix) %>% rownames_to_column(id_var) %>% 
    inner_join(subject_data, by = id_var) %>% 
    pivot_longer(starts_with('7Networks'), names_to = 'region', values_to = dep_var)
  
  
  new_age <- subject_data %>% pull(age) %>% quantile(probs = c(0.05, 0.95), na.rm = TRUE) %>% seq_range(100)
  #similarity_fits <- list()
  pat_predicted <- c()
  age_predicted <- c()
  pat_derivs <- c()
  age_derivs <- c()
  ticker = 0
  if (print_ticker) pb = txtProgressBar(min = 0, max = length(roi_names), initial = 0) 
  for (roi in roi_names) {
    gam_fit <- gam(model_formula, data = fc_matrix_long %>% filter(region == roi))
    #similarity_fits[[roi]] <- gam_fit
    
    pat_newdata <-  expand_grid(age = 66.65896, pathology_ad = seq(0, 1, length.out = 100), rsqa__MeanFD = 0.17, sex = c(0, 1))
    pat_newdata$pred <- predict(gam_fit , newdata = pat_newdata)
    pat_pred <- pat_newdata %>% group_by(pathology_ad) %>% summarise(predicted_fc = mean(pred))
    
    age_newdata <-  expand_grid(age = new_age, pathology_ad = 0.1, rsqa__MeanFD = 0.17, sex = c(0, 1))
    age_newdata$pred <- predict(gam_fit , newdata = age_newdata)
    age_pred <- age_newdata %>% group_by(age) %>% summarise(predicted_fc = mean(pred))
    
    pat_deriv <- derivatives(gam_fit , type = "central" , select = "pathology_ad", partial_match = TRUE, data = pat_newdata)
    pat_deriv <- pat_deriv %>% group_by(pathology_ad) %>% summarise(averaged_derivative = mean(.derivative))
    age_deriv <- derivatives(gam_fit , type = "central" , select = "age", partial_match = TRUE, data = age_newdata)
    age_deriv <- age_deriv %>% group_by(age) %>% summarise(averaged_derivative = mean(.derivative))
    pat_derivs <- cbind(pat_derivs, pat_deriv$averaged_derivative)
    age_derivs <- cbind(age_derivs, age_deriv$averaged_derivative)
    
    pat_predicted <- cbind(pat_predicted, pat_pred$predicted_fc)
    age_predicted <- cbind(age_predicted, age_pred$predicted_fc)
    
    if (print_ticker) {
      ticker <- ticker + 1
      setTxtProgressBar(pb,ticker)
    }
  }
  close(pb)
  colnames(pat_predicted) <- roi_names
  colnames(age_predicted) <- roi_names
  pat_preds <- as_tibble(pat_predicted) %>% mutate(pathology_ad = pat_pred$pathology_ad)
  age_preds <- as_tibble(age_predicted) %>% mutate(age = age_pred$age)
  
  
  colnames(pat_derivs) <- roi_names
  colnames(age_derivs) <- roi_names
  pat_derivs <- as_tibble(pat_derivs) %>% mutate(pathology_ad = pat_deriv$pathology_ad)
  age_derivs <- as_tibble(age_derivs) %>% mutate(age = age_deriv$age)
  
  return(list(pat_pred =pat_preds, age_pred = age_preds, pat_derivs =pat_derivs, age_derivs = age_derivs))
}

calc_ests <- function(subject_data,
                      list_of_parcel_data, 
                      id_var = "image_file",
                      mod_formula = formula(paste0("fc_derivative ~ age")),
                      covariates = c("sex", "rsqa__MeanFD"),
                      filter_criteria = quo(),
                      vectorised = TRUE,
                      longitudinal = FALSE,
                      sub_id = "sid",
                      longitudinal_formula = formula(paste0("fc_derivative ~ age + (1 | ", sub_id, ")"))) {
  
  source("src/util.R")
  require(tidyverse)
  
  
  
  analysis_name <- names(list_of_parcel_data)
  
  if (!longitudinal) {
    list_of_fits <- list()
    for (analysis in analysis_name) {
      print(paste0("Running linear models for ", analysis))
      list_of_fits[[analysis]] <- 
        nodal_regression_fits(
          subject_data %>% filter(!!filter_criteria), 
          list_of_parcel_data[[analysis]], 
          roi_names = rois,
          id_var = id_var,
          model_formula = mod_formula,
          vectorised = vectorised)
    }
    
    list_of_ests <- list()
    print("Getting model estimates")
    for (analysis in names(list_of_fits)){
      list_of_ests[[analysis]] <- get_nodal_ests(list_of_fits[[analysis]], vectorised = vectorised, mc = TRUE) %>% 
        select(term, region, statistic)
    }
    
    
  }
  
  if (longitudinal) {
    list_of_ests <- list()
    print("fitting longitudinal")
    for (analysis in analysis_name){
      print(analysis)
      list_of_ests[[analysis]] <- nodal_lmm_ests(
        subject_data %>% filter(!!filter_criteria),
        list_of_parcel_data[[analysis]], 
        roi_names = rois,
        id_var = id_var,
        subject_id = sub_id,
        model_formula = longitudinal_formula
      )
    }
  }
  
  ests <- list_of_ests[[1]] %>% rename_with(~ names(list_of_ests[1]), statistic)
  if (length(list_of_ests) > 1) {
    for(i in 2:length(list_of_ests)) {
      ests <- inner_join(ests, list_of_ests[[i]] %>% rename_with(~ names(list_of_ests[i]), statistic), by = c("term", "region"))
    }
  }
  
  ests
  
}

calc_grad_rel <- function(ests, gradient_data, calc_r2 = FALSE, gradients = c(1, 3),
                          perms = readRDS("data/atlas_data/permutations_1000_hungarian.rds"),
                          covariates = c("sex", "rsqa__MeanFD")) {
  
  grad_char <- paste0("gradient", gradients)
  terms_of_interest <- unique(ests$term)[!(unique(ests$term) %in% c(covariates, "(Intercept)"))]
  res_data <- data.frame()
  
  for (g in grad_char) {
    for (term_ in terms_of_interest){
      
      if (calc_r2) {
        plot_data  <- gradient_data %>% inner_join(ests %>% filter(term == term_), by = "region") 
        mod_form <- formula(paste(g, "~ nodal_affinity"))
        fit <- lm(mod_form, data = plot_data)
        tryCatch( 
          {
            ci <- confintr::ci_rsquared(fit)[["interval"]]
          },
          error = function(cond) {
            ci <- c(NA, NA)
          })
        
        results_row <- data.frame(term = term_, gradient = g, r2 = summary(fit)[["r.squared"]], lo_ci = ci[1],  hi_ci = ci[2])
        res_data <- rbind(res_data, results_row)
      } else {
        plot_data  <- gradient_data %>% inner_join(ests %>% filter(term == term_), by = "region") 
        results_row <- data.frame(term = term_, gradient = g, r = cor(plot_data[[g]], plot_data[["nodal_affinity"]]), 
                                  p_spin = perm_sphere_p(plot_data[[g]], plot_data[["nodal_affinity"]], perm.id = perms, corr.type = "pearson") )
        res_data <- rbind(res_data, results_row)
        
      }

      
    }
  }
  
  res_data
}




write_gradient_supp_table <- function(params){
  cs_bf <- calc_ests(biofinder_df %>% filter(fmri_bl), 
                     list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                     mod_formula = formula(paste0(" ~ age + pathology_ad + sex + rsqa__MeanFD")),
                     covariates = c("sex", "rsqa__MeanFD"))
  
  
  cs_ad <- calc_ests(adni_df %>% filter(fmri_bl), 
                     list_of_parcel_data = list(nodal_affinity = fc_measures_adni$affinity),
                     mod_formula = formula(paste0(" ~ age + pathology_ad + sex + rsqa__MeanFD")),
                     covariates = c("sex", "rsqa__MeanFD"),
                     id_var = "id_ses")
  
  
  health <- calc_ests(biofinder_df %>% filter(fmri_bl, diagnosis=="Normal" | diagnosis=="SCD", abnorm_ab==0, !apoe4),
                      mod_formula = formula(paste0(" ~ scale(age) * scale(-mPACC_v1) + pathology_ad + sex + rsqa__MeanFD")), 
                      list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                      covariates = c("sex", "rsqa__MeanFD"),
                      filter_criteria = quo())
  
  
  clin <- calc_ests(biofinder_df %>% filter(fmri_bl, diagnosis=="MCI" | diagnosis=="AD", !is.na(mPACC_v1)) %>% 
                      mutate(`-mPACC` = -mPACC_v1), 
                    list_of_parcel_data = list(nodal_affinity = fc_measures_bf$affinity),
                    mod_formula = formula(paste0(" ~ age + pathology_ad + `-mPACC` +  sex + rsqa__MeanFD")),
                    covariates = c("sex", "rsqa__MeanFD"),
                    filter_criteria = quo())
  
  
  
  supp_table_grads <- data.frame()
  for (i in 1:nrow(params)) {
    param_i <- params[i, ]
    gradient_data_i <- gradient_over_params %>% filter(study %in% c("biofinder", "adni"), method == param_i$method, affinity == param_i$affinity, threshold == param_i$threshold)
    
    cs_bf_grad_rel <- calc_grad_rel(cs_bf, gradient_data_i %>% filter(study=="biofinder"))
    cs_bf_grad_rel <- cs_bf_grad_rel %>% mutate(
      cohort = "BioFINDER",
      filter = "Full sample", 
      model = "FC ~ age + pathology",
      method = param_i$method,
      threshold = param_i$threshold,
      affinity = param_i$affinity,
      order = 1
    )
    
    
    cs_ad_grad_rel <- calc_grad_rel(cs_ad, gradient_data_i %>% filter(study=="adni"))
    cs_ad_grad_rel <- cs_ad_grad_rel %>% mutate(
      cohort = "ADNI",
      filter = "Full sample", 
      model = "FC ~ age + pathology",
      method = param_i$method,
      threshold = param_i$threshold,
      affinity = param_i$affinity,
      order = 2
    )
    
    health_grad_rel <- calc_grad_rel(health, gradient_data_i %>% filter(study=="biofinder"))
    health_grad_rel <- health_grad_rel %>% mutate(
      cohort = "BioFINDER",
      filter = "cognitively healthy", 
      model = "FC ~ age*(-)mPACC + pathology",
      method = param_i$method,
      threshold = param_i$threshold,
      affinity = param_i$affinity,
      order = 3
    )
    
    clin_grad_rel <- calc_grad_rel(clin, gradient_data_i %>% filter(study=="biofinder"))
    clin_grad_rel <- clin_grad_rel %>% mutate(
      cohort = "BioFINDER",
      filter = "impaired (MCI and AD)", 
      model = "FC ~ age + pathology + -mPACC",
      method = param_i$method,
      threshold = param_i$threshold,
      affinity = param_i$affinity,
      order = 4
    )
    
    all <- rbind(cs_bf_grad_rel, cs_ad_grad_rel, health_grad_rel, clin_grad_rel)
    supp_table_grads <- rbind(supp_table_grads, all)
  }
  supp_table_grads
}


rotate.parcellation = function(coord.l,coord.r,nrot=10000,method='hungarian') {
  
  # Function from Váša F.,. (2017). Cerebral Cortex, 28(1):281–294.
  # and https://github.com/frantisekvasa/rotate_parcellation
  
  require(matrixStats)
  require(clue)
  # check that coordinate dimensions are correct
  if (!all(dim(coord.l)[2]==3,dim(coord.r)[2]==3)) {
    if (all(dim(coord.l)[1]==3,dim(coord.r)[1]==3)) {
      print('transposing coordinates to be of dimension nROI x 3')
      coord.l <- t(coord.l)
      coord.r <- t(coord.r)
    }
  }
  
  nroi.l <- dim(coord.l)[1]   # n(regions) in the left hemisphere
  nroi.r <- dim(coord.r)[1]   # n(regions) in the right hemisphere
  nroi <- nroi.l+nroi.r       # total n(regions)
  
  perm.id <- array(0,dim=c(nroi,nrot)); # initialise output array
  r <- 0; c <- 0; # count successful (r) and unsuccessful (c) iterations
  
  # UPDATED 16/10/2019 - set up updated permutation scheme 
  I1 <- diag(3); I1[1,1] = -1;
  # main loop -  use of "while" is to ensure any rotation that maps to itself is excluded (this is rare, but can happen)
  while (r < nrot) {
    
    # UPDATED 16/10/2019
    A <- matrix(rnorm(9, mean = 0, sd = 1), nrow = 3, ncol = 3)
    qrdec <- qr(A)       # QR decomposition
    TL <- qr.Q(qrdec)    # Q matrix
    temp <- qr.R(qrdec)  # R matrix
    TL <- TL%*%diag(sign(diag(temp)))
    if (det(TL)<0) {
      TL[,1] <- -TL[,1]
    }
    # reflect across the Y-Z plane for right hemisphere
    TR <- I1 %*% TL %*% I1;
    coord.l.rot <- coord.l %*% TL; # transformed (rotated) left coordinates
    coord.r.rot <- coord.r %*% TR; # transformed (rotated) right coordinates
    
    
    
    dist.l <- array(0,dim=c(nroi.l,nroi.l));
    dist.r <- array(0,dim=c(nroi.r,nroi.r));
    
    for (i in 1:nroi.l) { # left
      for (j in 1:nroi.l) {
        dist.l[i,j] <- sqrt( sum( (coord.l[i,]-coord.l.rot[j,])^2 ) )
      }
    }
    for (i in 1:nroi.r) { # right
      for (j in 1:nroi.r) {
        dist.r[i,j] <- sqrt( sum( (coord.r[i,]-coord.r.rot[j,])^2 ) )
      }
    }
    
    # UPDATED 5/5/2023
    if (method == 'vasa') {
      
      # LEFT
      # calculate distances, proceed in order of "most distant minimum"
      # -> for each unrotated region find closest rotated region (minimum), then assign the most distant pair (maximum of the minima), 
      # as this region is the hardest to match and would only become harder as other regions are assigned
      temp.dist.l <- dist.l
      rot.l <- c(); ref.l <- c();
      #tba.r = tba.c = 1:nroi.l # rows and columns that are yet "to be assigned"
      for (i in 1:nroi.l) {
        # max(min) (described above)
        ref.ix <- which( rowMins(temp.dist.l,na.rm=T) == max(rowMins(temp.dist.l,na.rm=T),na.rm=T) )   # "furthest" row
        rot.ix <- which( temp.dist.l[ref.ix,] == min(temp.dist.l[ref.ix,],na.rm=T) ) # closest region
        
        # # alternative option: mean of row - take the closest match for unrotated region that is on average furthest from rotated regions
        # ref.ix = which(nanmean(temp.dist.l,2)==nanmax(nanmean(temp.dist.l,2)))    # "furthest" row
        # rot.ix = which(temp.dist.l(ref.ix,:)==nanmin(temp.dist.l(ref.ix,:)))      # closest region    
        ref.l <- c(ref.l,ref.ix) # store reference and rotated indices
        rot.l <- c(rot.l,rot.ix)
        temp.dist.l[,rot.ix] = NA # set temporary column indices to NaN, to be disregarded in next iteration
        temp.dist.l[ref.ix,] = 0 # because in the above form of the code, R doesn't deal well with whole rows and columns of NaN, set row to low value (which won't matter as furthest rows are assigned first)
        #temp.dist.l[,rot.ix] = NA # set temporary indices to NaN, to be disregarded in next iteration
        #temp.dist.l[ref.ix,] = NA
      }
      
      # RIGHT
      # calculate distances, proceed in order of "most distant minimum"
      # -> for each unrotated region find closest rotated region (minimum), then assign the most distant pair (maximum of the minima), 
      # as this region is the hardest to match and would only become harder as other regions are assigned
      temp.dist.r <- dist.r;
      rot.r <- c(); ref.r <- c();
      for (i in 1:nroi.r) {
        # max(min) (described above)
        ref.ix <- which( rowMins(temp.dist.r,na.rm=T) == max(rowMins(temp.dist.r,na.rm=T),na.rm=T) )   # "furthest" row
        rot.ix <- which( temp.dist.r[ref.ix,] == min(temp.dist.r[ref.ix,],na.rm=T) )             # closest region
        
        # # alternative option: mean of row - take the closest match for unrotated region that is on average furthest from rotated regions
        # ref.ix = which(nanmean(temp.dist.r,2)==nanmax(nanmean(temp.dist.r,2)))    # "furthest" row
        # rot.ix = which(temp.dist.r(ref.ix,:)==nanmin(temp.dist.r(ref.ix,:)))      # closest region
        ref.r = c(ref.r,ref.ix) # store reference and rotated indices
        rot.r = c(rot.r,rot.ix)
        temp.dist.r[,rot.ix] = NA # set temporary column indices to NaN, to be disregarded in next iteration
        temp.dist.r[ref.ix,] = 0 # because in the above form of the code, R doesn't deal well with whole rows and columns of NaN, set row to low value (which won't matter as furthest rows are assigned first)
        #temp.dist.l[,rot.ix] = NA # set temporary indices to NaN, to be disregarded in next iteration
        #temp.dist.l[ref.ix,] = NA
      }
      
      # UPDATED 5/5/2023
    } else if (method == 'hungarian') {
      
      # solve_LSAP(), from library(clue): "Solve the linear sum assignment problem using the Hungarian method."
      
      # LEFT
      rot.l <- as.vector(solve_LSAP(dist.l, maximum=FALSE))
      ref.l <- c(1:nroi.l) # for compatibility with method='vasa'
      
      # RIGHT
      rot.r <- as.vector(solve_LSAP(dist.r, maximum=FALSE))
      ref.r <- c(1:nroi.r) # for compatibility with method='vasa'
      
    } else {
      
      # if an non-valid permutation method is used
      stop(paste('\'',method,'\' is not an accepted permutation method; valid options are \'vasa\' or \'hungarian\'',sep=''))
      
    } # end of method choice
    
    # mapping is x->y
    # collate vectors from both hemispheres + sort mapping according to "reference" vector
    ref.lr <- c(ref.l,nroi.l+ref.r); rot.lr <- c(rot.l,nroi.l+rot.r);
    b <- sort(ref.lr,index.return=T); 
    ref.lr.sort <- ref.lr[b$ix]; rot.lr.sort <- rot.lr[b$ix];
    
    # verify that permutation worked (output should be vector with values 1:nroi = 1:(nroi_l+nroi_r))
    if (!all(sort(rot.lr.sort,decreasing=F)==c(1:nroi))) {
      #save.image('~/Desktop/perm_error.RData')
      browser("permutation error")
    }
    
    # verify that permutation does not map to itself
    if (!all(rot.lr.sort==c(1:nroi))) {
      r <- r+1
      perm.id[,r] = rot.lr.sort # if it doesn't, store it
    } else {
      c <- c+1
      print(paste('map to itself n. ',toString(c),sep=''))
    }
    
    # track progress
    if (r%%10==0) print(paste('permutation ',toString(r),' of ',toString(nrot),sep=''))
    
  }
  
  return(perm.id)
  
}
#coords <- read_csv("data/atlas_data/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_1mm.Centroid_RAS.csv")


perm_sphere_p = function(x, y, perm.id , corr.type='pearson') {
  
  nroi <- dim(perm.id)[1]  # number of regions
  nperm <- dim(perm.id)[2] # number of permutations
  
  rho.emp <- cor(x,y, method=corr.type)  # empirical correlation
  
  # permutation of measures
  x.perm <- y.perm <- array(NA,dim=c(nroi,nperm))
  for (r in 1:nperm) {
    for (i in 1:nroi) {
      x.perm[i,r] <- x[perm.id[i,r]]
      y.perm[i,r] <- y[perm.id[i,r]]
    }
  }
  
  # correlation to unpermuted measures
  rho.null.xy = rho.null.yx = vector(length=nperm)
  for (r in 1:nperm) {
    rho.null.xy[r] <- cor(x.perm[,r],y,method=corr.type)
    rho.null.yx[r] <- cor(y.perm[,r],x,method=corr.type)
  }
  
  # p-value definition depends on the sign of the empirical correlation
  if (rho.emp>0) {
    p.perm.xy <- sum(rho.null.xy>rho.emp)/nperm
    p.perm.yx <- sum(rho.null.yx>rho.emp)/nperm
  } else { 
    p.perm.xy <- sum(rho.null.xy<rho.emp)/nperm
    p.perm.yx <- sum(rho.null.yx<rho.emp)/nperm
  } 
  
  # return average p-value
  return((p.perm.xy+p.perm.yx)/2)
  
}
