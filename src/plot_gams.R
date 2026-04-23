
plot_gams_v1 <- function(gam_predictions, grad_df, 
                         atlas_geometry = readRDS("data/atlas_data/schaef_ggseg2.rds"),
                         gradient_cols = data.frame(gradient1 = c("#3F596D", "#D38A4E"), 
                                                    gradient2 =  c("#4682B4", "#781286"),  
                                                    gradient3 =  c("#8A6081", "#738518")),
                         add_shade = FALSE,
                         shade_alpha = 0.01,
                         shade_size = 0.01,
                         grad_name = TRUE,
                         biofinder_data, b_size = 7,
                         spintest = TRUE,
                         perms = readRDS("data/atlas_data/permutations_1000_hungarian.rds"),
                         figure_pat = "paper/figures") {
  
  require(ggside)
  require(ggpmisc)
  
  old <- theme_set(theme_bw(base_size = b_size))
  theme_update(panel.background = element_rect(fill = "transparent", colour = NA),
               plot.background = element_rect(fill = "transparent", colour = NA),
               legend.background = element_rect(fill = "transparent", colour = NA),
               legend.box.background = element_rect(fill = "transparent", colour = NA))
  
  
  legend_labs <-  c("Ab42/40", "Braak12", "Braak34", "Braak56")
  biof_path_plot <- biofinder_df |> filter(fmri_bl, !is.na(age), !is.na(pathology_ad)) |> 
    select(pathology_ad, ab_ratio, 
           starts_with("braak"), starts_with("cho")) |> 
    mutate(across(where(is.numeric) & !pathology_ad, scale)) 
  ab_mean <- attr(biof_path_plot$ab_ratio, c("scaled:center"))
  ab_sd <- attr(biof_path_plot$ab_ratio, c("scaled:scale"))
  ab_pos <- (0.08 - ab_mean)/ab_sd
  
  biof_path_plot <- biof_path_plot  |> 
    pivot_longer(-pathology_ad, names_to = "scaled_pat_measures", values_to = "value") %>% 
    ggplot(aes(pathology_ad, value, color = scaled_pat_measures)) +
    geom_smooth() +
    guides(color = guide_legend(
      ncol=1, byrow=TRUE,
      title.position="top", 
      title.hjust = 0)) +
    geom_hline(yintercept = ab_pos, linetype = 6) +
    geom_text(inherit.aes = FALSE, 
              data = data.frame(x = 0.85, y = ab_pos),
              aes(x = x, y = y),
              label = "Aβ+", nudge_y = 0.5, size = 5) +
    labs(x = "Pathology score", y = "Scaled value")+
    ggsci::scale_color_nejm(name = "Pathology", labels = legend_labs) +
    theme(legend.position = "right",
          legend.title.position = "top",
          #ggside.panel.scale = 0.2,
          legend.text = element_text(size = rel(0.6)),
          legend.title = element_text(size = rel(0.8)))  
  # geom_xsidehistogram(aes(x = pathology_ad),
  #                     color = "black",
  #                     fill = "lightgray",
  #                     bins = 20,
  #                     data = biofinder_data %>% filter(fmri_bl, !is.na(age), !is.na(pathology_ad)),
  #                      inherit.aes = FALSE) +
  # scale_xsidey_continuous(breaks = c(0, 100, 200), position = "left") +
  # ggside(x.pos = "bottom")
  
  pathology_plot <- biof_path_plot +
    theme(legend.position = "")  
  
  pat_leg <- ggpubr::get_legend(biof_path_plot)
  pat_leg <- ggpubr::as_ggplot(pat_leg)
  
  
  gradient_plots <- list()
  i = 1
  grad_char <- c("gradient1", "gradient2", "gradient3")
  for (grad in grad_char) {
    gradient_plots[[paste0(grad)]] <- grad_df %>% filter(study=="biofinder") %>% 
      #mutate(segregation = ifelse(segregation<0, 0, segregation)) %>% 
      right_join(atlas_geometry$atlas, by = "region") %>%
      ggplot() +
      geom_sf(aes(
        fill = .data[[grad]],
        geometry = geometry), linewidth= 0.1,
        show.legend = FALSE)+
      # geom_sf(data = network_geometry %>% drop_na(),
      #         aes(
      #           #fill = region,
      #           color = name,
      #           geometry = geometry
      #         ), alpha = 0, linewidth = 0.5,
      #         show.legend = FALSE) +
      (if (add_shade)
        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
       else NULL) +
      theme_void()+
      (if (grad_name)
        labs(title = case_when(
          grad == "gradient1" ~ "Sens-Assoc",
          grad == "gradient2" ~ "Vis-Mot",
          grad == "gradient3" ~ "Rep-Exec",
          TRUE ~ grad
        ))  else NULL) +
      theme(legend.position = "",
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(hjust = 0.5)
      ) +
      #guides(color = guide_legend(override.aes = list(size = 1))) +
      scale_fill_gradient2(
        low = gradient_cols[[grad]][1],
        mid = "white",
        high = gradient_cols[[grad]][2] 
      ) 
    i = i +1
  }
  
  
  brain_pat <- list()
  for (pat_grp in unique(cut(gam_predictions$pat_derivs$pathology_ad, 4))) {
    brain_pat[[pat_grp]] <- gam_predictions$pat_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(pathology_ad = cut(pathology_ad, 4)) %>% 
      filter(pathology_ad == pat_grp) %>% 
      group_by(region, pathology_ad) %>% 
      summarise(value=mean(value)) %>% 
      right_join(atlas_geometry$atlas, by = "region") %>%
      ggplot() +
      geom_sf(aes(
        fill = value,
        geometry = geometry), linewidth= 0.1,
        show.legend = FALSE)+
      (if (add_shade)
        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
       else NULL) +
      theme_void()+
      theme(legend.position = "",
            plot.margin = unit(c(0, 0, 0, 0), "npc"),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(color = "black", hjust = 0.5, size = rel(2)))+
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red"),
        # limits=c(-0.15, 0.15), oob=squish
      ) 
  }
  
  scatter_pat <- list()
  for (pat_grp in unique(cut(gam_predictions$pat_derivs$pathology_ad, 4))) {
    gam_predictions$pat_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(pathology_ad = cut(pathology_ad, 4)) %>% 
      filter(pathology_ad == pat_grp) %>% 
      group_by(region, pathology_ad) %>% 
      summarise(value=mean(value)) %>% 
      inner_join(grad_df) -> plot_df
    
    x_range_df <- data.frame(x = seq(min(plot_df$gradient1), max(plot_df$gradient1), length.out = 100))
    colnames(x_range_df)[1] <- "gradient1"
    y_range_df <- data.frame(y = seq(min(plot_df$value), max(plot_df$value), length.out = 100))
    colnames(y_range_df)[1] <- "value"
    
    x_axis_colorbar <- gradient_cols[[1]]
    
    scatter_pat[[pat_grp]] <- plot_df %>% 
      ggplot(aes(x = value, y = gradient1,
                 color = name)) +
      geom_point(alpha = 0.1, show.legend = FALSE) +
      stat_poly_line(se = FALSE, color = "#323232") +
      #xlim(-0.2, 0.05) +
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE, angle = 30), position = "bottom") +
      labs(
        #title = "Pathology AD {current_frame}",
        color = "Network",
        x = "FCS slopes averaged over pathology quartiles", 
        y = "") +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe())+
      scale_y_continuous(position = "right", limits = c(min(grad_df$gradient1), max(grad_df$gradient1))) +
      theme(axis.title.y = element_blank(),
            axis.text = element_text(size = rel(0.6))) +
      geom_xsidetile(data = y_range_df, aes(x = value, y = 0, fill = value), 
                     show.legend = FALSE, 
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red")
      ) +
      ggnewscale::new_scale_fill() +
      geom_ysidetile(data = x_range_df, aes(y = gradient1, x = 0, fill = gradient1), 
                     alpha = ifelse(pat_grp == unique(cut(gam_predictions$pat_derivs$pathology_ad, 4))[4], 1, 0),
                     show.legend = FALSE,
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = x_axis_colorbar[1],
        mid = "white",
        high = x_axis_colorbar[2] 
      ) +
      theme_ggside_void() +
      ggside(x.pos = "bottom", y.pos = "right") +
      theme(ggside.panel.scale = 0.05)
    
    if (spintest) {
      
      r <- cor(plot_df$gradient1, plot_df$value, method = "pearson")
      p_val <- perm_sphere_p(plot_df$gradient1, plot_df$value, perm.id = perms, corr.type='pearson')
      
      p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
      
      label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
      
      
      scatter_pat[[pat_grp]] <- scatter_pat[[pat_grp]] + ggpp::annotate(geom = "text_npc", label = label_corr,
                                                                        size =  3.6,
                                                                        npcx = "left", npcy = "top",
                                                                        family = "sans",
                                                                        parse = TRUE)
    } else {
      scatter_pat[[pat_grp]] <- scatter_pat[[pat_grp]] +       
        stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                       str_remove(after_stat(rr.confint.label), "95% CI "),
                                       sep = "*\" \"*")),
                     parse = TRUE, color = "#323232", label.x = "left", label.y = "top", size = 3.6) 
    }
    
  }
  
  
  brain_age <- list()
  for (age_grp in unique(cut(gam_predictions$age_derivs$age, 4))) {
    brain_age[[age_grp]] <- gam_predictions$age_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(age = cut(age, 4)) %>% 
      filter(age == age_grp) %>% 
      group_by(region, age) %>% 
      summarise(value=mean(value)) %>% 
      right_join(atlas_geometry$atlas, by = "region") %>%
      ggplot() +
      geom_sf(aes(
        #frame = pathology_ad,
        fill = value,
        geometry = geometry), linewidth= 0.1,
        show.legend = FALSE)+
      (if (add_shade)
        geom_sf(data = atlas_geometry$shade, size = shade_size, alpha = shade_alpha)
       else NULL) +
      theme_void()+
      #labs(title = "Pathology AD {current_frame}")+
      theme(legend.position = "",
            plot.margin = unit(c(0, 0, 0, 0), "npc"),
            panel.background = element_rect(fill = "transparent", colour = NA),
            plot.background = element_rect(fill = "transparent", colour = NA),
            legend.background = element_rect(fill = "transparent", colour = NA),
            legend.box.background = element_rect(fill = "transparent", colour = NA),
            plot.title = element_text(color = "black", hjust = 0.5, size = rel(2)))+
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red"),
        # limits=c(-0.15, 0.15), oob=squish
      ) 
  }
  
  scatter_age <- list()
  for (age_grp in unique(cut(gam_predictions$age_derivs$age, 4))) {
    gam_predictions$age_derivs %>% 
      pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>% 
      mutate(age = cut(age, 4)) %>% 
      filter(age == age_grp) %>% 
      group_by(region, age) %>% 
      summarise(value=mean(value)) %>% 
      inner_join(grad_df) -> plot_df
    
    x_range_df <- data.frame(x = seq(min(plot_df$gradient3), max(plot_df$gradient3), length.out = 100))
    colnames(x_range_df)[1] <- "gradient3"
    y_range_df <- data.frame(y = seq(min(plot_df$value), max(plot_df$value), length.out = 100))
    colnames(y_range_df)[1] <- "value"
    
    x_axis_colorbar <- gradient_cols[[3]]
    
    scatter_age[[age_grp]] <-  plot_df %>% 
      ggplot(aes(x = value, y = gradient3,
                 color = name)) +
      geom_point(alpha = 0.1, show.legend = FALSE) +
      stat_poly_line(se = FALSE, color = "#323232") +
      ylim(min(grad_df$gradient3), max(grad_df$gradient3)) +
      scale_x_continuous(guide = guide_axis(check.overlap = TRUE, angle = 30)#, position = "top"
      ) +
      labs(
        color = "Network",
        x = "FCS slopes averaged over age quartiles", 
        y = "") +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe()) +
      theme(axis.title.y = element_blank(),
            axis.text = element_text(size = rel(0.6))) +
      geom_xsidetile(data = y_range_df, aes(x = value, y = 0, fill = value), 
                     show.legend = FALSE, 
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = muted("blue"),
        mid = "white",
        high = muted("red")
      ) +
      ggnewscale::new_scale_fill() +
      geom_ysidetile(data = x_range_df, aes(y = gradient3, x = 0, fill = gradient3), 
                     alpha = ifelse(age_grp == unique(cut(gam_predictions$age_derivs$age, 4))[1], 1, 0),
                     show.legend = FALSE,
                     inherit.aes = FALSE) +
      scale_fill_gradient2(
        low = x_axis_colorbar[1],
        mid = "white",
        high = x_axis_colorbar[2] 
      ) +
      theme_ggside_void() +
      ggside(x.pos = "bottom", y.pos = "left") +
      theme(ggside.panel.scale = 0.05)
    
    
    if (spintest) {
      
      r <- cor(plot_df$gradient3, plot_df$value, method = "pearson")
      p_val <- perm_sphere_p(plot_df$gradient3, plot_df$value, perm.id = perms, corr.type='pearson')
      p_lab <- scales::label_pvalue(accuracy = 0.001, prefix = c("italic(p[spin]) < ", "italic(p[spin]) == ", "italic(p[spin]) > "))(p_val)
      label_corr <- paste0("italic(r) == ", r |> round(2),"*','~", p_lab)
      
      
      scatter_age[[age_grp]] <- scatter_age[[age_grp]] + ggpp::annotate(geom = "text_npc", label = label_corr,
                                                                        size = 3.6,
                                                                        npcx = "left", npcy = "top",
                                                                        family = "sans",
                                                                        parse = TRUE)
    } else {
      scatter_age[[age_grp]] <- scatter_age[[age_grp]] + 
        stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                       str_remove(after_stat(rr.confint.label), "95% CI "),
                                       sep = "*\" \"*")),
                     parse = TRUE, color = "#323232", label.x = "left", label.y = "top", size = 3.6) 
    }
    
  }
  
  
  quantile_trajectories <- gam_predictions$pat_pred %>% pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "predicted_fc") %>% 
    inner_join(grad_df %>% select(region, gradient1)) %>% 
    mutate(grad_grp = cut(gradient1, breaks = quantile(gradient1, probs = seq(0, 1, length.out = 21)))) %>% 
    group_by(pathology_ad, grad_grp) %>% 
    summarise(mean_pred_fc = mean(predicted_fc),
              mean_grad_value = mean(gradient1)) %>% 
    ggplot(aes(pathology_ad, mean_pred_fc, group = mean_grad_value, color = mean_grad_value)) +
    geom_line(linewidth = 1) +
    labs(x = "Pathology score", y = "Predicted FCS", color = "Mean SA\nin ventile") +
    #guides(colour = guide_colorbar(title.position="top", title.hjust = 0.0))+
    scale_color_gradient2(
      high = gradient_cols[[1]][2], mid = "white", low = gradient_cols[[1]][1]) +
    theme(
      legend.key.height = unit(1, "null"),
      legend.key.width = unit(0.02, "npc"),
      legend.box.spacing = unit(0.0025, "npc"),
      legend.title = element_blank(),
      legend.text = element_text(margin = margin(l = 0.8), size = rel(0.7))
    ) 
  
  
  r2_pat <- gam_predictions$pat_derivs %>% pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>%
    inner_join(grad_df) %>%
    group_by(pathology_ad) %>%
    summarise(grad_R2 = cor(value, gradient1)) %>%
    ggplot(aes(pathology_ad, grad_R2)) +
    geom_line() +
    labs(y = "Correlation (r)", x = "Pathology score") +
    theme(
      ggside.panel.scale = 0.2
    ) 
  # +
  #   geom_xsidehistogram(aes(x = pathology_ad), data = biofinder_data, inherit.aes = FALSE) +
  #   ggside(x.pos = "top")
  
  
  
  for(x_int in seq_range(gam_predictions$pat_pred$pathology_ad, 5)){
    quantile_trajectories <- quantile_trajectories +
      geom_vline(xintercept = x_int, linetype = "dashed")
    r2_pat <- r2_pat +
      geom_vline(xintercept = x_int, linetype = "dashed")
    pathology_plot <- pathology_plot +
      geom_vline(xintercept = x_int, linetype = "dashed")
  }
  
  age_range <- range(gam_predictions$age_pred$age)
  
  quantile_trajectories_term2 <- gam_predictions$age_pred %>% 
    pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "predicted_fc") %>% 
    inner_join(grad_df %>% select(region, gradient3)) %>% 
    mutate(grad_grp = cut(gradient3, breaks = quantile(gradient3, probs = seq(0, 1, length.out = 21)))) %>% 
    group_by(age, grad_grp) %>% 
    summarise(mean_pred_fc = mean(predicted_fc),
              mean_grad_value = mean(gradient3)) %>% 
    ggplot(aes(age, mean_pred_fc, group = mean_grad_value, color = mean_grad_value)) +
    geom_line(linewidth = 1) +
    labs(y = "Predicted FCS", x = "Age", color = "Mean RE\nin ventile") +
    scale_color_gradient2(
      high = gradient_cols[[3]][2], mid = "white", low = gradient_cols[[3]][1]) +
    #guides(colour = guide_colorbar(title.position="top", title.hjust = 0.0))+
    scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
    theme(
      legend.key.height = unit(1, "null"),
      legend.key.width = unit(0.02, "npc"),
      legend.box.spacing = unit(0.005, "npc"),
      legend.title = element_blank(),
      legend.text = element_text(margin = margin(l = 0.8), size = rel(0.7))
    ) 
  
  r2_term2 <- gam_predictions$age_derivs %>% pivot_longer(starts_with("7Networks"), names_to = "region", values_to = "value") %>%
    inner_join(grad_df ) %>%
    group_by(age) %>%
    summarise(grad_R2 = cor(value, gradient3)) %>%
    ggplot(aes(age, grad_R2)) +
    geom_line() +
    theme(ggside.panel.scale = 0.2) +
    labs(y = "Correlation (r)", x = "Age") 
  # geom_xsidehistogram(aes(x = age), data = biofinder_data |> 
  #                       filter(age > age_range[1], age < age_range[2]), 
  #                     inherit.aes = FALSE) +
  # ggside(x.pos = "bottom")
  
  
  for(x_int in seq_range(gam_predictions$age_pred$age, 5)){
    quantile_trajectories_term2 <- quantile_trajectories_term2 +
      geom_vline(xintercept = x_int, linetype = "dashed")
    r2_term2 <- r2_term2 +
      geom_vline(xintercept = x_int, linetype = "dashed")
  }
  
  
  
  a <- wrap_plots(c(scatter_pat,
                    list(gradient_plots[[1]]),
                    list(plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer()),
                    brain_pat, 
                    list(plot_spacer())
  ), 
  nrow = 3) +
    plot_layout(
      guides = "collect",           
      axis_titles = "collect",      
      axes = "collect",
      heights = c(0.5, -0.15, 0.5)          
    )
  
  a <- ggdraw() +
    draw_plot(a) +
    draw_figure_label("A")
  
  
  b <- 
    wrap_plots(list(pathology_plot + labs(tag = "B.1"), plot_spacer(),
                    r2_pat + labs(tag = "B.2"), r2_term2,
                    quantile_trajectories+ labs(tag = "B.3"), quantile_trajectories_term2),
               ncol = 2, byrow = TRUE) +
    plot_layout(
      axis_titles = "collect", axes = "collect"
      #guides = "collect" 
    )  & 
    theme(plot.tag.position  = c(0.175, 0.9),
          plot.tag = element_text(size = rel(0.6), hjust = 0, vjust = 0),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "transparent", colour = NA)
          #legend.box = "horizontal"
    ) 
  
  b <- ggdraw() +
    draw_plot(b) +
    draw_figure_label("B")
  
  c <- wrap_plots(c(list(plot_spacer()),
                    brain_age, 
                    list(plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer(),plot_spacer()),
                    list(gradient_plots[[3]]),
                    scatter_age), nrow = 3) +
    plot_layout(
      guides = "collect",           
      axis_titles = "collect",      
      axes = "collect",
      heights = c(0.5, -0.2, 0.5)
    )
  
  c <- ggdraw() +
    draw_plot(c) +
    draw_figure_label("C")
  
  get_net_legend <- function(b_size = 11, point_size = 3){
    #scale_factor <- 5
    x <- grad_df %>% filter(study == "biofinder") |> 
      ggplot(aes(gradient1, gradient3, color = name)) +
      geom_point(alpha = 0.5) +
      theme_bw(base_size = b_size) +
      guides(color = guide_legend(
        label.hjust=0,
        byrow = TRUE,
        nrow = 2, 
        override.aes = list(size = rel(point_size))
      )) +
      scale_color_manual(values = net_names %>% select(name, col) %>% deframe) +
      theme(
        legend.position = "bottom",
        #legend.key.size = unit(0.0, "cm"),
        legend.key.spacing.x = unit(0.2, "cm"),
        legend.key.spacing.y = unit(-1.25, "cm"),
        legend.direction = "horizontal",
        #legend.title.position = "",
        legend.text.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = rel(0.8), margin = margin(l = 0, r = 2, unit = "pt")),
        legend.background = element_blank()
      )
    
    leg <- ggpubr::get_legend(x)
    leg <- ggpubr::as_ggplot(leg)
    leg 
  }
  
  net_legend <- get_net_legend()
  
  age_hist <- biofinder_data |> filter(age > age_range[1], age < age_range[2]) |>  
    ggplot(aes(x = age)) +
    geom_histogram(bins = 20, color = "black", fill = "lightgray") +
    labs(y = "", x = "") +
    scale_y_continuous(breaks = c(0, 100), position = "right") +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
    )
  
  path_hist <- biofinder_data  |>
    ggplot(aes(x = pathology_ad)) +
    geom_histogram(bins = 20, color = "black", fill = "lightgray") +
    labs(y = "", x = "") +
    scale_y_continuous(breaks = c(0, 200), position = "left") +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank()
    )
  
  # 
  nonlin_p <-   ggdraw() +
    draw_plot(a, x = 0, y = 0.7, width = 1, height = 0.30) +
    draw_plot(b, x = 0, y = 0.30, width = 1, height = 0.40) +
    draw_plot(age_hist, x = 0.55, y = 0.55, width = 0.44, height = 0.065) +
    draw_plot(path_hist, x = 0.015, y = 0.66, width = 0.44, height = 0.065) +
    draw_plot(pat_leg, y = 0.585, x = 0.4, width = 0.2, height = 0.1) +
    draw_plot(c, x = 0, y = 0.0, width = 1, height = 0.30) +
    draw_plot(net_legend, x = 0.05, y = -0.025, width = 0.2, height = 0.1)+
    draw_text(text = c("Mean SA\nin ventile", 
                       "Mean RE\nin ventile"), x = c(0.485, 0.965), y = c(0.465, 0.465), size = 11) +
    draw_line(
      x = c(0.099, 0.01, 0.01),
      y = c(0.685, 0.725, 0.81), 
      linetype = 2
    ) +  
    draw_line(
      x = c(0.099 + 0.082, 0.01 + 0.19, 0.01 + 0.19),
      y = c(0.685, 0.725, 0.81), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.099 + 0.082*2, 0.01 + 0.19*2, 0.01 + 0.19*2),
      y = c(0.685, 0.725, 0.81), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.099 + 0.082*3, 0.01 + 0.19*3, 0.01 + 0.19*3),
      y = c(0.685, 0.725, 0.81), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.099 + 0.082*4, 0.01 + 0.19*4, 0.01 + 0.19*4),
      y = c(0.685, 0.725, 0.81), 
      linetype = 2
    ) + # NEDRE DEL
    draw_line(
      x = c(0.495 + 0.082, 0.04 + 0.19, 0.04 + 0.19),
      y = c(0.35, 0.285, 0.19), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.495 + 0.082*2, 0.04 + 0.19*2, 0.04 + 0.19*2),
      y = c(0.35, 0.285, 0.19), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.495 + 0.082*3, 0.04 + 0.19*3, 0.04 + 0.19*3),
      y = c(0.35, 0.285, 0.19), 
      linetype = 2
    ) +
    draw_line(
      x = c(0.495 + 0.082*4, 0.04 + 0.19*4, 0.04 + 0.19*4),
      y = c(0.35, 0.285, 0.19), 
      linetype = 2
    )+
    draw_line(
      x = c(0.495 + 0.082*5, 0.04 + 0.19*5, 0.04 + 0.19*5),
      y = c(0.35, 0.285, 0.19), 
      linetype = 2
    )
  nonlin_p
}
