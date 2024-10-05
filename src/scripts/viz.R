####╔════    ════╗####
####💠 Data Viz 💠####
####╚════    ════╝####

#--------------------#
####🔺Correlation ####
#--------------------#

corr_matrix_plot <- function(dat, vars, title = "") {
  return(
    dat
    |> mutate(
      across(where(is.character), factor),
      across(where(is.factor), label_encoding)
    )
    |> correlation(select = vars, include_factors = TRUE, redundant = TRUE, method = "auto")
    |> rename(R = matches("^r$|^rho$"))
    |> mutate(across(matches("Parameter[1-2]"), \(x) factor(x, levels = vars)))
    |> ggplot(aes(x = Parameter1, y = Parameter2))
      + geom_tile(aes(fill = R), colour = "white", linewidth = 1.2, stat = "identity")
      + geom_text(aes(label = round(R, 2), colour = abs(R) > 0.5), size = rel(4.5))
      + scale_color_manual(values = c("black", "white"))
      + scale_fill_gradient2(na.value = "white", breaks = seq(-1, 1, 0.2), limits = c(-1, 1))
      + scale_x_discrete(position = "top")
      + scale_y_discrete(limits = rev)
      + guides(fill = guide_colourbar(title = "R", barheight = rel(17), title.hjust = 0.15), colour = "none")
      + labs(title = title)
      + theme(
        plot.title = element_markdown(hjust = 0.5)
        , axis.title.x = element_blank()
        , axis.title.y = element_blank()
        , axis.text.x = element_text(face = "bold", angle = 30, hjust = 0, size = 8)
        , axis.text.y = element_text(face = "bold", angle = 45, hjust = 1, size = 8)
      )
  )
}

#-----------------#
####🔺Boxplots ####
#-----------------#

## Generating a boxplot for an individual gene, showing the main effect of a predictor (using the model fitted to this gene's data as input)
make_signif_boxplot <- function(
    mod, xaxis = "condition", facet = NULL, cluster = "mouse", add_cluster_averages = TRUE, subtitle = NULL, caption = NULL, 
    invert_DCq = TRUE, scale = "link", adjust = "none", method = "pairwise", resp_name = NULL
) {
  
  get_n_units <- function(df) {
    if(!is.null(cluster) && cluster %in% colnames(df)) return(length(unique(df[[cluster]])))
    else return(dplyr::tally(df))
  }
  
  dat <- insight::get_data(mod)
  
  if (!is.null(cluster) && cluster %ni% colnames(dat)) {
    cluster <- NULL
    add_cluster_averages <- FALSE
  }
  
  if (!is.null(cluster) 
      && cluster %in% colnames(dat) 
      && dat |> group_by(across(any_of(c(xaxis, facet, cluster)))) |> count() |> filter(n > 1) |> nrow() == 0
  ) {
    cluster <- NULL
    add_cluster_averages <- FALSE
  }
  
  resp <- insight::find_response(mod)
  if (is.null(resp_name)) resp_name <- get_response_name(resp)
  
  is_DCq <- resp %in% c("DCq", "DCt", "dct", "dcq")
  if (is_DCq && invert_DCq) {
    dat <- dat |> mutate(DCq = -1 * DCq)
    resp_name <- "-1 * DCq"
  }
  
  ## Making sure the variables of interest are contrasts for emmeans
  dat <- dat |> mutate(across(c(any_of(c(xaxis, facet)) & where(\(c) !is.factor(c))), as.factor))
  
  extra_dat <- dat |> group_by(across(any_of(c(xaxis, facet)))) |> summarize(N = str_glue("N = {get_n_units(pick(everything()))}")) |> ungroup()
  
  max <- max(dat[[resp]])
  min <- min(dat[[resp]])
  amp <- abs(max - min)
  
  if(adjust == "none") correction <- "(uncorrected)"
  else correction <- str_glue("({adjust} corrected)")
  
  # -----------[ Contrasts ]----------- #
  
  specs <- paste0(" ~ ", xaxis)
  if(!is.null(facet)) specs <- paste0(specs, " | ", facet)
  specs <- as.formula(specs)
  
  emms <- emmeans::emmeans(mod, specs = specs, type = "response", data = insight::get_data(mod))
  if (tolower(scale) %in% c("response", "resp")) emm <- regrid(emm, transform = "response")
  
  contrasts <- emmeans::contrast(emms, method = method, adjust = adjust, infer = TRUE) |> 
    as_tibble() |> 
    rename(Contrast = contrast) |> 
    tidyr::extract(col = Contrast, into = c("X1", "X2"), regex = "(.*) [- | /] (.*)", remove = FALSE)
  
  p_data_contrasts <- (
    contrasts
    |> group_by(across(any_of(c(facet))))
    |> mutate(
      x1 = match(X1, levels(dat[[xaxis]])),
      x2 = match(X2, levels(dat[[xaxis]])),
      p.signif = label_pval(p_value)
    ) 
    |> arrange(x.diff := abs(x2 - x1))
    |> mutate(
      step = 1:n(),
      pos.x = (x2 + x1) * 0.5,
      pos.y = max + step * 0.1 * (max - min)
    ) 
    |> ungroup()
    |> filter(p.signif <= .05)
  )
  
  # -----------[ Plot ]----------- #
  
  plot <- (
    ggplot(dat, aes(x = .data[[xaxis]], y = .data[[resp]], color = .data[[xaxis]], fill = .data[[xaxis]]))
    + geom_boxplot(outlier.alpha = 0, size = 1.1, fill = NA)
    + stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), width = 0.75, linewidth = 1.1, linetype = "dotted")
    + { if (!is.null(cluster)) geom_jitter(
      data = \(x) x |> group_by(across(any_of(c(xaxis, facet)))) |> group_modify(\(d, g) slice_sample(d, n = min(nrow(d), 50))) |> ungroup(), 
      size = 1.5, width = 0.1, alpha = 0.3
    )
      else geom_jitter(
        data = \(x) x |> group_by(across(any_of(c(xaxis, facet)))) |> group_modify(\(d, g) slice_sample(d, n = min(nrow(d), 50))) |> ungroup(), 
        mapping = aes(fill = .data[[xaxis]]), shape = 23, color = color_text, size = 3, width = 0.1, alpha = 0.9
      )
    }
    + {if (add_cluster_averages) stat_summary(
      aes(group = .data[[cluster]], fill = .data[[xaxis]]), geom = "point", fun = mean, 
      size = ifelse(is.null(facet), 4, 3), shape = 23, color = color_text, alpha = 0.9, position = position_dodge(0.2)
    )}
    #+ ggrepel::geom_text_repel(aes(label = mouse), color = "black")
    + geom_errorbarh(
      data = p_data_contrasts, aes(xmin = x1, xmax = x2, y = pos.y), inherit.aes = FALSE, 
      color = "black", height = 0.03 * amp, linewidth = 0.5
    )
    + geom_text(
      data = p_data_contrasts, aes(x = pos.x, y = pos.y, label = p.signif), inherit.aes = FALSE,
      size = 5, color = "black", fontface = "bold", vjust = 0, hjust = 0.5, position = position_nudge(y = 0.02 * amp)
    )
    + geom_label(
      aes(y = min - 0.05 * amp, fontface = "bold", label = N, color = .data[[xaxis]]),
      data = extra_dat, fill = NA, size = 5, alpha = 0.7
    )
    + theme(
      legend.position = "none", 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      axis.title.x = element_blank(),
      plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "plain"),
      plot.caption = element_text(hjust = 0.5, face = "plain", size = 9)
    )
    + labs(y = resp_name)
    + {if(!is.null(subtitle)) labs(subtitle = subtitle)}
    + {if(!is.null(caption)) labs(caption = caption)}
    + {if (!is.null(facet)) facet_wrap( ~ .data[[facet]])}
    + {if (add_cluster_averages) labs(caption = str_glue("Small round points are individual measurements\n Diamonds represent {cluster}-averages"))}
  )
  
  return(plot)
}

## Generating a boxplot for an individual gene, showing interaction effects between two predictors (using the model fitted to this gene's data as input)
make_signif_boxplot_inter <- function(
    mod, pred1 = "condition", pred2, cluster = NULL, add_cluster_averages = FALSE, invert_DCq = TRUE, stage = NULL,
    scale = "link", adjust = "none", resp_name = NULL
) {
  
  get_n_units <- function(df) {
    if(!is.null(cluster) && cluster %in% colnames(df)) return(length(unique(df[[cluster]])))
    else return(dplyr::tally(df))
  }
  
  dat <- insight::get_data(mod)
  resp <- insight::find_response(mod)
  
  is_DCq <- resp %in% c("DCq", "DCt", "dct", "dcq")
  
  if (is_DCq && invert_DCq) dat <- dat |> mutate({{ resp }} := -1 * .data[[resp]]) # -1 * DCq for better legibility
  
  if (is.null(resp_name)) {
    if (is_DCq) resp_name <- ifelse(invert_DCq, str_glue("-1 * {resp}"), str_glue("{resp}"))
    else resp_name <- get_response_name(resp)
  }
  
  ## Making sure the variables of interest are contrasts for emmeans
  dat <- dat |> mutate(across(c(any_of(c(pred1, pred2)) & where(\(c) !is.factor(c))), as.factor))
  
  extra_dat <- dat |> group_by(across(any_of(c(pred1, pred2)))) |> summarize(N = str_glue("N = {get_n_units(pick(everything()))}")) |> ungroup()
  
  max <- max(dat[[resp]])
  min <- min(dat[[resp]])
  amp <- abs(max - min)
  
  # -----------[ Contrasts ]----------- #
  
  specs <- paste0(" ~ ", pred1)
  if(!is.null(pred2)) specs <- paste0(specs, " | ", pred2)
  specs <- as.formula(specs)
  
  emmeans <- emmeans::emmeans(mod, specs = specs, type = "response", data = insight::get_data(mod))
  if (tolower(scale) %in% c("response", "resp")) emmeans <- regrid(emmeans, transform = "response")
  
  contrasts <- emmeans::contrast(emmeans, method = "pairwise", adjust = adjust, infer = TRUE) |> 
    as.data.frame() |> 
    rename(Contrast = contrast) |> 
    tidyr::extract(col = Contrast, into = c("X1", "X2"), regex = "(.*) [- | /] (.*)", remove = FALSE)
  
  p_data_contrasts <- contrasts |>
    group_by(across(any_of(c(pred2)))) |>
    mutate(
      x1 = (match(.data[[pred2]], levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X1, levels(dat[[pred1]])),
      x2 = (match(.data[[pred2]], levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X2, levels(dat[[pred1]])),
      p.signif = label_pval(p_value)
    ) |>
    arrange(x.diff := abs(x2 - x1)) |>
    mutate(
      step = 1:n(),
      pos.x = (x2 + x1) * 0.5,
      pos.y = max + step * 0.1 * (max - min)
    ) |>
    ungroup()
  
  contrasts_interactions <- emmeans::contrast(emmeans, interaction = c("pairwise"), by = NULL, adjust = "none", infer = TRUE) |> 
    as.data.frame() |> 
    tidyr::extract(col = paste0(pred1, "_pairwise"), into = c("X1", "X2"), regex = "(.*) [- | /] (.*)", remove = FALSE) |> 
    tidyr::extract(col = paste0(pred2, "_pairwise"), into = c("F1", "F2"), regex = "(.*) [- | /] (.*)", remove = FALSE)
  
  p_data_interactions <- contrasts_interactions |>
    mutate(
      x1 = 0.5 * ((match(F1, levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X1, levels(dat[[pred1]])) +
                    (match(F1, levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X2, levels(dat[[pred1]]))),
      x2 = 0.5 * ((match(F2, levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X1, levels(dat[[pred1]])) +
                    (match(F2, levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) + match(X2, levels(dat[[pred1]]))),
      p.signif = label_pval(p_value)
    ) |>
    arrange(x.diff := abs(x2 - x1)) |>
    mutate(
      step = 1:n() + choose(length(unique(dat[[pred1]])), 2),
      pos.x = (x2 + x1) * 0.5,
      pos.y = max + step * 0.1 * (max - min)
    )
  
  # -----------[ Plot ]----------- #
  
  plot <- (
    ggplot(dat, aes(x = interaction(.data[[pred1]], .data[[pred2]], sep = "_"), y = .data[[resp]], color = .data[[pred1]]))
    + geom_boxplot(outlier.alpha = 0, size = 1.1, fill = NA)
    + stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), width = 0.75, size = 1.1, linetype = "dotted")
    + { 
      if (!is.null(cluster)) geom_jitter(
        data = \(x) x |> group_by(across(any_of(c(pred1, pred2)))) |> group_modify(\(d, g) slice_sample(d, n = min(nrow(d), 50))) |> ungroup(), 
        size = 1.5, width = 0.1, alpha = 0.3
      )
      else geom_jitter(
        data = \(x) x |> group_by(across(any_of(c(pred1, pred2)))) |> group_modify(\(d, g) slice_sample(d, n = min(nrow(d), 50))) |> ungroup(), 
        mapping = aes(fill = .data[[pred1]]), shape = 23, color = color_text, size = 3, width = 0.1, alpha = 0.9
      )
    }
    + {
      if (add_cluster_averages) stat_summary(
        aes(group = .data[[cluster]], fill = .data[[pred1]]), geom = "point", fun = mean, 
        size = 3, shape = 23, color = color_text, alpha = 0.9, position = position_dodge(0.2)
      )
    }
    #+ ggrepel::geom_text_repel(aes(label = mouse), color = "black")
    + geom_errorbarh(
      data = p_data_contrasts, aes(xmin = paste(X1, .data[[pred2]], sep = "_"), xmax = paste(X2, .data[[pred2]], sep = "_"), y = pos.y), inherit.aes = FALSE,
      color = "black", height = 0.02 * amp, size = 0.5
    )
    + geom_text(
      data = p_data_contrasts, aes(x = pos.x, y = pos.y, label = p.signif), inherit.aes = FALSE,
      size = 5, color = "black", fontface = "bold", vjust = 0, hjust = 0.5, position = position_nudge(y = 0.02 * amp)
    )
    + geom_label(
      aes(y = min - 0.05 * amp, fontface = "bold", label = N, color = .data[[pred1]]), 
      data = extra_dat, fill = NA, size = 5, alpha = 0.7
    )
    ## Interactions
    + geom_errorbarh(
      data = p_data_interactions, aes(xmin = x1, xmax = x2, y = pos.y), inherit.aes = FALSE,
      color = "black", height = 0.02 * amp, size = 0.5
    )
    + geom_text(
      data = p_data_interactions, aes(x = pos.x, y = pos.y, label = p.signif), inherit.aes = FALSE,
      size = 5, color = "black", fontface = "bold", vjust = 0, hjust = 0.5, position = position_nudge(y = 0.02 * amp)
    )
    + theme(
      legend.position = "none",
      plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "plain")
    )
    + labs(y = resp_name, x = str_c(pred1, " by ", pred2))
    + {if(!is.null(stage)) labs(subtitle = str_glue("{stage}"))}
    + {if (add_cluster_averages) labs(caption = str_glue("Small round points are individual measurements\n Diamonds represent {cluster}-averages"))}
    + scale_x_discrete(labels = \(l) str_replace(l, "_", "\n"))
  )
  
  return(plot)
}

#------------------#
####🔺Timelines ####
#------------------#

make_fold_timeline_plot <- function(
    dat, facet_rows = "Pathway", trans = "identity", 
    color_by = NULL, colors = colors_effect, size_boost = 1
) {
  
  origin <- do.call(trans, list(1))
  
  dat <- (
    dat
    |> mutate(fold_trans = do.call(trans, list(fold)))
    |> mutate(fold_amp = ifelse(
      max(fold_trans, na.rm = TRUE) - min(fold_trans, na.rm = TRUE) != 0, 
      max(fold_trans, na.rm = TRUE) - min(fold_trans, na.rm = TRUE), 
      mean(fold_trans, na.rm = TRUE)) * 0.1,
      .by = all_of(c(facet_rows, "stage"))
    )
  )
  
  timeline <- (
    ggplot(dat)
    + { if(is.null(color_by)) aes(x = gene, color = fold >= 1) else aes(x = gene, color = .data[[color_by]]) }
    + geom_linerange(aes(ymax = fold_trans), ymin = origin, linewidth = 2 + (size_boost * 0.5))
    + geom_hline(yintercept = origin, linewidth = 0.3, linetype = "dotted")
    + geom_text(aes(
        label = str_c(round(fold, 2), stars.pval(p_value) |> str_replace(fixed("."), ""), sep = " "), 
        y = ifelse(fold_trans > origin, fold_trans + fold_amp, fold_trans - fold_amp),
        hjust = ifelse(fold > 1, 0, 1)
      ),
      vjust = 0.5, angle = 0, size = 2 + (size_boost * 0.25), check_overlap = TRUE
    )
    + scale_color_manual(" ", values = colors)
    + scale_y_continuous(breaks = c(0,1,2,3), expand = expansion(mult = 1.01 * (1 + (size_boost/100))))
    + scale_x_discrete(expand = expansion(add = 1 * size_boost), limits = \(x) rev(x))
    + labs(
      x = "",
      y = ifelse(trans != "identity", str_glue("Fold Change *({trans} scale)*"), "Fold Change")
    )
    + coord_flip()
    + facet_grid(
      vars(.data[[facet_rows]]), vars(stage), 
      scales = "free_y", space = "free_y", labeller = label_wrap_gen(width = 12, multi_line = TRUE)
    )
    + { if (!is.null(color_by)) guides(color = guide_legend(title = color_by)) }
    + theme(
      legend.position = ifelse(is.null(color_by), "none", "bottom")
      , axis.text.x = element_blank()
      , axis.title.x = element_markdown(size = 9)
      , axis.text.y = element_text(size = 7)
      , strip.text = element_text(size = 5 * size_boost)
      , plot.title = element_markdown(size = 9, face = "plain", vjust = 1, hjust = 0.5)
    )
  )
  
  return(timeline)
}

#----------------#
####🔺Heatmap ####
#----------------#

make_heatmap <- function(data, xaxis, facet) {
  
  max_upreg <- data |> filter(fold >= 1) |> pull(fold) |> max()
  
  return(
    ggplot(data, aes(x = .data[[xaxis]], y = gene))
    + scale_fill_gradient(
      name = regulation_type$UPREG, 
      low = "seagreen2", high = "seagreen4", limits = c(1, max_upreg), 
      trans = scales::log10_trans(), labels = \(x) round(x, 2)
    )
    + geom_tile(data = \(d) filter(d, fold >= 1), aes(fill = fold), colour = "white")
    + new_scale("fill")
    + scale_fill_gradient(
      name = regulation_type$DOWNREG, 
      low = "firebrick4", high = "firebrick2", limits = c(0, 1), 
      labels = \(x) round(x, 2)
    )
    + geom_tile(data = \(d) filter(d, fold < 1), aes(fill = fold), colour = "white")
    + geom_text(
      aes(label = paste(round(fold, 2), stars.pval(p_value), sep = " ")), size = 3, colour = "white", fontface = "bold", 
      check_overlap = TRUE
    )
    + facet_grid(cols = vars(.data[[facet]]), scales = "free_x", space = "free_x")
    + theme_light_mar
    + theme(
      legend.title = element_markdown(face = "bold", vjust = 0.80),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
      axis.text.y = element_text(hjust = 0),
      axis.ticks = element_line(linewidth = 0.4)
    )
    + labs(x = xaxis, y = "")
  )
}
