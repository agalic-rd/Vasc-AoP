####╔════    ════╗####
####💠 Data Viz 💠####
####╚════    ════╝####

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
      p.signif = label_pval(p.value)
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
      p.signif = label_pval(p.value)
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
      p.signif = label_pval(p.value)
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