cli::cli_h2("┗ [Vasc-AoP] Loading visualizations")

#--------------------#
####🔺Correlation ####
#--------------------#

corr_matrix_plot <- function(dat, vars, title = "") {
    return(
        dat |>
            dplyr::mutate(
                dplyr::across(where(is.character), factor),
                dplyr::across(where(is.factor), label_encoding)
            ) |>
            correlation(select = vars, include_factors = TRUE, redundant = TRUE, method = "auto") |>
            dplyr::rename(R = matches("^r$|^rho$")) |>
            dplyr::mutate(dplyr::across(matches("Parameter[1-2]"), \(x) {
                factor(x, levels = vars)
            })) |>
            ggplot(aes(x = Parameter1, y = Parameter2)) +
            geom_tile(aes(fill = R), colour = "white", linewidth = 1.2, stat = "identity") +
            geom_text(aes(label = round(R, 2), colour = abs(R) > 0.5), size = rel(4.5)) +
            scale_color_manual(values = c("black", "white")) +
            scale_fill_gradient2(na.value = "white", breaks = seq(-1, 1, 0.2), limits = c(-1, 1)) +
            scale_x_discrete(position = "top") +
            scale_y_discrete(limits = rev) +
            guides(fill = guide_colourbar(title = "R", barheight = rel(17), title.hjust = 0.15), colour = "none") +
            labs(title = title) +
            theme(
                plot.title = element_markdown(hjust = 0.5),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_text(face = "bold", angle = 30, hjust = 0, size = 8),
                axis.text.y = element_text(face = "bold", angle = 45, hjust = 1, size = 8)
            )
    )
}

#-----------------#
####🔺Boxplots ####
#-----------------#

## Generating a boxplot for an individual gene, showing the main effect of a predictor (using the model fitted to this gene's data as input)
make_signif_boxplot <- function(
    mod,
    xaxis = "condition",
    facet = NULL,
    cluster = "mouse",
    add_cluster_averages = TRUE,
    subtitle = NULL,
    caption = NULL,
    invert_DCq = TRUE,
    scale = "link",
    adjust = "none",
    method = "pairwise",
    resp_name = NULL,
    ncol = 2,
    only_facets = NULL
) {
    get_n_units <- function(df) {
        if (!is.null(cluster) && cluster %in% colnames(df)) {
            return(length(unique(df[[cluster]])))
        } else {
            return(dplyr::tally(df))
        }
    }

    dat <- insight::get_data(mod)

    if (!is.null(cluster) && cluster %ni% colnames(dat)) {
        cluster <- NULL
        add_cluster_averages <- FALSE
    }

    if (
        !is.null(cluster) &&
            cluster %in% colnames(dat) &&
            dat |>
                dplyr::group_by(dplyr::across(any_of(c(xaxis, facet, cluster)))) |>
                dplyr::count() |>
                dplyr::filter(n > 1) |>
                nrow() ==
                0
    ) {
        cluster <- NULL
        add_cluster_averages <- FALSE
    }

    resp <- insight::find_response(mod)
    if (is.null(resp_name)) {
        resp_name <- get_response_name(resp)
    }

    is_DCq <- resp %in% c("DCq", "DCt", "dct", "dcq")
    if (is_DCq && invert_DCq) {
        dat <- dat |> dplyr::mutate(DCq = -1 * DCq)
        resp_name <- "-1 * DCq"
    }

    dat <- dat |> dplyr::arrange(dplyr::across(any_of(c(facet, xaxis))))

    extra_dat <- dat |>
        dplyr::group_by(dplyr::across(any_of(c(xaxis, facet)))) |>
        dplyr::summarize(N = stringr::str_glue("N = {get_n_units(dplyr::pick(everything()))}")) |>
        dplyr::ungroup()

    if (!is.null(only_facets) && !is.null(facet)) {
        dat <- dat |> dplyr::filter(.data[[facet]] %in% only_facets) |> droplevels()
        extra_dat <- extra_dat |> dplyr::filter(.data[[facet]] %in% only_facets) |> droplevels()
    }

    max <- max(dat[[resp]])
    min <- min(dat[[resp]])
    amp <- abs(max - min)

    if (adjust == "none") {
        correction <- "(uncorrected)"
    } else {
        correction <- stringr::str_glue("({adjust} corrected)")
    }

    # -----------[ Contrasts ]----------- #

    specs <- paste0(" ~ ", xaxis)
    if (!is.null(facet)) {
        specs <- paste0(specs, " | ", facet)
    }
    specs <- as.formula(specs)

    emms <- emmeans::emmeans(mod, specs = specs, type = "response", data = droplevels(insight::get_data(mod)))
    if (tolower(scale) %in% c("response", "resp")) {
        emms <- regrid(emms, transform = "response")
    }

    contrasts <- emmeans::contrast(emms, method = method, adjust = adjust, infer = TRUE) |>
        tibble::as_tibble() |>
        dplyr::rename(Contrast = contrast) |>
        tidyr::extract(col = Contrast, into = c("X1", "X2"), regex = "(.*) [- | /] (.*)", remove = FALSE) |>
        dplyr::arrange(dplyr::across(any_of(facet)))

    if (!is.null(only_facets) && !is.null(facet)) {
        contrasts <- contrasts |> dplyr::filter(.data[[facet]] %in% only_facets) |> droplevels()
    }

    p_data_contrasts <- (contrasts |>
        dplyr::group_by(dplyr::across(any_of(c(facet)))) |>
        dplyr::mutate(
            x1 = match(X1, levels(dat[[xaxis]])),
            x2 = match(X2, levels(dat[[xaxis]])),
            p.signif = label_pval(p.value)
        ) |>
        dplyr::arrange(x.diff := abs(x2 - x1)) |>
        dplyr::mutate(step = 1:n(), pos.x = (x2 + x1) * 0.5, pos.y = max + step * 0.1 * (max - min)) |>
        dplyr::ungroup() |>
        dplyr::filter(p.signif <= .05))
    # -----------[ Plot ]----------- #

    plot <- (ggplot(dat, aes(x = .data[[xaxis]], y = .data[[resp]], color = .data[[xaxis]], fill = .data[[xaxis]])) +
        geom_boxplot(outlier.alpha = 0, size = 1.1, fill = NA) +
        stat_summary(
            fun = mean,
            geom = "errorbar",
            aes(ymax = after_stat(y), ymin = after_stat(y)),
            width = 0.75,
            linewidth = 1.1,
            linetype = "dotted"
        ) +
        {
            if (!is.null(cluster)) {
                geom_jitter(
                    data = \(x) {
                        x |>
                            dplyr::group_by(dplyr::across(any_of(c(xaxis, facet)))) |>
                            dplyr::group_modify(\(d, g) {
                                dplyr::slice_sample(d, n = min(nrow(d), 50))
                            }) |>
                            dplyr::ungroup()
                    },
                    size = 1.5,
                    width = 0.1,
                    alpha = 0.3
                )
            } else {
                geom_jitter(
                    data = \(x) {
                        x |>
                            dplyr::group_by(dplyr::across(any_of(c(xaxis, facet)))) |>
                            dplyr::group_modify(\(d, g) {
                                dplyr::slice_sample(d, n = min(nrow(d), 50))
                            }) |>
                            dplyr::ungroup()
                    },
                    mapping = aes(fill = .data[[xaxis]]),
                    shape = 23,
                    color = color_text,
                    size = 3,
                    width = 0.1,
                    alpha = 0.9
                )
            }
        } +
        {
            if (add_cluster_averages) {
                stat_summary(
                    aes(group = .data[[cluster]], fill = .data[[xaxis]]),
                    geom = "point",
                    fun = mean,
                    size = ifelse(is.null(facet), 4, 3),
                    shape = 23,
                    color = color_text,
                    alpha = 0.9,
                    position = position_dodge(0.2)
                )
            }
        } +
        #+ ggrepel::geom_text_repel(aes(label = mouse), color = "black")
        geom_errorbar(
            data = p_data_contrasts,
            aes(xmin = x1, xmax = x2, y = pos.y),
            inherit.aes = FALSE,
            color = "black",
            height = 0.03 * amp,
            linewidth = 0.5,
            orientation = "y"
        ) +
        geom_text(
            data = p_data_contrasts,
            aes(x = pos.x, y = pos.y, label = p.signif),
            inherit.aes = FALSE,
            size = 5,
            color = "black",
            fontface = "bold",
            vjust = 0,
            hjust = 0.5,
            position = position_nudge(y = 0.02 * amp)
        ) +
        geom_label(
            aes(y = min - 0.05 * amp, fontface = "bold", label = N, color = .data[[xaxis]]),
            data = extra_dat,
            fill = NA,
            size = 5,
            alpha = 0.7
        ) +
        scale_y_continuous(labels = scales::scientific) +
        theme(
            legend.position = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "plain"),
            plot.caption = element_text(hjust = 0.5, face = "plain", size = 9)
        ) +
        labs(x = get_response_name(xaxis), y = resp_name) +
        {
            if (!is.null(subtitle)) labs(subtitle = subtitle)
        } +
        {
            if (!is.null(caption)) labs(caption = caption)
        } +
        {
            if (!is.null(facet)) facet_wrap(~ .data[[facet]], ncol = ncol)
        } +
        {
            if (add_cluster_averages) {
                labs(
                    caption = stringr::str_glue(
                        "Small round points are individual measurements\n Diamonds represent {cluster}-averages"
                    )
                )
            }
        })

    return(plot)
}

## Generating a boxplot for an individual gene, showing interaction effects between two predictors (using the model fitted to this gene's data as input)
make_signif_boxplot_inter <- function(
    mod,
    pred1 = "condition",
    pred2,
    facet = NULL,
    cluster = NULL,
    add_cluster_averages = FALSE,
    invert_DCq = TRUE,
    stage = NULL,
    add_interactions = FALSE,
    scale = "link",
    adjust = "none",
    resp_name = NULL,
    max_points = 50,
    ncol = 2,
    only_signif = FALSE,
    only_facets = NULL,
    save_to_subfolder = NULL
) {
    get_n_units <- function(df) {
        if (!is.null(cluster) && cluster %in% colnames(df)) {
            return(length(unique(df[[cluster]])))
        } else {
            return(dplyr::tally(df))
        }
    }

    dat <- insight::get_data(mod)
    resp <- insight::find_response(mod)

    is_DCq <- resp %in% c("DCq", "DCt", "dct", "dcq")

    if (is_DCq && invert_DCq) {
        dat <- dat |> dplyr::mutate({{ resp }} := -1 * .data[[resp]])
    } # -1 * DCq for better legibility

    if (is.null(resp_name)) {
        if (is_DCq) {
            resp_name <- ifelse(invert_DCq, stringr::str_glue("-1 * {resp}"), stringr::str_glue("{resp}"))
        } else {
            resp_name <- get_response_name(resp)
        }
    }

    dat <- dat |> dplyr::arrange(dplyr::across(any_of(c(facet, pred2, pred1))))

    amp_data <- dat |>
        dplyr::group_by(dplyr::across(any_of(c(facet)))) |>
        dplyr::summarize(max = max(.data[[resp]]), min = min(.data[[resp]]), amp = abs(max - min)) |>
        dplyr::ungroup()

    extra_dat <- dat |>
        dplyr::group_by(dplyr::across(any_of(c(pred1, pred2, facet))), .drop = TRUE) |>
        dplyr::summarize(N = stringr::str_glue("N = {get_n_units(dplyr::pick(everything()))}")) |>
        # left_join(amp_data) |>
        dplyr::ungroup()

    if (!is.null(only_facets)) {
        dat <- dat |> dplyr::filter(.data[[facet]] %in% only_facets) |> droplevels()
        extra_dat <- extra_dat |> dplyr::filter(.data[[facet]] %in% only_facets) |> droplevels()
    }

    # Global min/max/amp
    max <- max(dat[[resp]])
    min <- min(dat[[resp]])
    amp <- abs(max - min)

    # -----------[ Contrasts ]----------- #

    specs <- paste0(" ~ ", pred1)
    if (!is.null(pred2)) {
        specs <- paste0(specs, " | ", pred2)
    }
    specs <- as.formula(specs)

    emms <- emmeans::emmeans(
        mod,
        specs = specs,
        type = "response",
        by = facet,
        data = droplevels(insight::get_data(mod))
    )
    if (tolower(scale) %in% c("response", "resp")) {
        emms <- regrid(emms, transform = "response")
    }

    contrasts <- emmeans::contrast(emms, method = "pairwise", adjust = adjust, infer = TRUE) |>
        as.data.frame() |>
        dplyr::rename(Contrast = contrast) |>
        tidyr::extract(col = Contrast, into = c("X1", "X2"), regex = "(.*) [- | /] (.*)", remove = FALSE) |>
        dplyr::arrange(dplyr::across(any_of(c(facet, pred2))))

    if (!is.null(only_facets)) {
        contrasts <- contrasts |> dplyr::filter(.data[[facet]] %in% only_facets) |> droplevels()
    }

    p_data_contrasts <- contrasts |>
        # left_join(amp_data) |>
        dplyr::group_by(dplyr::across(any_of(c(pred2, facet)))) |>
        dplyr::mutate(
            x1 = (match(.data[[pred2]], levels(dat[[pred2]])) - 1) *
                length(unique(dat[[pred1]])) +
                match(X1, levels(dat[[pred1]])),
            x2 = (match(.data[[pred2]], levels(dat[[pred2]])) - 1) *
                length(unique(dat[[pred1]])) +
                match(X2, levels(dat[[pred1]])),
            p.signif = label_pval(p.value)
        ) |>
        dplyr::arrange(x.diff := abs(x2 - x1)) |>
        dplyr::mutate(step = 1:n(), pos.x = (x2 + x1) * 0.5, pos.y = max + step * 0.1 * (max - min)) |>
        dplyr::ungroup()

    if (only_signif) {
        p_data_contrasts <- p_data_contrasts |>
            dplyr::group_by(dplyr::across(any_of(c(facet)))) |>
            dplyr::filter(any(p.signif <= .05)) |>
            dplyr::ungroup() |>
            droplevels()
    }

    if (add_interactions) {
        contrasts_interactions <- emmeans::contrast(
            emms,
            interaction = c("pairwise"),
            by = facet,
            adjust = "none",
            infer = TRUE
        ) |>
            as.data.frame() |>
            tidyr::extract(
                col = paste0(pred1, "_pairwise"),
                into = c("pred1_1", "pred1_2"),
                regex = "(.*) [- | /] (.*)",
                remove = FALSE
            ) |>
            tidyr::extract(
                col = paste0(pred2, "_pairwise"),
                into = c("pred2_1", "pred2_2"),
                regex = "(.*) [- | /] (.*)",
                remove = FALSE
            ) |>
            dplyr::arrange(dplyr::across(any_of(c(facet))))

        if (!is.null(only_facets)) {
            contrasts_interactions <- contrasts_interactions |>
                dplyr::filter(.data[[facet]] %in% only_facets) |>
                droplevels()
        }

        p_data_interactions <- contrasts_interactions |>
            # left_join(amp_data) |>
            dplyr::group_by(dplyr::across(any_of(c(facet)))) |>
            dplyr::mutate(
                x1 = 0.5 *
                    ((match(pred2_1, levels(dat[[pred2]])) - 1) *
                        length(unique(dat[[pred1]])) +
                        match(pred1_1, levels(dat[[pred1]])) +
                        (match(pred2_1, levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) +
                        match(pred1_2, levels(dat[[pred1]]))),
                x2 = 0.5 *
                    ((match(pred2_2, levels(dat[[pred2]])) - 1) *
                        length(unique(dat[[pred1]])) +
                        match(pred1_1, levels(dat[[pred1]])) +
                        (match(pred2_2, levels(dat[[pred2]])) - 1) * length(unique(dat[[pred1]])) +
                        match(pred1_2, levels(dat[[pred1]]))),
                p.signif = label_pval(p.value)
            ) |>
            dplyr::arrange(x.diff := abs(x2 - x1)) |>
            dplyr::mutate(
                step = 1:n() + choose(length(unique(dat[[pred1]])), 2),
                pos.x = (x2 + x1) * 0.5,
                pos.y = max + step * 0.1 * (max - min)
            ) |>
            dplyr::ungroup()

        if (only_signif) {
            p_data_interactions <- p_data_interactions |>
                dplyr::group_by(dplyr::across(any_of(c(facet)))) |>
                dplyr::filter(any(p.signif <= .05)) |>
                dplyr::ungroup() |>
                droplevels()
        }
    }

    if (only_signif) {
        dat <- dat |> dplyr::filter(.data[[facet]] %in% p_data_contrasts[[facet]]) |> droplevels()

        extra_dat <- extra_dat |> dplyr::filter(.data[[facet]] %in% p_data_contrasts[[facet]]) |> droplevels()
    }

    # -----------[ Plot ]----------- #

    plot <- (ggplot(
        dat,
        aes(x = interaction(.data[[pred1]], .data[[pred2]], sep = "_"), y = .data[[resp]], color = .data[[pred1]])
    ) +
        geom_boxplot(outlier.alpha = 0, size = 1.1, fill = NA) +
        stat_summary(
            fun = mean,
            geom = "errorbar",
            aes(ymax = after_stat(y), ymin = after_stat(y)),
            width = 0.75,
            size = 1.1,
            linetype = "dotted"
        ) +
        {
            if (!is.null(cluster)) {
                geom_jitter(
                    data = \(x) {
                        x |>
                            dplyr::group_by(dplyr::across(any_of(c(pred1, pred2)))) |>
                            dplyr::group_modify(\(d, g) {
                                dplyr::slice_sample(d, n = min(nrow(d), max_points))
                            }) |>
                            dplyr::ungroup()
                    },
                    size = 1.5,
                    width = 0.1,
                    alpha = 0.3
                )
            } else {
                geom_jitter(
                    data = \(x) {
                        x |>
                            dplyr::group_by(dplyr::across(any_of(c(pred1, pred2)))) |>
                            dplyr::group_modify(\(d, g) {
                                dplyr::slice_sample(d, n = min(nrow(d), max_points))
                            }) |>
                            dplyr::ungroup()
                    },
                    mapping = aes(fill = .data[[pred1]]),
                    shape = 23,
                    color = color_text,
                    size = 3,
                    width = 0.1,
                    alpha = 0.9
                )
            }
        } +
        {
            if (add_cluster_averages) {
                stat_summary(
                    aes(group = .data[[cluster]], fill = .data[[pred1]]),
                    geom = "point",
                    fun = mean,
                    size = 3,
                    shape = 23,
                    color = color_text,
                    alpha = 0.9,
                    position = position_dodge(0.2)
                )
            }
        } +
        geom_errorbar(
            data = p_data_contrasts,
            aes(xmin = paste(X1, .data[[pred2]], sep = "_"), xmax = paste(X2, .data[[pred2]], sep = "_"), y = pos.y),
            height = 0.02 * amp,
            inherit.aes = FALSE,
            color = "black",
            size = 0.5,
            orientation = "y"
        ) +
        geom_text(
            data = p_data_contrasts,
            aes(x = pos.x, y = pos.y, label = p.signif),
            inherit.aes = FALSE,
            size = 5,
            color = "black",
            fontface = "bold",
            vjust = 0,
            hjust = 0.5,
            position = position_nudge(y = 0.02 * amp)
        ) +
        geom_label(
            aes(y = min - 0.05 * amp, fontface = "bold", label = N, color = .data[[pred1]]),
            data = extra_dat,
            fill = NA,
            size = 5,
            alpha = 0.7
        ) +
        ## Interactions
        {
            if (add_interactions) {
                geom_errorbar(
                    data = p_data_interactions,
                    aes(xmin = x1, xmax = x2, y = pos.y),
                    height = 0.02 * amp,
                    inherit.aes = FALSE,
                    color = "black",
                    size = 0.5,
                    orientation = "y"
                )
            }
        } +
        {
            if (add_interactions) {
                geom_text(
                    data = p_data_interactions,
                    aes(x = pos.x, y = pos.y, label = p.signif),
                    inherit.aes = FALSE,
                    size = 5,
                    color = "black",
                    fontface = "bold",
                    vjust = 0,
                    hjust = 0.5,
                    position = position_nudge(y = 0.02 * amp)
                )
            }
        } +
        theme(legend.position = "none", plot.subtitle = ggtext::element_markdown(hjust = 0.5, face = "plain")) +
        scale_y_continuous(labels = scales::scientific) +
        labs(y = resp_name, x = stringr::str_c(get_response_name(pred1), " by ", get_response_name(pred2))) +
        {
            if (!is.null(stage)) labs(subtitle = stringr::str_glue("{stage}"))
        } +
        {
            if (!is.null(facet)) facet_wrap(~ .data[[facet]], ncol = ncol)
        } +
        {
            if (add_cluster_averages) {
                labs(
                    caption = stringr::str_glue(
                        "Small round points are individual measurements\n Diamonds represent {cluster}-averages"
                    )
                )
            }
        } +
        scale_x_discrete(labels = \(l) stringr::str_replace(l, "_", "\n")))

    if (!is.null(save_to_subfolder)) {
        n_facet <- length(unique(dat[[facet]]))
        fig_width <- 1 + 4 * n_facet
        fig_height <- 6 * ceiling(n_facet / ncol)
        model_name <- deparse(substitute(mod))
        fig_name <- stringr::str_glue("{model_name}_inter")
        if (!is.null(only_facets)) {
            fig_name <- stringr::str_glue("{fig_name}_{paste(only_facets, collapse = '_')}")
        }
        save_png(plot, fig_name, subfolder = save_to_subfolder, width = fig_width, height = fig_height)
    }

    return(plot)
}

#------------------#
####🔺Timelines ####
#------------------#

make_fold_timeline_plot <- function(
    dat,
    facet_rows = "Pathway",
    trans = "identity",
    color_by = NULL,
    colors = colors_effect,
    size_boost = 1
) {
    origin <- do.call(trans, list(1))

    dat <- (dat |>
        dplyr::mutate(fold_trans = do.call(trans, list(fold))) |>
        dplyr::mutate(
            fold_amp = ifelse(
                max(fold_trans, na.rm = TRUE) - min(fold_trans, na.rm = TRUE) != 0,
                max(fold_trans, na.rm = TRUE) - min(fold_trans, na.rm = TRUE),
                mean(fold_trans, na.rm = TRUE)
            ) *
                0.1,
            .by = all_of(c(facet_rows, "stage"))
        ))

    timeline <- (ggplot(dat) +
        {
            if (is.null(color_by)) {
                aes(x = gene, color = fold >= 1)
            } else {
                aes(x = gene, color = .data[[color_by]])
            }
        } +
        geom_linerange(aes(ymax = fold_trans), ymin = origin, linewidth = 2 + (size_boost * 0.5)) +
        geom_hline(yintercept = origin, linewidth = 0.3, linetype = "dotted") +
        geom_text(
            aes(
                label = stringr::str_c(round(fold, 2), stars_pval(p_value), sep = " "),
                y = ifelse(fold_trans > origin, fold_trans + fold_amp, fold_trans - fold_amp),
                hjust = ifelse(fold > 1, 0, 1)
            ),
            vjust = 0.5,
            angle = 0,
            size = 2 + (size_boost * 0.25),
            check_overlap = TRUE
        ) +
        scale_color_manual(" ", values = colors) +
        scale_y_continuous(breaks = c(0, 1, 2, 3), expand = expansion(mult = 1.01 * (1 + (size_boost / 100)))) +
        scale_x_discrete(expand = expansion(add = 1 * size_boost), limits = \(x) {
            rev(x)
        }) +
        labs(
            x = "",
            y = ifelse(trans != "identity", stringr::str_glue("Fold Change *({trans} scale)*"), "Fold Change")
        ) +
        coord_flip() +
        facet_grid(
            vars(.data[[facet_rows]]),
            vars(stage),
            scales = "free_y",
            space = "free_y",
            labeller = label_wrap_gen(width = 12, multi_line = TRUE)
        ) +
        {
            if (!is.null(color_by)) {
                guides(color = guide_legend(title = color_by))
            }
        } +
        theme(
            legend.position = ifelse(is.null(color_by), "none", "bottom"),
            axis.text.x = element_blank(),
            axis.title.x = element_markdown(size = 9),
            axis.text.y = element_text(size = 7),
            strip.text = element_text(size = 5 * size_boost),
            plot.title = element_markdown(size = 9, face = "plain", vjust = 1, hjust = 0.5)
        ))

    return(timeline)
}

#----------------#
####🔺Heatmap ####
#----------------#

make_heatmap <- function(data, xaxis, yaxis = "gene", facet) {
    max_upreg <- data |> dplyr::filter(fold >= 1) |> dplyr::pull(fold) |> max()

    return(
        ggplot(data, aes(x = .data[[xaxis]], y = .data[[yaxis]])) +
            scale_fill_gradient(
                name = regulation_type$UPREG,
                low = "#6fedb1",
                high = "#016937",
                limits = c(1, max_upreg),
                trans = scales::log10_trans(),
                labels = \(x) round(x, 2)
            ) +
            scale_color_identity() +
            geom_tile(data = \(d) dplyr::filter(d, fold >= 1), aes(fill = fold), colour = "white") +
            new_scale("fill") +
            scale_fill_gradient(
                name = regulation_type$DOWNREG,
                low = "#820318",
                high = "#d95d72",
                limits = c(0, 1),
                labels = \(x) round(x, 2)
            ) +
            geom_tile(data = \(d) dplyr::filter(d, fold < 1), aes(fill = fold), colour = "white") +
            geom_text(aes(label = tile_label), size = 3, colour = "white", fontface = "bold", check_overlap = TRUE) +
            facet_grid(cols = vars(.data[[facet]]), scales = "free_x", space = "free_x") +
            theme_light_mar +
            theme(
                legend.title = element_markdown(face = "bold", vjust = 0.80),
                legend.position = "bottom",
                legend.text = element_text(size = 11),
                strip.text = element_text(face = "bold"),
                axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
                axis.text.y = element_text(hjust = 0),
                axis.ticks = element_line(linewidth = 0.4)
            ) +
            labs(x = stringr::str_to_title(xaxis), y = "")
    )
}

#--------------------------#
####🔺Model diagnostics ####
#--------------------------#

make_acf_plot <- function(mod) {
    forecast::ggAcf(residuals(mod, type = "response", retype = "normalized"), color = "#1b6ca8") +
        geom_point() +
        labs(title = "Autocorrelation of residuals", subtitle = "Data (lines) should be inside the dashed area") +
        see::theme_lucid() +
        theme(
            plot.title = element_markdown(size = 15, margin = margin(0, 0, 0, 1)),
            plot.subtitle = element_markdown(size = 12, margin = margin(1, 0, 0, 0))
        )
}

ppc_plots <- function(mod, simulations, term = "condition", type = "fixed", is_count = NULL, max_cols_per_plot = 3) {
    Y <- insight::get_response(mod)
    n_unique <- dplyr::n_distinct(insight::get_data(mod)[[term]])

    if (is.null(is_count)) {
        is_count <- ifelse(insight::get_family(mod)$family |> stringr::str_detect("binom|poiss"), TRUE, FALSE)
    }

    # ppc_fun <- ifelse(is_count, bayesplot::ppc_bars, bayesplot::ppc_dens_overlay)
    ppc_fun_grouped <- ifelse(is_count, bayesplot::ppc_bars_grouped, bayesplot::ppc_dens_overlay_grouped)
    ppc_fun_pred_grouped <- bayesplot::ppc_intervals_grouped

    if (type %in% c("fixed", "fe")) {
        .term <- insight::get_predictors(mod)[[term]]

        # ppc_global <- ppc_fun(y = Y, yrep = simulations)
        if (is_count) {
            ppc_root <- bayesplot::ppc_rootogram(Y, simulations, style = "suspended")
        }

        ppc_grouped <- ppc_fun_grouped(Y, simulations, group = .term) +
            facet_wrap(~group, ncol = min(max_cols_per_plot, n_unique), scales = "free")
        ppc_pred_grouped <- ppc_fun_pred_grouped(Y, simulations, group = .term, prob_outer = 0.95) +
            facet_wrap(~group, ncol = min(max_cols_per_plot, n_unique), scales = "free")
    } else if (type %in% c("random", "re")) {
        .term <- insight::get_random(mod)[[term]]

        ppc_grouped <- ppc_fun_grouped(Y, simulations, group = .term) +
            facet_wrap(~group, ncol = min(max_cols_per_plot, n_unique), scales = "free")
        ppc_pred_grouped <- ppc_fun_pred_grouped(Y, simulations, group = .term, prob_outer = 0.95) +
            facet_wrap(~group, ncol = min(max_cols_per_plot, n_unique), scales = "free")
    }

    return(
        if (type %in% c("fixed", "fe")) {
            if (is_count) {
                (ppc_root / ppc_grouped / ppc_pred_grouped) +
                    plot_layout(guides = 'collect', ncol = 1, nrow = 3) +
                    # plot_annotation(title = "Simulation-based Posterior Predictive Checks", subtitle = stringr::str_glue("For [{term}]")) &
                    theme(legend.position = 'right', axis.title.x = element_blank())
            } else {
                (ppc_grouped / (ppc_pred_grouped + theme(axis.title.x = element_blank()))) +
                    plot_layout(ncol = 1, nrow = 2) +
                    # plot_annotation(title = "Simulation-based Posterior Predictive Checks", subtitle = stringr::str_glue("For [{term}]")) &
                    theme(legend.position = 'right')
            }
        } else {
            list(ppc_grouped, ppc_pred_grouped)
        }
    )
}

ppc_stat_plots <- function(
    mod,
    simulations,
    term = "condition",
    type = "fixed",
    stats = c("min", "max", "mean", "sd"),
    n_cols = 2,
    max_cols_per_plot = 5
) {
    n_unique <- dplyr::n_distinct(insight::get_data(mod)[[term]])

    if (type %in% c("fixed", "fe")) {
        .term <- insight::get_predictors(mod)[[term]]
    } else if (type %in% c("random", "re")) {
        .term <- insight::get_random(mod)[[term]]
    }

    return(
        patchwork::wrap_plots(
            purrr::map(stats, \(.x) {
                bayesplot::ppc_stat_grouped(
                    insight::get_response(mod),
                    simulations,
                    group = .term,
                    stat = .x,
                    facet_args = list(ncol = min(max_cols_per_plot, n_unique))
                ) +
                    scale_x_continuous(labels = \(l) signif(l, digits = 2))
            }),
            ncol = n_cols,
            guides = 'auto'
        ) +
            # plot_annotation(title = "Simulation-based Predictive Checks (on statistics)", subtitle = stringr::str_glue("For [{term}]")) &
            theme(legend.position = 'right', axis.text.x = element_text(size = rel(1.5), angle = 30, hjust = 1))
    )
}

#------------------------------#
####🔺Bar plots of Normoxia ####
#------------------------------#

## Bar + individual-point plot for a single metric normalized by cerebellar_volume.
## X axis = stage, bars = mean per stage, error bars = SEM, points = individual animals colored by sex
make_cerv_normalized_plot <- function(
    dat,
    metric,
    norm_var = "cerebellar_volume",
    xaxis = "stage",
    color_by = "sex",
    y_label = NULL,
    title = NULL
) {
    norm_col <- paste0(metric, "_norm")
    dat <- dat |> dplyr::mutate(!!norm_col := .data[[metric]] / .data[[norm_var]])

    if (is.null(y_label)) {
        y_label <- paste0(get_response_name(metric), " / ", get_response_name(norm_var))
    }

    ggplot(dat, aes(x = .data[[xaxis]], y = .data[[norm_col]])) +
        stat_summary(geom = "col", fun = "mean", width = 0.6, linewidth = 0.7, color = "black", fill = NA) +
        stat_summary(geom = "errorbar", fun.data = "mean_se", width = 0.2, linewidth = 0.5) +
        geom_point(
            aes(fill = .data[[color_by]], shape = .data[[color_by]]),
            color = "black",
            size = 2.5,
            stroke = 0.8,
            alpha = 0.95,
            position = position_dodge2(0.5)
        ) +
        scale_fill_manual(values = c("M" = "black", "F" = "white"), name = get_response_name(color_by)) +
        scale_shape_manual(values = c("M" = 24, "F" = 21), name = get_response_name(color_by)) +
        labs(x = get_response_name(xaxis), y = y_label, title = title) +
        theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "inside",
            legend.position.inside = c(0.99, 0.99),
            legend.justification = c("right", "top"),
            legend.direction = "horizontal",
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.key = element_rect(fill = "transparent", color = NA)
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
}

#------------------#
####🔺PCA plots ####
#------------------#

make_pca_plot <- function(dat, stage = "P4") {
    filtered_dat <- dat |> dplyr::filter(stage == !!stage & level == "Total")
    data_numeric <- filtered_dat |> dplyr::select(where(is.numeric))
    group <- filtered_dat |>
        dplyr::mutate(
            group = factor(paste0(stage, " ", condition), levels = unique(paste0(stage, " ", levels(condition))))
        ) |>
        dplyr::pull(group)

    pca_res <- PCA(data_numeric, scale.unit = TRUE, ncp = 4, graph = FALSE)

    fviz_pca_biplot(
        pca_res,
        geom.ind = "point",
        habillage = group,
        addEllipses = TRUE,
        ellipse.type = "convex",
        palette = colors_cond,
        col.var = "black",
        labelsize = 3,
        arrow.size = 0.1,
        pointsize = 1,
        repel = TRUE
    )
}
