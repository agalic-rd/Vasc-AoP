cli::cli_h2("┗ [Vasc-AoP] Loading stats helpers")

#------------------------------------#
####🔺Summarizing data or a model ####
#------------------------------------#

distribution_summary <- function(data, dvs, between = "Condition") {
    data |>
        dplyr::select(all_of(between), all_of(dvs)) |>
        tidyr::pivot_longer(all_of(dvs), names_to = "DV", values_to = "Value") |>
        dplyr::group_by(dplyr::across(any_of(between))) |>
        dplyr::group_map(\(d, g) {
            datawizard::describe_distribution(dplyr::group_by(d, DV), verbose = FALSE) |>
                dplyr::select(-"Variable") |>
                dplyr::rename(Variable = "DV") |>
                dplyr::mutate(Variance = SD^2, CoV = ifelse(SD / Mean > 1e4, NA_real_, SD / Mean)) |>
                tibble::add_column(g, .after = 1) |>
                dplyr::select(
                    "Variable",
                    all_of(between),
                    "Mean",
                    "SD",
                    "Variance",
                    "CoV",
                    "IQR",
                    "Min",
                    "Max",
                    "Skewness",
                    "Kurtosis",
                    "n"
                )
        }) |>
        purrr::reduce(
            dplyr::full_join,
            by = c(
                "Variable",
                between,
                "Mean",
                "SD",
                "Variance",
                "CoV",
                "IQR",
                "Min",
                "Max",
                "Skewness",
                "Kurtosis",
                "n"
            )
        ) |>
        dplyr::arrange(Variable, dplyr::across(any_of(between)))
}

get_model_based_outliers <- function(data, mod, mod_dharma, responses) {
    outliers <- get_data(mod) |>
        tibble::rownames_to_column("ID") |>
        dplyr::filter(ID %in% DHARMa::outliers(mod_dharma)) |>
        utils::type.convert(as.is = TRUE)

    if (nrow(outliers) > 0) {
        outliers <- dplyr::semi_join(data, y = outliers) |> dplyr::select(-setdiff(responses, find_response(mod)))
    }

    return(outliers)
}

#--------------------------------------------#
####🔺Extracting information from a model ####
#--------------------------------------------#

# Should we exponentiate the coefficients of a model (based on its link function)
should_exp <- \(mod) insight::get_family(mod)$link %in% c("log", "logit")

## Check which (if any) models have any NA as fixed effect coefficients (which signals that the model fitting failed silently)
has_na_coefs <- function(mods) {
    purrr::map_lgl(mods, \(mod) {
        as.data.frame(summary(mod)$coefficients$cond) |>
            dplyr::mutate(dplyr::across(where(is.character), \(x) dplyr::na_if(x, "NaN"))) |>
            lapply(\(x) anyNA(x)) |>
            purrr::flatten_lgl() |>
            any()
    })
}
