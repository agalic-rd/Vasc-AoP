####╔═════     ═════╗####
####💠Stats helpers💠####
####╚═════     ═════╝####

cli_h2("┗ [SCRIPTS] Loading stats functions")

#------------------------------------#
####🔺Summarizing data or a model ####
#------------------------------------#

distribution_summary <- function(data, dvs, between = "Condition") {
  data |> select(all_of(between), all_of(dvs)) |> 
    pivot_longer(all_of(dvs), names_to = "DV", values_to = "Value") |> 
    group_by(across(any_of(between)))  |> 
    group_map(
      \(d, g) datawizard::describe_distribution(group_by(d, DV), verbose = FALSE) |>
        mutate(
          Variance = SD^2,
          CoV = ifelse(SD / Mean > 1e4, NA_real_, SD / Mean),
          Variable = str_remove(.group, fixed("DV="))
        ) |> 
        add_column(g, .after = 1) |> 
        select("Variable", all_of(between), "Mean", "SD", "Variance", "CoV", "IQR", "Min", "Max", "Skewness", "Kurtosis", "n")
    ) |> 
    reduce(full_join, by = c("Variable", between, "Mean", "SD", "Variance", "CoV", "IQR", "Min", "Max", "Skewness", "Kurtosis", "n")) |> 
    arrange(Variable, across(any_of(between)))
}

get_model_based_outliers <- function(data, mod, mod_dharma, responses) {
  
  outliers <- get_data(mod) |> 
    rownames_to_column("ID") |> 
    filter(ID %in% DHARMa::outliers(mod_dharma)) |> 
    utils::type.convert(as.is = TRUE) 
  
  if (nrow(outliers) > 0) outliers <- semi_join(data, y = outliers) |> 
    select(-setdiff(responses, find_response(mod)))
  
  return(outliers)
}

#--------------------------------------------#
####🔺Extracting information from a model ####
#--------------------------------------------#

# Should we exponentiate the coefficients of a model (based on its link function)
should_exp <- \(mod) insight::get_family(mod)$link %in% c("log", "logit")

## Check which (if any) models have any NA as fixed effect coefficients (which signals that the model fitting failed silently)
has_na_coefs <- function(mods) {
  map_lgl(
    mods,
    \(mod) as.data.frame(summary(mod)$coefficients$cond) |> 
      mutate(across(where(is.character), \(x) na_if(x, "NaN"))) |> 
      lapply(\(x) any(is.na(x))) |> 
      flatten_lgl() |> 
      any()
  )
}
