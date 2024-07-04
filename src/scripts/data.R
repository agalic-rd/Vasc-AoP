####╔═════   ═════╗####
####💠Loading Data💠####
####╚═════   ═════╝####

cli_h2("┗ [SCRIPTS] Loading data ingestion functions")

# Loading and shaping Transparization data
load_vasc_data <- function(path = configs$data$vasc_raw_data) {
  
  col_names <- (
    openxlsx2::read_xlsx(path, rows = 1:2, col_names = FALSE, fill_merged_cells = TRUE)
    |> summarize(across(everything(), \(x) paste(na.omit(x), collapse = "__")))
    |> unlist()
  )
  
  vasc_data <- openxlsx2::read_xlsx(path, start_row = 3, col_names = FALSE) |> 
    filter(!if_all(everything(), is.na))
  
  colnames(vasc_data) <- col_names
  
  ## Adding totals (all depths)
  vasc_data <- (
    bind_rows(
      vasc_data,
      vasc_data |> 
        summarize(across(where(is.numeric), sum), .by = c(Stage, Mouse, Condition)) |> 
        mutate(Level = "Total")
    )
    |> filter(Mouse != "KA4N") # TEMP
    |> arrange(Stage, Condition, Level, Mouse)
  )
  
  # Binned variables
  vasc_binned_data <- vasc_data |> 
    select(Level, Stage, Mouse, Condition, contains("__")) |> 
    pivot_longer(contains("__"), names_sep = "__", names_to = c(".value", "Bins")) |> 
    janitor::clean_names()
  
  binned_col_names <- stringr::str_subset(colnames(vasc_binned_data), ".*_bin$")
  
  vasc_binned_data <- (
    purrr::map(
      binned_col_names,
      \(col) select(vasc_binned_data, level, stage, mouse, condition, bins, any_of(col))
    )
    |> set_names(binned_col_names)
  )
  
  # Unbinned variables
  vasc_data <- vasc_data |>
    select(-contains("__")) |>
    janitor::clean_names()
  
  # Data checks (TODO: move to their own function)
  cerebellar_values_not_equal <- vasc_data |> 
    filter(level != "Total") |> 
    pivot_wider(names_from = level, values_from = -c(level, stage, condition, mouse)) |> 
    filter(
      cerebellar_volume_Deep != cerebellar_volume_Superficial 
      | cerebellar_area_Deep != cerebellar_area_Superficial, 
      .by = stage
    ) |> 
    select(stage, mouse, condition, starts_with("cerebellar_volume"), starts_with("cerebellar_area"))
  
  if (nrow(cerebellar_values_not_equal) > 0) cli::cli_alert_warning("Cerebellar values are not equal for some mice")
  
  return(list(data = list(normal = vasc_data, binned = vasc_binned_data)))
}


