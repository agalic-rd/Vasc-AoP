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
  
  return(list(data = list(normal = vasc_data, binned = vasc_binned_data)))
}


