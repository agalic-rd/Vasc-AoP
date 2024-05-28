####╔═══════     ══════╗####
####💠 Data Wrangling 💠####
####╚═══════     ══════╝####

cli_h2("┗ [SCRIPTS] Loading data ingestion functions")


read_vasc_raw <- function(path = configs$data$vasc_raw_data) {
  
  col_names <- (
    openxlsx2::read_xlsx(path, rows = 1:2, col_names = FALSE, fill_merged_cells = TRUE)
    |> summarize(across(everything(), \(x) paste(na.omit(x), collapse = "__")))
    |> unlist()
  )
  
  vasc_data <- openxlsx2::read_xlsx(path, start_row = 3, col_names = FALSE)
  
  colnames(vasc_data) <- col_names
  
  return(
    vasc_data
    |> pivot_longer(contains("__"), names_sep = "__", names_to = c(".value", "Bins"))
    |> janitor::clean_names()
  )
}
