####╔═════   ═════╗####
####💠Loading Data💠####
####╚═════   ═════╝####

cli_h2("┗ [SCRIPTS] Loading data ingestion functions")

#-------------------------#
####🔺Helper functions ####
#-------------------------#

## Turn a variable into a factor, ordered based on the level order defined within the relevant sheet of the provided data dictionary
to_factor <- function(var, dict_path = configs$data$data_dict) {
  return(factor(var, read_excel(dict_path, sheet = deparse(substitute(var)))$Name))
}

## Helper functions to associate a regulation status (down or up-regulated) to each gene

regulation_type <- list(
  NOT_REG = "Not Regulated", 
  MAYBE_UPREG = "Maybe Upregulated", 
  UPREG = "Upregulated",
  MAYBE_DOWNREG = "Maybe Downregulated",
  DOWNREG = "Downregulated"
)

get_regulation_type <- function(fold, p_value) {
  case_when(
    p_value <= .05 & fold < 1 ~ regulation_type$DOWNREG,
    p_value <= .05 & fold > 1 ~ regulation_type$UPREG,
    is.na(p_value) | is.na(fold) ~ NA_character_,
    .default = regulation_type$NOT_REG
  )
}


#-------------------------------#
####🔺Data loading functions ####
#-------------------------------#

## Loading the supplementary data
load_supplementary_data <- function() {
  
  res <- list()
  
  res$gene_data <- map(
    set_names(excel_sheets(configs$data$gene_data)), 
    \(sheet) read_excel(configs$data$gene_data, sheet) |> janitor::clean_names()
  )

  res$gene_data$fx <- res$gene_data$fx |> 
    tidyr::pivot_longer(ends_with("_effect"), names_pattern = "(.+)_effect", names_to = "effect", values_to = "gene") |> 
    tidyr::separate_longer_delim(gene, delim = "; ")
  
  return(res)
}

# Loading and shaping Transparization data
load_clearing_data <- function(path = configs$data$IHC$clearing_raw) {
  
  col_names <- (
    openxlsx2::read_xlsx(path, rows = 1:2, col_names = FALSE, fill_merged_cells = TRUE)
    |> summarize(across(everything(), \(x) paste(na.omit(x), collapse = "__")))
    |> unlist()
  )
  
  clearing_data <- openxlsx2::read_xlsx(path, start_row = 3, col_names = FALSE) |> 
    filter(!if_all(everything(), is.na))
  
  colnames(clearing_data) <- col_names
  
  ## Adding totals (all depths)
  clearing_data <- (
    bind_rows(
      clearing_data,
      clearing_data |> 
        summarize(across(where(is.numeric), sum), .by = c(Stage, Mouse, Condition)) |> 
        mutate(Level = "Total")
    )
    |> filter(Mouse != "KA4N")
    |> arrange(Stage, Condition, Level, Mouse)
  )
  
  # Binned variables
  clearing_binned_data <- clearing_data |> 
    select(Level, Stage, Mouse, Condition, contains("__")) |> 
    pivot_longer(contains("__"), names_sep = "__", names_to = c(".value", "Bins")) |> 
    janitor::clean_names()
  
  binned_col_names <- stringr::str_subset(colnames(clearing_binned_data), ".*_bin$")
  
  clearing_binned_data <- (
    purrr::map(
      rlang::set_names(binned_col_names),
      \(col) select(clearing_binned_data, level, stage, mouse, condition, bins, any_of(col)) |> 
        mutate(across(c(condition, stage), \(x) to_factor(x)))
    )
    |> purrr::reduce(\(x, acc) left_join(acc, x, by = c("level", "stage", "mouse", "condition", "bins")))
  )
  
  # Unbinned variables
  clearing_data <- clearing_data |>
    select(-contains("__")) |>
    janitor::clean_names() |> 
    mutate(across(c(condition, stage), \(x) to_factor(x))) |> 
    mutate(
      volume = volume / 1e6,
      
    )
  
  # Data checks (TODO: move to their own function)
  cerebellar_values_not_equal <- clearing_data |> 
    filter(level != "Total") |> 
    pivot_wider(names_from = level, values_from = -c(level, stage, condition, mouse)) |> 
    filter(
      cerebellar_volume_Deep != cerebellar_volume_Superficial 
      | cerebellar_area_Deep != cerebellar_area_Superficial, 
      .by = stage
    ) |> 
    select(stage, mouse, condition, starts_with("cerebellar_volume"), starts_with("cerebellar_area"))
  
  if (nrow(cerebellar_values_not_equal) > 0) cli::cli_alert_warning("Cerebellar values are not equal for some mice")
  
  return(list(normal = clearing_data, binned = clearing_binned_data))
}

## Loading the PCR data
# @target: which data to load (i.e. OS or ND)
# @reprocess: if TRUE, re-process the raw data, else simply load it from the processed data
# @refit: if TRUE, refit the models and extract their predictions (careful, refitting takes upward to 10 minutes)
# @max_cq_clean: maximum Cq value allowed for the clean data
# @model: the model to fit to the data. Only required if reprocess = TRUE
load_pcr_data <- function(target, reprocess = FALSE, max_cq_clean = 33, refit = FALSE, model) {
  
  possible_targets <- configs$data$PCR |> 
    names() |> 
    keep(\(x) str_detect(x, "_raw$|_processed$")) |> 
    str_split_i("_", 1) |> 
    unique()
  
  if (missing(target) || target %ni% possible_targets) 
    cli_abort("[DATA] Incorrect {.var target} provided: possible values are {.pkg {possible_targets}}")
  
  if (reprocess) {
    
    path <- configs$data$PCR[[str_glue("{target}_raw")]]
    
    res <- list()
    
    ## Raw data
    res$raw <- (
      map(set_names(excel_sheets(path)), \(x) read_excel(path, sheet = x)) 
      |> list_rbind(names_to = "stage")
      |> janitor::clean_names()
      |> tidyr::extract(
        mouse,
        into = c("bloodline", "pup", "condition"),
        regex = "^(\\w{1,2})(\\d{1,2})([hH]+|[nN]+)$",
        convert = TRUE, remove = FALSE
      )
      |> mutate(condition = if_else(condition == "H", "IH", condition))
      |> mutate(across(c(condition, stage), \(x) to_factor(x)))
      |> arrange(stage, gene, mouse)
      |> select(mouse, experiment, stage, gene, condition, mean_cq, dcq, outlier)
      |> mutate(id = row_number(), .before = 1)
    )
    
    ## Processed data (i.e. without unusable data points/groups)
    res$clean <- (
      res$raw
      ## Filtering useless data points
      |> filter(!outlier)
      |> filter(experiment != "B") # !!!
      |> drop_na(mean_cq, dcq)
      |> filter(mean_cq <= max_cq_clean)
      |> filter(n() >= 3, .by = c(stage, gene, condition)) # Checking that there are at least 3 values per Condition
      |> filter(n_distinct(condition) == 2, .by = c(stage, gene)) # Checking that each Gene has been measured in both Conditions
      ## Arranging the data
      |> select(id, mouse, experiment, stage, gene, condition, mean_cq, dcq)
      |> mutate(across(c(condition, stage), \(x) to_factor(x)))
      |> arrange(stage, gene, mouse)
    )
    
    if (refit) {
      
      if (missing(model)) 
        cli_abort("[DATA] Please provide a model to fit the data on.")
      
      compute_fold_change <- function(mod) {
        return(
          insight::get_data(mod)
          |> select(condition, dcq)
          |> pivot_wider(names_from = condition, values_from = dcq, values_fn = \(x) mean(x, na.rm = TRUE))
          |> summarize(fold = 2**(-1 * (IH - N)))
          |> pull(fold)
        )
      }
      
      ## Fitting the provided model to each Gene, for each Layer and Stage
      res$models <- (
        res$clean
        |> group_split(stage, gene)
        |> map_dfr(
          \(d) suppressMessages({summarize(d, mod = list(model(pick(everything()))), .by = c(stage, gene))}), 
          .progress = "Fitting models:"
        )
        |> filter(!has_na_coefs(mod))
        |> mutate(fold = map_dbl(mod, compute_fold_change))
        |> select(stage, gene, fold, mod)
      )
      
      get_emmeans_data <- function(mod) {
        return(
          emmeans(mod, specs = "condition", type = "response")
          |> contrast(method = "pairwise", adjust = "none", infer = TRUE)
          |> as.data.frame()
          |> pivot_wider(names_from = contrast, values_from = estimate)
          |> select(last_col(), LCB = lower.CL, UCB = upper.CL, p_value = p.value)
          |> mutate(across(where(is.character), \(x) na_if(x, "NaN")))
        )
      }
      
      ## Extracting model predictions
      res$predictions <- (
        res$models
        |> group_split(stage, gene)
        |> map_dfr(\(d) mutate(d, get_emmeans_data(mod[[1]])), .progress = "Extracting model predictions:")
        |> filter(!is.na(p_value))
        |> mutate(expression = get_regulation_type(fold, p_value))
        |> select(stage, gene, fold, expression, matches("-|/"), LCB, UCB, p_value)
      )
    }
    
    return(res)
    
  } else {
    readRDS(configs$data$PCR[[str_glue("{target}_processed")]])
  }
}
