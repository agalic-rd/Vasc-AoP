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

compute_fold_change <- function(mod) {
    return(
        insight::get_data(mod)
        |> dplyr::select(condition, dcq)
        |> tidyr::pivot_wider(names_from = condition, values_from = dcq, values_fn = \(x) mean(x, na.rm = TRUE))
        |> dplyr::summarize(fold = 2**(-1 * (IH - N)))
        |> dplyr::pull(fold)
    )
}

get_emmeans_data <- function(mod) {
    return(
        emmeans(mod, specs = "condition", type = "response")
        |> emmeans::contrast(method = "pairwise", adjust = "none", infer = TRUE)
        |> as.data.frame()
        |> tidyr::pivot_wider(names_from = contrast, values_from = estimate)
        |> dplyr::select(last_col(), LCB = lower.CL, UCB = upper.CL, p_value = p.value)
        |> dplyr::mutate(across(where(is.character), \(x) na_if(x, "NaN")))
    )
}

#-------------------------------#
####🔺Data loading functions ####
#-------------------------------#

## Loading the supplementary data
load_supplementary_data <- function() {
    
    res <- list()
    
    res$gene_data <- purrr::map(
        purrr::set_names(readxl::excel_sheets(configs$data$gene_data)), 
        \(sheet) readxl::read_excel(configs$data$gene_data, sheet) |> janitor::clean_names()
    )

    res$gene_data$fx <- res$gene_data$fx |> 
        tidyr::pivot_longer(ends_with("_effect"), names_pattern = "(.+)_effect", names_to = "effect", values_to = "gene") |> 
        tidyr::separate_longer_delim(gene, delim = "; ")
    
    return(res)
}

# Loading and shaping Transparization data
load_clearing_data <- function(path = configs$data$clearing_raw) {
    
    col_names <- (
        openxlsx2::read_xlsx(path, rows = 1:2, col_names = FALSE, fill_merged_cells = TRUE)
        |> dplyr::summarize(across(everything(), \(x) paste(na.omit(x), collapse = "__")))
        |> unlist()
    )
    
    clearing_data <- openxlsx2::read_xlsx(path, start_row = 3, col_names = FALSE) |> 
        filter(!if_all(everything(), is.na))
    
    colnames(clearing_data) <- col_names
    
    ## Adding totals (all depths)
    clearing_data <- (
        dplyr::bind_rows(
            clearing_data,
            clearing_data |> 
                dplyr::summarize(across(where(is.numeric), sum), .by = c(Stage, Mouse, Condition)) |> 
                dplyr::mutate(Level = "Total")
        )
        |> dplyr::filter(Mouse != "KA4N")
        |> dplyr::arrange(Stage, Condition, Level, Mouse)
    )
    
    # Binned variables
    clearing_binned_data <- clearing_data |> 
        dplyr::select(Level, Stage, Mouse, Condition, contains("__")) |> 
        tidyr::pivot_longer(contains("__"), names_sep = "__", names_to = c(".value", "Bins")) |> 
        janitor::clean_names()
    
    binned_col_names <- stringr::str_subset(colnames(clearing_binned_data), ".*_bin$")
    
    clearing_binned_data <- (
        purrr::map(
            rlang::set_names(binned_col_names),
            \(col) dplyr::select(clearing_binned_data, level, stage, mouse, condition, bins, any_of(col)) |> 
                dplyr::mutate(across(c(condition, stage, level), \(x) to_factor(x)))
        )
        |> purrr::reduce(\(x, acc) dplyr::left_join(acc, x, by = c("level", "stage", "mouse", "condition", "bins")))
        |> tibble::as_tibble()
    )
    
    # Unbinned variables
    clearing_data <- clearing_data |>
        dplyr::select(-contains("__")) |>
        janitor::clean_names() |> 
        dplyr::mutate(across(c(condition, stage, level), \(x) to_factor(x))) |> 
        tibble::as_tibble()
    
    # Data checks
    cerebellar_values_not_equal <- clearing_data |> 
        filter(level != "Total") |> 
        tidyr::pivot_wider(names_from = level, values_from = -c(level, stage, condition, mouse)) |> 
        dplyr::filter(
            cerebellar_volume_Deep != cerebellar_volume_Superficial 
            | cerebellar_area_Deep != cerebellar_area_Superficial, 
            .by = stage
        ) |> 
        dplyr::select(stage, mouse, condition, starts_with("cerebellar_volume"), starts_with("cerebellar_area"))
    
    if (nrow(cerebellar_values_not_equal) > 0) cli::cli_alert_warning("Cerebellar values are not equal for some mice")
    
    return(list(normal = clearing_data, binned = clearing_binned_data))
}

## Loading the PCR data
# @reprocess: if TRUE, re-process the raw data, else simply load it from the processed data
# @refit: if TRUE, refit the models and extract their predictions (careful, refitting takes upward to 10 minutes)
# @max_cq_clean: maximum Cq value allowed for the clean data
# @model: the model to fit to the data. Only required if reprocess = TRUE
load_pcr_data <- function(max_cq_clean = 33, reprocess = FALSE, model = NULL) {

    if (missing(model) && reprocess) {
        cli_abort("[DATA] Please provide a model to fit the data on.")
    }
    
    if (reprocess) {
        
        path <- configs$data$pcr_raw
        
        res <- list()
        
        ## Raw data
        res$raw <- (
            purrr::map(rlang::set_names(readxl::excel_sheets(path)), \(x) readxl::read_excel(path, sheet = x)) 
            |> purrr::list_rbind(names_to = "stage")
            |> janitor::clean_names()
            |> tidyr::extract(
                mouse,
                into = c("bloodline", "pup", "condition"),
                regex = "^(\\w{1,2})(\\d{1,2})([hH]+|[nN]+)$",
                convert = TRUE, remove = FALSE
            )
            |> dplyr::mutate(condition = if_else(condition == "H", "IH", condition))
            |> dplyr::mutate(across(c(condition, stage), \(x) to_factor(x)))
            |> dplyr::arrange(stage, gene, mouse)
            |> dplyr::select(mouse, experiment, stage, gene, condition, mean_cq, dcq, outlier)
            |> dplyr::mutate(id = row_number(), .before = 1)
        )
        
        ## Processed data (i.e. without unusable data points/groups)
        res$clean <- (
            res$raw
            ## Filtering useless data points
            |> dplyr::filter(!outlier)
            |> dplyr::filter(experiment != "B") # !!!
            |> tidyr::drop_na(mean_cq, dcq)
            |> dplyr::filter(mean_cq <= max_cq_clean)
            |> dplyr::filter(n() >= 3, .by = c(stage, gene, condition)) # At least 3 values per Condition
            |> dplyr::filter(n_distinct(condition) == 2, .by = c(stage, gene)) # Each Gene measured in both Conditions
            ## Arranging the data
            |> dplyr::select(id, mouse, experiment, stage, gene, condition, mean_cq, dcq)
            |> dplyr::mutate(across(c(condition, stage), \(x) to_factor(x)))
            |> dplyr::arrange(stage, gene, mouse)
        )
        
        ## Fitting the provided model to each Gene, for each Layer and Stage
        res$models <- (
            res$clean
            |> dplyr::group_split(stage, gene)
            |> purrr::map_dfr(
                \(d) dplyr::summarize(d, mod = list(model(pick(everything()))), .by = c(stage, gene)), 
                .progress = "Fitting models:"
            )
            |> dplyr::filter(!has_na_coefs(mod))
            |> dplyr::mutate(fold = purrr::map_dbl(mod, compute_fold_change))
            |> dplyr::select(stage, gene, fold, mod)
        )
        
        ## Extracting model predictions
        res$predictions <- (
            res$models
            |> dplyr::group_split(stage, gene)
            |> purrr::map_dfr(
                \(d) dplyr::mutate(d, get_emmeans_data(mod[[1]])), 
                .progress = "Extracting model predictions:"
            )
            |> dplyr::filter(!is.na(p_value))
            |> dplyr::mutate(expression = get_regulation_type(fold, p_value))
            |> dplyr::select(stage, gene, fold, expression, matches("-|/"), LCB, UCB, p_value)
        )

        saveRDS(res, configs$data$pcr_processed)
        
        return(res)
        
    } else {
        readRDS(configs$data$pcr_processed)
    }
}
