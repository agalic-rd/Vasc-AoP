####╔════    ════╗####
####💠 Data Viz 💠####
####╚════    ════╝####

cli_h2("┗ [SCRIPTS] Loading vizualisation functions")

## Generating a boxplot for an individual gene, showing the main effect of a predictor (using the model fitted to this gene's data as input)
make_signif_boxplot <- function(
    mod, xaxis = "Condition", facet = NULL, cluster = "Mouse", add_cluster_averages = TRUE, subtitle = NULL, caption = NULL,
    scale = "link", adjust = "none", method = "pairwise", resp_name = NULL
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
  if (is.null(resp_name)) resp_name <- get_response_name(resp, "IHC")
  
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