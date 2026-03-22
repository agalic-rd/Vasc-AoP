cli::cli_h2("┗ [Vasc-AoP] Setting project configs")

#----------------#
####🔺Options ####
#----------------#

options(
    verbose = FALSE,
    scipen = 999L,
    digits = 4L,
    na.action = "na.omit",
    contrasts = c("contr.sum", "contr.poly"),
    seed = 256,
    dplyr.summarise.inform = FALSE
)

set.seed(getOption("seed"))

#------------------------#
####🔺Project configs ####
#------------------------#

configs <- yaml::read_yaml(here::here("_configs.yml"), eval.expr = TRUE)

#----------------#
####🔺Masking ####
#----------------#

get <- base::get
filter <- dplyr::filter
