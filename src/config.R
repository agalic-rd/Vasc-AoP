cli::cli_h2("┗ [Vasc-AoP] Setting project configs")

#----------------#
####🔺Options ####
#----------------#

options(
    verbose = FALSE,
    bitmapType = "cairo",
    scipen = 999L,
    digits = 4L,
    na.action = "na.omit",
    contrasts = c("contr.sum", "contr.poly"),
    seed = 256,
    dplyr.summarise.inform = FALSE
)

if (nzchar(Sys.getenv("DISPLAY")) && capabilities("cairo")) {
    options(device = function(...) x11(type = "cairo", ...))
}

set.seed(getOption("seed"))

#------------------------#
####🔺Project configs ####
#------------------------#

configs <- yaml::read_yaml(here::here("_configs.yml"), eval.expr = TRUE)

#----------------#
####🔺Masking ####
#----------------#

get <- base::get
