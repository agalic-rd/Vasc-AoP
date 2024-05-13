####╔═════    ══════╗####
####💠Project Config💠####
####╚═════    ══════╝####

cli_h2("┗ [CONFIG] Loading external {.file _configs.yml} file")

configs <- load_yml(here("_configs.yml"))

options(
  scipen = 999L, 
  digits = 4L,
  na.action = "na.omit",
  contrasts = c("contr.sum", "contr.poly"),
  seed = 256
)

set.seed(getOption("seed"))

#------------------------#
####🔺Package options ####
#------------------------#

cli_h2("┗ [CONFIG] Setting packages options")


#----------------#
####🔺Masking ####
#----------------#

get <- base::get