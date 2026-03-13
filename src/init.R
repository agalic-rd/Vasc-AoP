#--------------#
####🔺infix ####
#--------------#

'%notin%' <- Negate('%in%')
`%ni%` <- `%notin%` # Alias
'%||%' <- function(x, y) if (is.null(x)) y else x

#-------------#
####🔺here ####
#-------------#

if (!"here" %in% as.data.frame(installed.packages())$Package) {
    renv::install("here", prompt = FALSE)
}

#-------------#
####🔺renv ####
#-------------#

source(here::here("src/renv/init.R"))

if (should_restore()) {
    safe_restore()
}

cli::cli_h1("[Vasc-AoP] Loading project")

#----------------------#
####🔺loading files ####
#----------------------#

cli::cli_h2("┗ [Vasc-AoP] Loading project's scripts")

scripts <- c("config.R", "theme.R", "viz.R", "utils.R", "stats.R", "data.R")
void_ <- lapply(scripts, \(script) source(here::here("src", script)))

rm(scripts, void_)

#-------------------------#
####🔺loading packages ####
#-------------------------#

cli::cli_h2("┗ [Vasc-AoP] Loading project's packages")

load_all_packages()

cli::cli_h1("[Vasc-AoP] Ready!")
