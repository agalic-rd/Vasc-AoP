####╔══════     ══════╗####
####💠 Project Setup 💠####
####╚══════     ══════╝####

load_package_if_installed <- function(pkg) {
  suppressPackageStartupMessages({ require(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE) })
}

install_if_missing_then_load <- function(pkgs) {
  for (pkg in pkgs) {
    if (!load_package_if_installed(pkg)) 
      install.packages(pkg, verbose = FALSE, prompt = FALSE)
    
    load_package_if_installed(pkg)
  }
}

load_yml <- function(fp) {
  if (file.exists(fp)) 
    return(yaml::read_yaml(fp, eval.expr = TRUE))
  else {
    cli::cli_alert_warning("[SETUP] YAML file {.path {fp}} not found !")
    return(list())
  }
}

#-------------------------------#
####🔺Setting up the project ####
#-------------------------------#

install_if_missing_then_load(c("here", "renv"))

if (is.null(renv::project())) 
  renv::init(project = here(), bare = TRUE, restart = FALSE)

install_if_missing_then_load(c("rlang", "cli", "yaml"))

#--------------------------------#
####🔺Installing the packages ####
#--------------------------------#

cli_h1("[PACKAGES] Setting up the project packages")

source(here("src", "pakman.R"), echo = FALSE)

if (file.exists(here::here("renv.lock"))) {
  restore_project_packages()
} else {
  install_project_packages()
}

#--------------------------------#
####🔺Configuring the project ####
#--------------------------------#

cli_h1("[CONFIG] Setting project and packages' options")

source(here("src", "config.R"), echo = FALSE)
source(here("src", "theme.R"), echo = FALSE)

#---------------------#
####🔺Source files ####
#---------------------#

cli_h1("[SCRIPTS] Loading project scripts")

scripts_dir <- here("src", "scripts")

list.files(scripts_dir, pattern = "*.R") |> 
  grep(pattern = 'setup|init|config|theme|pakman', invert = TRUE, value = TRUE) |> 
  sapply(\(file) source(here(scripts_dir, file), echo = FALSE)) |> 
  invisible()
