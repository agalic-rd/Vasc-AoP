get_r_version <- function() {
    return(paste0(version$major, ".", sub("\\..*", "", version$minor)))
}

is_installed <- function(pkg) {
    suppressMessages({
        require(pkg, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)
    })
}

read_packages <- function(packages_file = "packages.txt") {
    packages <- readLines(here::here(packages_file))
    packages <- trimws(packages)
    packages <- packages[packages != ""]
    packages <- gsub(packages, pattern = ",", replacement = "")
    return(packages)
}

get_pkg_name <- function(remotes_string) {
    res <- renv:::renv_remotes_parse(remotes_string)
    return(res$package %||% res$repo)
}

should_restore <- function(packages_file = "packages.txt") {
    void_ <- capture.output({
        status <- renv::status()
    })
    packages_installed <- names(status$library$Packages)

    packages_needed <- read_packages(packages_file)
    packages_needed <- sapply(packages_needed, get_pkg_name)

    if (length(setdiff(packages_needed, packages_installed)) > 0) {
        return(TRUE)
    }

    return(FALSE)
}

safe_restore <- function() {
    lockfile_path <- renv::paths$lockfile()

    if (file.exists(lockfile_path)) {
        renv::restore(prompt = FALSE, repos = getOption("repos"))
    }
    invisible(lockfile_path)
}

install_all_packages <- function(packages_file = "packages.txt") {
    packages <- read_packages(packages_file)

    renv::install(packages, prompt = FALSE, rebuild = TRUE, repos = getOption("repos"))
    renv::snapshot(packages = sapply(packages, get_pkg_name), prompt = FALSE, force = TRUE)
}

load_all_packages <- function(packages_file = "packages.txt") {
    packages <- read_packages(packages_file)
    package_names <- sapply(packages, get_pkg_name)

    void_ <- lapply(package_names, \(pkg_name) {
        suppressPackageStartupMessages({
            library(pkg_name, character.only = TRUE)
        })
    })
}
