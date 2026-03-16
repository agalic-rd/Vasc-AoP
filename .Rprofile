if (file.exists("renv/activate.R")) {
    renv::activate()
}

source("renv/activate.R")

if (!"here" %in% as.data.frame(installed.packages())$Package) {
    renv::install("here", prompt = FALSE)
}
