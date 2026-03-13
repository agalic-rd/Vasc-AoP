source(here::here("src/renv/helpers.R"))

# See:
## - https://rstudio.github.io/renv/reference/config.html
## - https://rstudio.github.io/renv/reference/snapshot.html
options(
    repos = c(PPM = "https://packagemanager.posit.co/cran/latest", CRAN = "https://cloud.r-project.org"),
    renv.config.pak.enabled = FALSE,
    renv.config.sandbox.enabled = FALSE,
    renv.config.snapshot.inference = FALSE,
    renv.config.snapshot.validate = FALSE,
    renv.config.synchronized.check = FALSE
)

# GitHub
if (!nzchar(Sys.getenv("GITHUB_PAT"))) {
    warning("[RENV] GITHUB PAT not found - package loading might fail due to Github API's download cap.")
}
