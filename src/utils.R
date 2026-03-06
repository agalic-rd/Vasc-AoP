cli::cli_h2("┗ [Vasc-AoP] Loading utils")

#--------------#
####🔺Pipes ####
#--------------#

"%ni%" <- Negate("%in%")

"%s+%" <- function(lhs, rhs) paste0(lhs, rhs)

"%ne%" <- function(lhs, rhs) {
    if (is.null(lhs) || rlang::is_empty(lhs) || (length(lhs) == 1 && lhs == "")) {
        return(rhs)
    } else {
        return(lhs)
    }
}

'%|e|%' <- function(x, y) {
    if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

#-------------#
####🔺Misc ####
#-------------#

## Get element by name from a list:
rmatch <- function(x, name) {
    pos <- match(name, names(x))
    if (!is.na(pos)) {
        return(x[[pos]])
    }
    for (el in x) {
        if (inherits(el, "list")) {
            out <- Recall(el, name)
            if (!is.null(out)) return(out)
        }
    }
}

#-------------------------#
####🔺Formatting utils ####
#-------------------------#

stars_pval <- function(p) {
    gtools::stars.pval(p) |> stringr::str_replace(fixed("."), "")
}

label_pval <- function(p) {
    str_c(scales::label_pvalue()(p) |> str_remove(">") |> str_trim(), stars_pval(p), sep = " ")
}

get_response_name <- function(var, col = "Label", dict = data_dict) {
    res <- dict[[1]] |> filter(Name == var) |> pull(col)

    return(res %ne% var)
}

get_var_level_name <- function(var, level, col = "Description", dict = data_dict) {
    res <- dict[[var]] |> filter(Name == level) |> pull(col)

    return(res %ne% level)
}

#---------------#
####🔺Images ####
#---------------#

save_png <- function(
    plot,
    filename = NULL,
    subfolder = NULL,
    device = "png",
    dpi = 600,
    width = 8,
    height = 8,
    display = TRUE
) {
    if (is.null(filename)) {
        filename <- as.list(match.call()[-1])$plot
    }

    file_path <- here("fig", paste0(filename, ".", device))
    if (!is.null(subfolder)) {
        if (!dir.exists(here::here("fig", subfolder))) {
            dir.create(here::here("fig", subfolder), recursive = TRUE)
        }
        file_path <- here("fig", subfolder, paste0(filename, ".", device))
    }

    ggsave(filename = file_path, plot = plot, device = device, scale = 1, dpi = dpi, width = width, height = height)

    if (display) return(plot)
}
