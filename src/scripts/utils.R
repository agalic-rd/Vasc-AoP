####╔═══    ══╗####
####💠 Utils 💠####
####╚═══    ══╝####

cli_h2("┗ [SCRIPTS] Loading helper functions")

#--------------#
####🔺Pipes ####
#--------------#

"%ni%" <- Negate("%in%")

"%s+%" <- \(lhs, rhs) paste0(lhs, rhs)

"%ne%" <- \(lhs, rhs) if(is.null(lhs) || rlang::is_empty(lhs) || (length(lhs) == 1 && lhs == "")) return(rhs) else return(lhs)

#-------------#
####🔺Misc ####
#-------------#

## Get element by name from a list:
rmatch <- function(x, name) {
  pos <- match(name, names(x))
  if (!is.na(pos)) return(x[[pos]])
  for (el in x) {
    if (class(el) == "list") {
      out <- Recall(el, name)
      if (!is.null(out)) return(out)
    }
  }
}