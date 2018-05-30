#' Wrappers for trajectory inference methods
#'
#' @import dplyr
#' @import tidyr
#' @import methods
#' @import ParamHelpers
#' @import tibble
#' @import ggplot2
#' @import dynwrap
#' @import dynutils
#' @importFrom dynplot process_dynplot plot_default
#' @importFrom stats cor dist kmeans median prcomp quantile runif setNames step time
#' @importFrom utils installed.packages head tail
#' @importFrom purrr %>% map map_df map_chr map_lgl map_int map_dbl keep discard set_names pmap map2
#' @importFrom magrittr %<>% %$% set_colnames set_rownames
#' @importFrom lhs randomLHS
# lhs is needed for ParamHelpers generateDesign
#'
#' @docType package
#' @name dynmethods
NULL
