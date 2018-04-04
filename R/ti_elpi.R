#' #' Description for elpi
#' #' @export
#' description_elpi <- function() abstract_aga_description("elpi")
#'
#' #' Description for elpifix
#' #' @export
#' description_elpifix <- function() abstract_aga_description("elpifix")
#'
#' abstract_elpi_description <- function(method) {
#'   par_set <- makeParamSet(
#'     makeNumericParam("Lambda", lower = 0.000001, upper = 2, default = )
#'   )
#'
#'   run_fun <- switch(
#'     method,
#'     aga = run_elpi,
#'     agapt = run_elpifix
#'   )
#'
#'   create_description(
#'     name = ifelse(method == "elpi", "ElPiGraph", "ElPiGraph fixed"),
#'     short_name = method,
#'     package_loaded = c(),
#'     package_required = c("ElPiGraph.R"),
#'     par_set = par_set,
#'     properties = c(),
#'     run_fun = run_fun,
#'     plot_fun = plot_elpi
#'   )
#' }
#'
#'
#' run_elpi <- function(
#'   expression,
#'   milestone_network = NULL,
#'   Lambda=0.01
#' ) {
#'   milestone_ids <- unique(c(milestone_network$from, milestone_network$to))
#'
#'   NumNodes <- length(task$milestone_ids)
#'   NumEdges <- length(task$milestone_network)
#'   InitEdges <- task$milestone_network %>%
#'     igraph::graph_from_data_frame() %>%
#'     igraph::get.data.frame() %>%
#'     mutate(from=as.numeric(factor(from, task$milestone_ids)), to=as.numeric(factor(to, task$milestone_ids))) %>%
#'     as.matrix()
#'
#'   remove_names <- function(x) {rownames(x) <- NULL;}
#'
#'   result <- computeElasticPrincipalTree(
#'     task$counts,
#'     NumNodes,
#'     NumEdges,
#'     InitNodePositions = task$counts[sample(rownames(task$counts), NumNodes),] %>% magrittr::set_rownames(NULL) %>% magrittr::set_colnames(NULL),
#'     InitEdges=InitEdges[, c(1, 2)],
#'     Lambda = Lambda
#'   )
#'
#' }
#'
#'
#' run_elpifix <- function(expression, milestone_network) {
#'   invoke(run_elpi, environment())
#' }
#' formals(run_elpifix) <- c(
#'   formals(run_elpifix),
#'   formals(run_elpi)[!names(formals(run_elpifix)) %in% names(formals(run_elpi))]
#' )
#'
