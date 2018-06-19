#' Inferring trajectories with PhenoPath
#'
#' @inherit ti_angle description
#'
#' @inheritParams phenopath::phenopath
#' @inheritParams phenopath::clvm
#'
#' @seealso [phenopath::phenopath()], [phenopath::clvm()]
#'
#' @export
ti_phenopath <- create_ti_method(
  name = "PhenoPath",
  short_name = "phenopath",
  package_required = c("phenopath", "dyndimred"),
  package_loaded = c(),
  doi = "10.1101/159913",
  trajectory_types = c("linear", "bifurcation", "convergence", "multifurcation"),
  topology_inference = "fixed",
  type = "algorithm",
  authors = list(
    list(
      given = "Kieran",
      family = "Campbell",
      email = "kicampbell@bccrc.ca",
      github = "kieranrcampbell"
    ),
    list(
      given = "Christopher",
      family = "Yau",
      email = "cyau@well.ox.ac.uk"
    )
  ),
  preprint_date = "2017-07-06",
  version = "1.1.1",
  code_url = "https://github.com/kieranrcampbell/phenopath",
  parameters = list(
    thin = list(
      type = "integer",
      default = 40L,
      upper = 500L,
      lower = 2L,
      description = "The number of iterations to wait each time before\nre-calculating the elbo"
    ),
    z_init = list(
      type = "discrete",
      default = "1",
      values = c("1", "2", "3", "4", "5", "random"),
      description = "The initialisation of the latent trajectory. Should be one of\n\\enumerate{\n\\item A positive integer describing which principal component of the data should\nbe used for initialisation (default 1), \\emph{or}\n\\item A numeric vector of length number of samples to be used \ndirectly for initialisation, \\emph{or}\n\\item The text character \\code{\"random\"}, for random initialisation \nfrom a standard normal distribution.\n}"
    ),
    model_mu = list(
      type = "logical",
      default = FALSE,
      description = "Logical - should a gene-specific intercept term be modelled?"
    ),
    scale_y = list(
      type = "logical",
      default = TRUE,
      description = "Logical - should the expression matrix be centre scaled?")
  ),
  run_fun = "dynmethods::run_phenopath",
  plot_fun = "dynmethods::plot_phenopath"
)

run_phenopath <- function(
  expression,
  thin = 40,
  z_init = 1,
  model_mu = FALSE,
  scale_y = TRUE
) {
  requireNamespace("phenopath")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run phenopath
  fit <- phenopath::phenopath(
    exprs_obj = expression,
    x = rep(1, nrow(expression)),
    elbo_tol = 1e-6,
    thin = thin,
    z_init = ifelse(z_init == "random", "random", as.numeric(z_init)),
    model_mu = model_mu,
    scale_y = scale_y
  )
  pseudotime <- phenopath::trajectory(fit) %>%
    setNames(rownames(expression))

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # run pca for visualisation purposes
  space <- dyndimred::dimred(expression, method = "pca", ndim = 2)

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotime = pseudotime
  ) %>% add_dimred(
    dimred = space
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom viridis scale_colour_viridis
plot_phenopath <- function(prediction) {
  space <- prediction$dimred %>%
    as.data.frame() %>%
    rownames_to_column("cell_id") %>%
    mutate(pseudotime = prediction$pseudotime[cell_id])
  g <- ggplot(space) +
    geom_point(aes(comp_1, comp_2, colour = pseudotime)) +
    viridis::scale_colour_viridis() +
    labs(colour = "Pseudotime") +
    theme(legend.position = c(.92, .12))
  process_dynplot(g, prediction$id)
}
