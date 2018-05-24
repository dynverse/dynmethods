#' Inferring trajectories with PhenoPath
#'
#' @inherit ti_angle description
#'
#' @export
#'
#' @include wrapper_create_ti_method.R
ti_phenopth <- create_ti_method(
  name = "PhenoPath",
  short_name = "phenopth",
  package_required = c("phenopath"),
  package_loaded = c(),
  par_set = makeParamSet(
    makeIntegerParam(id = "thin", lower = 2L, upper = 500L, default = 40L),
    makeDiscreteParam(id = "z_init", default = 1, values = list(1, 2, 3, 4, 5, "random")),
    makeLogicalParam(id = "model_mu", default = FALSE),
    makeLogicalParam(id = "scale_y", default = TRUE)
  ),
  run_fun = "run_phenopath",
  plot_fun = "plot_phenopath"
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
    z_init = z_init,
    model_mu = model_mu,
    scale_y = scale_y
  )
  pseudotimes <- phenopath::trajectory(fit) %>%
    setNames(rownames(expression))

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # run pca for visualisation purposes
  space <- dimred(expression, method = "pca", ndim = 2)

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotimes = pseudotimes
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
    mutate(pseudotime = prediction$pseudotimes[cell_id])
  g <- ggplot(space) +
    geom_point(aes(Comp1, Comp2, colour = pseudotime)) +
    viridis::scale_colour_viridis() +
    labs(colour = "Pseudotime") +
    theme(legend.position = c(.92, .12))
  process_dynplot(g, prediction$id)
}
