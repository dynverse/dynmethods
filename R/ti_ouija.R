#' Inferring trajectories with ouija
#'
#' @inherit ti_angle description
#'
#' @inheritParams ouija::ouija
#' @param iter Number of iterations
#'
#' @seealso [ouija::ouija()]
#'
#' @export
ti_ouija <- create_ti_method(
  name = "ouija",
  short_name = "ouija",
  package_required = c("ouija", "rstan"),
  package_loaded = c("coda"),
  parameters = list(
    iter = list(
      type = "numeric",
      default = 100,
      upper = 1000,
      lower = 10,
      distribution = "exponential",
      rate = 0.01,
      description = "Number of iterations"
    ),
    response_type = list(
      type = "discrete",
      default = "switch",
      values = c("switch", "transient"),
      description = "A vector declaring whether each gene exhibits \"switch\" or \"transient\"\nexpression. Defaults to \"switch\" for all genes"
    ),
    inference_type = list(
      type = "discrete",
      default = "hmc",
      values = c("hmc", "vb"),
      description = "The type of inference to be performed, either \\code{hmc} for Hamiltonian\nMonte Carlo or \\code{vb} for ADVI (Variational Bayes). Note that HMC is typically more accurate\nbut VB will be orders of magnitude faster."
    ),
    normalise_expression = list(
      type = "logical",
      default = TRUE,
      values = c("TRUE", "FALSE"),
      description = "Logical, default TRUE. If TRUE the data is pre-normalised\nso the average peak expression is approximately 1. This makes the strength parameters\napproximately comparable between genes."
    )
  ),
  run_fun = "dynmethods::run_ouija",
  plot_fun = "dynmethods::plot_ouija"
)

run_ouija <- function(
    expression,
    features_id,
    iter = 1000, # default is actually 10'000.
    response_type = "switch",
    inference_type = "hmc",
    normalise_expression = TRUE
  ) {
  requireNamespace("ouija")
  requireNamespace("rstan")
  requireNamespace("coda")

  # ouija assumes that a small number of marker genes is used, otherwise the method is verrry slow
  expression <- expression[, features_id]

  # write compiled instance of the stanmodel to HDD
  rstan::rstan_options(auto_write = TRUE)

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run ouija
  oui <- ouija::ouija(
    x = expression,
    iter = iter,
    response_type = response_type,
    inference_type = inference_type,
    normalise_expression = normalise_expression
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # obtain the pseudotime
  pseudotime <- ouija::map_pseudotime(oui) %>%
    setNames(rownames(expression))

  # run pca for visualisation purposes
  space <- dimred(expression, method = "pca", ndim = 2)

  # extract data for visualisation
  # adapted from ouija::plot_switch_times(oui)
  # to avoid saving the whole oui object
  k_trace <- rstan::extract(oui$fit, "k")$k
  kmean <- colMeans(k_trace)
  t0 <- rstan::extract(oui$fit, "t0")$t0
  t0_means <- colMeans(t0)
  t0_interval <- coda::HPDinterval(coda::mcmc(t0))
  t0_df <- data_frame(t0_mean = t0_means, lower = t0_interval[, 1], upper = t0_interval[, 2], kmean = kmean)
  t0_df$Gene <- colnames(oui$Y[, oui$response_type == "switch"])
  t0_df$Gene <- factor(t0_df$Gene, t0_df$Gene[order(t0_means)])

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotime = pseudotime,
    t0_df = t0_df
  ) %>% add_dimred(
    dimred = space
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_ouija <- function(prediction) {
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
