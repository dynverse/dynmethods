#' Description for Ouijia
#' @export
description_ouija <- function() create_description(
  name = "ouija",
  short_name = "ouija",
  package_required = c("ouija", "rstan"),
  package_loaded = c("coda"),
  par_set = makeParamSet(
    makeNumericParam(id = "iter", lower = log(2), default = log(100), upper = log(50000), trafo = function(x) round(exp(x))), # default 10000
    makeDiscreteParam(id = "response_type", default = "switch", values = c("switch", "transient")),
    makeDiscreteParam(id = "inference_type", default = "hmc", values = c("hmc", "vb")),
    makeLogicalParam(id = "normalise_expression", default = TRUE)
  ),
  properties = c(),
  run_fun = run_ouija,
  plot_fun = plot_ouija
)

run_ouija <- function(
    expression,
    marker_feature_ids,
    iter = 1000, # default is actually 10'000.
    response_type = "switch",
    inference_type = "hmc",
    normalise_expression = TRUE
  ) {
  requireNamespace("ouija")
  requireNamespace("rstan")
  requireNamespace("coda")

  # ouija assumes that a small number of marker genes is used, otherwise the method is verrry slow
  expression <- expression[, marker_feature_ids]

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

  # obtain the pseudotimes
  pseudotimes <- ouija::map_pseudotime(oui)

  # run pca for visualisation purposes
  space <- stats::prcomp(expression)$x[,1:2] %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    as_data_frame() %>%
    mutate(pseudotime = pseudotimes)

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

  # TIMING: after postproc
  tl <- tl %>% add_timing_checkpoint("method_afterpostproc")

  # return output
  wrap_prediction_model_linear(
    cell_ids = rownames(expression),
    pseudotimes = pseudotimes,
    t0_df = t0_df,
    space = space
  ) %>% attach_timings_attribute(tl)
}

plot_ouija <- function(prediction) {
  g <- ggplot(prediction$space) +
    geom_point(aes(PC1, PC2, colour = pseudotime)) +
    viridis::scale_colour_viridis() +
    labs(colour = "Pseudotime") +
    theme(legend.position = c(.92, .12))
  process_dynplot(g, prediction$id)
}






#' Description for Ouijiaflow
#' @export
description_ouijaflow <- function() create_description(
  name = "ouijaflow",
  short_name = "ouijaf",
  package_required = c("ouija", "rstan"),
  package_loaded = c("coda"),
  par_set = makeParamSet(
    makeNumericParam(id = "iter", lower = log(2), default = log(100), upper = log(50000), trafo = function(x) round(exp(x))), # default 10000
    makeDiscreteParam(id = "response_type", default = "switch", values = c("switch", "transient")),
    makeDiscreteParam(id = "inference_type", default = "hmc", values = c("hmc", "vb")),
    makeLogicalParam(id = "normalise_expression", default = TRUE)
  ),
  properties = c(),
  run_fun = run_ouija,
  plot_fun = plot_ouija
)
