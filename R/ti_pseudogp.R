#' Inferring trajectories with pseudogp
#'
#' @inherit ti_angle description
#'
#' @inheritParams pseudogp::fitPseudotime
#' @param dimreds A character vector specifying which dimensionality reduction methods to use.
#'   See \code{\link{list_dimred_methods}} for the list of available dimensionality reduction methods.
#'
#' @seealso [pseudogp::fitPseudotime()]
#'
#' @export
#'
#' @include wrapper_create_ti_method.R
ti_pseudogp <- create_ti_method(
  name = "pseudogp",
  short_name = "pseudogp",
  package_loaded = c("pseudogp"),
  package_required = c("rstan", "coda", "MCMCglmm"),
  par_set = makeParamSet(
    makeNumericParam(id = "smoothing_alpha", lower = 1, upper = 20, default = 10),
    makeNumericParam(id = "smoothing_beta", lower = 1, upper = 20, default = 3),
    makeNumericParam(id = "pseudotime_mean", lower = 0, upper = 1, default = 0.5),
    makeNumericParam(id = "pseudotime_var", lower = 0.01, upper = 1, default = 1),
    makeIntegerParam(id = "chains", lower = 1L, default = 3L, upper = 20L),
    makeNumericParam(id = "iter", lower = log(100), default = log(100), upper = log(1000), trafo = function(x) round(exp(x))), # default is 1000
    makeLogicalVectorParam(id = "dimreds", len = length(list_dimred_methods()), default = names(list_dimred_methods()) %in% c("pca", "mds")),
    makeDiscreteParam(id = "initialise_from", values = c("random", "principal_curve", "pca"), default = "random")
  ),
  run_fun = "run_pseudogp",
  plot_fun = "plot_pseudogp"
)

run_pseudogp <- function(
  expression,
  dimreds = names(list_dimred_methods()) == "pca",
  chains = 1,
  iter = 1000,
  smoothing_alpha = 10,
  smoothing_beta = 3,
  pseudotime_mean = 0.5,
  pseudotime_var = 1,
  initialise_from = "random"
) {
  requireNamespace("pseudogp")
  requireNamespace("rstan")
  requireNamespace("coda")
  requireNamespace("MCMCglmm")

  # perform dimreds
  dimred_names <- names(list_dimred_methods())[dimreds]
  spaces <- map(dimred_names, ~ dimred(expression, method = ., ndim = 2)) # only 2 dimensions per dimred are allowed

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # fit probabilistic pseudotime model
  fit <- pseudogp::fitPseudotime(
    X = spaces,
    smoothing_alpha = smoothing_alpha,
    smoothing_beta = smoothing_beta,
    iter = iter,
    chains = chains,
    initialise_from = initialise_from,
    pseudotime_var = pseudotime_var,
    pseudotime_mean = pseudotime_mean
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # extract pseudotime
  pst <- rstan::extract(fit, pars = "t")$t
  tmcmc <- coda::mcmc(pst)
  pseudotimes <- MCMCglmm::posterior.mode(tmcmc) %>%
    setNames(rownames(expression))

  # collect data for visualisation purposes
  # code is adapted from pseudogp::posteriorCurvePlot
  pst <- rstan::extract(fit, pars = "t", permute = FALSE)
  lambda <- rstan::extract(fit, pars = "lambda", permute = FALSE)
  sigma <- rstan::extract(fit, pars = "sigma", permute = FALSE)

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_linear_trajectory(
    pseudotimes = pseudotimes,
    spaces = spaces,
    chains = chains,
    pst = pst,
    lambda = lambda,
    sigma = sigma
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_pseudogp <- function(prediction) {
  # code is adapted from pseudogp::posteriorCurvePlot
  requireNamespace("pseudogp")
  requireNamespace("MCMCglmm")
  requireNamespace("coda")

  spaces <- prediction$spaces
  chains <- prediction$chains
  pst <- prediction$pst
  lambda <- prediction$lambda
  sigma <- prediction$sigma

  Ns <- length(spaces)
  cols <- ceiling(sqrt(Ns))
  rows <- ceiling(Ns/cols)

  offset_mult <- 1.1
  space_ann <- data_frame(
    i = seq_len(Ns),
    offset_x = (((i-1) %% cols) + 1) * offset_mult,
    offset_y = (floor((i - 1) / cols) + 1) * offset_mult,
    name = names(spaces)
  )

  plot_outs <- lapply(seq_len(Ns), function(i) {
    lams <- lambda[, , (2 * i - 1):(2 * i), drop = FALSE]
    sigs <- sigma[, , (2 * i - 1):(2 * i), drop = FALSE]
    x <- spaces[[i]]

    plot_offset <- c(space_ann$offset_x[[i]], space_ann$offset_y[[i]])

    xsc <- dynutils::scale_uniform(x) %>% sweep(2, plot_offset, "+")

    xsc_df <- data.frame(space = names(spaces)[[i]], xsc, stringsAsFactors = FALSE)

    ncurves <- min(50, nrow(x))
    n_posterior_samples <- dim(pst)[1]
    curve_samples <- sample(n_posterior_samples, ncurves)
    pmcs <- map_df(seq_len(chains), function(chain) {
      tmap <- MCMCglmm::posterior.mode(coda::mcmc(pst[, chain, ]))
      lmap <- MCMCglmm::posterior.mode(coda::mcmc(lams[, chain, ]))
      smap <- MCMCglmm::posterior.mode(coda::mcmc(sigs[, chain, ]))

      pmc <- pseudogp:::posterior_mean_curve(x, tmap, lmap, smap, nnt = 80)
      pmcsc <- apply_uniform_scale(pmc$mu, attr(xsc, "addend"), attr(xsc, "multiplier")) %>%
        sweep(2, plot_offset, "+")

      data.frame(space = names(spaces)[[i]], chain, pmcsc, t = pmc$t, alpha = .5 * exp(1) * exp(-ncurves) + .5, stringsAsFactors = FALSE)
    })

    lst(
      space = xsc_df,
      curves = pmcs
    )
  })

  dimreds <- plot_outs %>% map_df(~ .$space)
  curves <- plot_outs %>% map_df(~ .$curves) %>% arrange(space, chain, t)

  # calculated_alpha <- .5 * exp(1) * exp(-ncurves) + .5

  g <- ggplot() +
    geom_point(aes(Comp1, Comp2), alpha = .5, size = 3, colour = "white", fill = "darkred", shape = 21, dimreds) +
    geom_path(aes(X1, X2, group = paste0(space, "_", chain), alpha = alpha), curves, size = 1.5) +
    geom_text(aes(offset_x, offset_y+.45, label = name), space_ann, size = 8) +
    scale_alpha_identity()

  process_dynplot(g, prediction$id)
}
