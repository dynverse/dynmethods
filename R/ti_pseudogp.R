#' Inferring trajectories with pseudogp
#'
#' @inherit ti_angle description
#'
#' @inheritParams pseudogp::fitPseudotime
#' @param dimreds A character vector specifying which dimensionality reduction method to use.
#'   See [dyndimred::dimred] for the list of available dimensionality reduction methods.
#'
#' @seealso [pseudogp::fitPseudotime()]
#'
#' @export
ti_pseudogp <- create_ti_method(
  name = "pseudogp",
  short_name = "pseudogp",
  package_loaded = c("pseudogp"),
  package_required = c("rstan", "coda", "MCMCglmm", "dyndimred"),
  doi = "10.1371/journal.pcbi.1005212",
  trajectory_types = "linear",
  topology_inference = "fixed",
  type = "algorithm",
  license = "MIT",
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
      email = "cyau@well.ox.ac.uk",
      ORCID = "0000-0001-7615-8523"
    )
  ),
  preprint_date = "2016-04-05",
  publication_date = "2016-11-21",
  version = "0.1",
  code_url = "https://github.com/kieranrcampbell/pseudogp",
  parameters = list(
    smoothing_alpha = list(
      type = "numeric",
      default = 10,
      upper = 20,
      lower = 1,
      description = "The hyperparameter for the Gamma distribution that controls arc-length"
    ),
    smoothing_beta = list(
      type = "numeric",
      default = 3,
      upper = 20,
      lower = 1,
      description = "The hyperparameter for the Gamma distribution that controls arc-length"
    ),
    pseudotime_mean = list(
      type = "numeric",
      default = 0.5,
      upper = 1,
      lower = 0,
      description = "The mean of the constrained normal prior on the pseudotimes"
    ),
    pseudotime_var = list(
      type = "numeric",
      default = 1,
      upper = 1,
      lower = 0.01,
      description = "The variance of the constrained normal prior on the pseudotimes"
    ),
    chains = list(
      type = "integer",
      default = 3L,
      upper = 20L,
      lower = 1L,
      description = "The number of chains for the MCMC trace"
    ),
    iter = list(
      type = "numeric",
      default = 100,
      upper = 1000,
      lower = 100,
      description = "The number of iterations for the MCMC trace"
    ),
    dimreds = list(
      type = "logical_vector",
      default = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
      length = 9,
      description = "A character vector specifying which dimensionality reduction methods to use.\nSee \\code{\\link[dyndimred:dimred]{dyndimred::dimred}} for the list of available dimensionality reduction methods."
    ),
    initialise_from = list(
      type = "discrete",
      default = "random",
      values = c("random", "principal_curve", "pca"),
      description = "How to initialise the MCMC chain. One of \"random\" (stan decides),\n\"principal_curve\", or \"pca\" (the first component of PCA rescaled is taken to be the pseudotimes).\nNote: if multiple representations are provided, \\code{pseudogp} will take the principal curve or\npca from the first rather than combining them. If a particular representation is required, it is\nup to the user to re-order them."
    )
  ),
  run_fun = "dynmethods::run_pseudogp",
  plot_fun = "dynmethods::plot_pseudogp"
)

run_pseudogp <- function(
  expression,
  dimreds = names(dyndimred::list_dimred_methods()) == "pca",
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
  dimred_names <- names(dyndimred::list_dimred_methods())[as.logical(dimreds)]
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
  pseudotime <- MCMCglmm::posterior.mode(tmcmc) %>%
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
    pseudotime = pseudotime,
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
    geom_point(aes(comp_1, comp_2), alpha = .5, size = 3, colour = "white", fill = "darkred", shape = 21, dimreds) +
    geom_path(aes(X1, X2, group = paste0(space, "_", chain), alpha = alpha), curves, size = 1.5) +
    geom_text(aes(offset_x, offset_y+.45, label = name), space_ann, size = 8) +
    scale_alpha_identity()

  process_dynplot(g, prediction$id)
}
