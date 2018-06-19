#' Inferring trajectories with mfa
#'
#' @inherit ti_angle description
#'
#' @inheritParams mfa::mfa
#'
#' @seealso [mfa::mfa()]
#'
#' @export
ti_mfa <- create_ti_method(
  name = "mfa",
  short_name = "mfa",
  package_loaded = c(),
  package_required = c("mfa", "dyndimred"),
  doi = "10.12688/wellcomeopenres.11087.1",
  trajectory_types = c("linear", "bifurcation", "convergence", "multifurcation", "binary_tree", "tree"),
  topology_inference = "parameter",
  type = "algorithm",
  license = "GPL (>= 2)",
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
  publication_date = "2017-03-15",
  version = "0.99.2",
  code_url = "https://github.com/kieranrcampbell/mfa",
  parameters = list(
    b = list(
      type = "integer",
      default = 2L,
      upper = 10L,
      lower = 1L,
      description = "Number of branches to model"
    ),
    iter = list(
      type = "integer",
      default = 2000L,
      upper = 5000L,
      lower = 20L,
      description = "Number of MCMC iterations"
    ),
    thin = list(
      type = "integer",
      default = 1L,
      upper = 20L,
      lower = 1L,
      description = "MCMC samples to thin"
    ),
    pc_initialise = list(
      type = "integer",
      default = 1L,
      upper = 5L,
      lower = 1L,
      description = "Which principal component to initialise pseudotimes to"
    ),
    prop_collapse = list(
      type = "numeric",
      default = 0,
      upper = 1,
      lower = 0,
      description = "Proportion of Gibbs samples which should marginalise over c"
    ),
    scale_input = list(
      type = "logical",
      default = TRUE,
      description = "Logical. If true, input is scaled to have mean 0 variance 1"
    ),
    zero_inflation = list(
      type = "logical",
      default = FALSE,
      description = "Logical, should zero inflation be enabled?"
    )
  ),
  run_fun = "dynmethods::run_mfa",
  plot_fun = "dynmethods::plot_mfa"
)

run_mfa <- function(
  expression,
  b = 2,
  iter = 2000,
  thin = 1,
  zero_inflation = FALSE,
  pc_initialise = 1,
  prop_collapse = 0,
  scale_input = TRUE
) {
  requireNamespace("mfa")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # perform MFA
  m <- mfa::mfa(
    y = expression,
    b = b,
    iter = iter,
    thin = thin,
    zero_inflation = zero_inflation,
    pc_initialise = pc_initialise,
    prop_collapse = prop_collapse,
    scale_input = scale_input
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # obtain results
  ms <- summary(m) %>%
    mutate(cell_id = rownames(expression)) %>%
    select(cell_id, everything())

  # create milestone network
  milestone_ids <- paste0("M", seq(0, b))
  milestone_network <- data_frame(
    from = "M0",
    to = paste0("M", seq_len(b)),
    length = 1,
    directed = TRUE
  )

  # create progressions
  progressions <- with(ms, data_frame(
    cell_id = rownames(expression),
    from = "M0",
    to = paste0("M", branch),
    percentage = dynutils::scale_minmax(pseudotime)
  ))

  # pca for visualisation only
  pca_out <- dyndimred::dimred(expression, method = "pca", ndim = 2)

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression),
    cell_info = ms
  ) %>% add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    divergence_regions = NULL
  ) %>% add_dimred(
    dimred = pca_out
  ) %>% add_timings(
    tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

plot_mfa <- function(prediction) {
  df <- data.frame(
    prediction$dimred,
    prediction$cell_info
  )
  g <- ggplot() +
    geom_point(aes(comp_1, comp_2, colour = branch), df) +
    labs(colour = "Branch") +
    theme(legend.position = c(.92, .1))
  process_dynplot(g, prediction$id)
}
