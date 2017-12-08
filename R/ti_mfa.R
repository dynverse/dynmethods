#' Description for mfa
#' @export
description_mfa <- function() create_description(
  name = "mfa",
  short_name = "mfa",
  package_loaded = c(),
  package_required = c("mfa"),
  par_set = makeParamSet(
    makeIntegerParam(id = "b", lower = 1L, upper = 10L, default = 2L),
    makeIntegerParam(id = "iter", lower = 20L, upper = 5000L, default = 2000L),
    makeIntegerParam(id = "thin", lower = 1L, upper = 20L, default = 1L),
    makeIntegerParam(id = "pc_initialise", lower = 1L, upper = 5L, default = 1L),
    makeNumericParam(id = "prop_collapse", lower = 0, upper = 1, default = 0),
    makeLogicalParam(id = "scale_input", default = TRUE),
    makeLogicalParam(id = "zero_inflation", default = FALSE)
  ),
  properties = c(),
  run_fun = run_mfa,
  plot_fun = plot_mfa
)

#' @importFrom stats prcomp
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

  # obtain results
  ms <- summary(m)

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
  pca_out <- stats::prcomp(expression)$x[,1:2]

  # return output
  wrap_prediction_model(
    trajectory_type = "multifurcating",
    cell_ids = rownames(expression),
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    ms = ms,
    pca_out = pca_out
  )
}

plot_mfa <- function(prediction) {
  df <- data.frame(
    prediction$pca_out,
    prediction$ms
  )
  g <- ggplot() +
    geom_point(aes(PC1, PC2, colour = branch), df) +
    labs(colour = "Branch") +
    theme(legend.position = c(.92, .1))
  process_dynplot(g, prediction$id)
}
