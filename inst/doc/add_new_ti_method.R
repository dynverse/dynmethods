## ---- message=FALSE------------------------------------------------------
library(dynwrap)
library(dynmethods)
library(tidyverse)
library(ParamHelpers)

## ------------------------------------------------------------------------
ncells <- 1000
pseudotime <- runif(ncells)

expression <- matrix(
  c(
    (pseudotime - 0.5) ** 2,
    sqrt(pseudotime + 20),
    pseudotime
  ),
  ncol = 3,
  dimnames = list(as.character(rep(seq_len(ncells))), as.character(c("A", "B", "C")))
)
expression <- expression + rnorm(length(expression), sd = 0.02)

## ------------------------------------------------------------------------
pca <- prcomp(expression)

## ------------------------------------------------------------------------
par_set <- makeParamSet(
  makeIntegerParam("component", lower = 1, upper = 10, default = 1)
)

## ------------------------------------------------------------------------
component <- 1
pseudotime <- pca$x[, component]
pseudotime <- (pseudotime - min(pseudotime)) / (max(pseudotime) - min(pseudotime))

## ---- echo=FALSE---------------------------------------------------------
data("priors", package = "dynwrap")
priors %>% 
  select(prior_id, prior_description) %>% 
  rename(`Name of parameter` = prior_id, `Description` = prior_description) %>% 
  knitr::kable("markdown")

## ------------------------------------------------------------------------
start_cell_ids <- as.character(which.min(pseudotime))

## ------------------------------------------------------------------------
if (!is.null(start_cell_ids)) {
  if(mean(pseudotime[start_cell_ids]) > 0.5) {
    pseudotime <- 1-pseudotime
  }
}

## ------------------------------------------------------------------------
milestone_ids <- c("A", "B")
milestone_network <- tibble(from = "A", to = "B", length = 1, directed = TRUE)
divergence_regions <- tibble()

## ------------------------------------------------------------------------
milestone_percentages <- bind_rows(
  tibble(
    milestone_id = "A",
    cell_id = names(pseudotime),
    percentage = 1-pseudotime
  ),
  tibble(
    milestone_id = "B",
    cell_id = names(pseudotime),
    percentage = pseudotime
  )
)

## ------------------------------------------------------------------------
progressions <- tibble(
  cell_id = names(pseudotime),
  from = "A",
  to = "B",
  percentage = pseudotime
)

## ------------------------------------------------------------------------
trajectory <- 
  wrap_data(
    cell_id = names(pseudotime)
  ) %>% 
  add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    divergence_regions = divergence_regions,
    progressions = progressions # either milestone_percentages or progressions have to be provided
  )

## ------------------------------------------------------------------------
trajectory <- 
  wrap_data(
    cell_id = names(pseudotime)
  ) %>% 
  add_linear_trajectory(pseudotime)

## ------------------------------------------------------------------------
dynplot::plot_dimred(trajectory, "pseudotime", expression_source = expression)

## ------------------------------------------------------------------------
run_fun <- function(expression, component, start_cell_id = NULL) {
  pca <- prcomp(expression)
  
  pseudotimes <- pca$x[, component]
  
  milestone_ids <- c("A", "B")
  milestone_network <- tibble(from="A", to="B", length=1, directed=TRUE)
  divergence_regions <- tibble()
  
  progressions <- tibble(
    cell_id = names(pseudotime),
    from = "A",
    to = "B",
    percentage = pseudotime
  )
  
  trajectory <- 
    wrap_data(
      cell_id = names(pseudotime)
    ) %>% 
    add_trajectory(
      milestone_ids = milestone_ids,
      milestone_network = milestone_network,
      divergence_regions = divergence_regions,
      progressions = progressions # either milestone_percentages or progressions have to be provided
    )
}

## ------------------------------------------------------------------------
ti_dummy <- create_ti_method(
  "dummy", 
  par_set,
  run_fun
)

## ------------------------------------------------------------------------
task <- wrap_data("", rownames(expression)) %>% add_expression(expression, expression)

model <- infer_trajectory(task, ti_dummy())

## ------------------------------------------------------------------------
library(dynplot)
dynplot::plot_dimred(model, color_cells = "pseudotime" , expression_source = task$expression)

