library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(SCOUP)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  expression,
  groups_id,
  start_id,
  end_n,
  ndim = 2,
  max_ite1 = 100,
  max_ite2 = 100,
  alpha_min = 0.1,
  alpha_max = 100,
  t_min = 0.001,
  t_max = 2,
  sigma_squared_min = 0.1,
  thresh = 0.01,
  verbose = FALSE
) {
  requireNamespace("SCOUP")

  # if the dataset is cyclic, pretend it isn't
  if (end_n == 0) {
    end_n <- 1
  }

  start_cell <- sample(start_id, 1)
  # figure out indices of starting population
  # from the groups_id and the start_cell
  start_ix <- groups_id %>%
    filter(cell_id %in% start_cell) %>%
    select(group_id) %>%
    left_join(groups_id, by = "group_id") %>%
    .$cell_id

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run SP and SCOUP
  model <- SCOUP::run_SCOUP(
    expr = expression,
    start_ix = start_ix,
    ndim = ndim,
    nbranch = end_n,
    max_ite1 = max_ite1,
    max_ite2 = max_ite2,
    alpha_min = alpha_min,
    alpha_max = alpha_max,
    t_min = t_min,
    t_max = t_max,
    sigma_squared_min = sigma_squared_min,
    thresh = thresh,
    verbose = verbose
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  if (any(is.na(model$ll))) {
    stop("SCOUP returned NaNs", call. = FALSE)
  }

  # create progressions
  milestone_percentages <- model$cpara %>%
    as.data.frame() %>%
    as_data_frame() %>%
    tibble::rownames_to_column("cell_id") %>%
    mutate(time = max(time) - time) %>%
    rename(M0 = time) %>%
    tidyr::gather(milestone_id, percentage, -cell_id) %>%
    group_by(milestone_id) %>%
    mutate(percentage = percentage / max(percentage)) %>%
    ungroup() %>%
    group_by(cell_id) %>%
    mutate(percentage = percentage / sum(percentage)) %>%
    filter(percentage > 0 | milestone_id == "M0")

  # create milestone ids
  milestone_ids <- c("M0", paste0("M", seq_len(end_n)))

  # create milestone network
  milestone_network <- data_frame(
    from = milestone_ids[[1]],
    to = milestone_ids[-1],
    length = 1,
    directed = TRUE
  )

  # create divergence_regions
  divergence_regions <- data_frame(
    divergence_id = "divergence",
    milestone_id = milestone_ids,
    is_start = milestone_id == milestone_ids[[1]]
  )

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression),
    cell_info = model$cpara %>% tibble::rownames_to_column("cell_id")
  ) %>% add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    divergence_regions = divergence_regions
  ) %>% add_dimred(
    dimred = model$dimred %>% as.matrix
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')