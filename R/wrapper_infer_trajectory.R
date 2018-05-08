#' Infer a trajectory on one (or more) tasks, using one (or more) methods
#' @param task A task, a tibble of tasks
#' @param method A method, list of methods, or a design tibble
#' @param parameters Parameters, or a list of parameters
infer_trajectory <- function(task, method, parameters = NULL) {
  # return one model as default
  mode <- "single"

  # process method
  # allow giving names of methods
  if(is.character(method)) {
    if(length(method) > 1) {
      mode <- "multiple"
    }
    method <- dynmethods::get_descriptions() %>% slice(match(method, short_name))
  # allow single description
  } else if ("dynmethod::description" %in% class(method)) {
    method <- list_as_tibble(list(method))
  } else if (is.data.frame(method)) {
    mode <- "multiple"
  } else if (is.list(method)) {
    mode <- "multiple"
    method <- list_as_tibble(method)
  } else {
    stop("Invalid method")
  }
  method <- map(seq_len(nrow(method)), extract_row_to_list, tib=method)

  # process task
  # allow single task
  if(dynwrap::is_data_wrapper(task)) {
    task <- list_as_tibble(list(task))
  } else if (is.data.frame(task)) {
    mode <- "multiple"
  } else if (is.list(task)) {
    task <- list_as_tibble(task)
    mode <- "multiple"
  } else {
    stop("Invalid task")
  }
  task <- map(seq_len(nrow(task)), extract_row_to_list, tib=task)

  # process parameters
  # allow null parameters
  if (is.null(parameters)) {
    parameters <- map(seq_along(method), ~list())
  }

  # construct overall design
  design <- tibble(
    task = task,
    method = method,
    parameters = parameters
  )

  # execute
  design$model <- design %>%
    pmap(execute_method_on_task)

  if (mode == "single") {
    design$model[[1]]
  } else {
    design
  }
}
