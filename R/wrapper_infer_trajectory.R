#' Infer a trajectory on one (or more) tasks, using one (or more) methods
#' @param task A task, a tibble of tasks
#' @param method A method, list of methods, or a design tibble
#' @param parameters Parameters, or a list of parameters
#' @param tasks One or more tasks (either as a list or tibble)
#' @param methods One or more methods (either as a list or tibble)
#'
#' @export
infer_trajectory <- function(task, method, parameters = NULL) {
  infer_trajectories(task, method, parameters) %>% pull(model) %>% first()
}

#' @rdname infer_trajectory
#' @export
infer_trajectories <- function(tasks, methods, parameters = NULL) {
  # process method
  # allow giving names of methods
  if(is.character(methods)) {
    methods <- dynmethods::get_descriptions() %>% slice(map_int(methods, agrep, dynmethods::get_descriptions()$short_name))
    # allow single description
  } else if ("dynmethod::description" %in% class(methods)) {
    methods <- list_as_tibble(list(methods))
  } else if (is.data.frame(methods)) {
  } else if (is.list(methods)) {
    methods <- list_as_tibble(methods)
  } else {
    stop("Invalid methods")
  }
  methods <- map(seq_len(nrow(methods)), extract_row_to_list, tib=methods)

  # process task
  # allow single task
  if(dynwrap::is_data_wrapper(tasks)) {
    tasks <- list_as_tibble(list(tasks))
  } else if (is.data.frame(tasks)) {
  } else if (is.list(tasks)) {
    tasks <- list_as_tibble(tasks)
  } else {
    stop("Invalid task")
  }
  tasks <- map(seq_len(nrow(tasks)), extract_row_to_list, tib=tasks)

  # process parameters
  # allow null parameters
  if (is.null(parameters)) {
    parameters <- map(seq_along(methods), ~list())
  }

  # construct overall design
  design <- tibble(
    task = tasks,
    method = methods,
    parameters = parameters
  )

  # execute
  design$model <- design %>%
    pmap(execute_method_on_task)

  design
}
