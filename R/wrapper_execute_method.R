#' Run a method on a set of tasks with a set of parameters
#' @param tasks The tasks on which to evaluate.
#' @param method The method to evaluate.
#' @param parameters The parameters to evaluate with.
#' @param give_priors All the priors a method is allowed to receive. Must be a subset of: \code{"start_milestones"},
#'  \code{"start_cells"}, \code{"end_milestones"}, \code{"end_cells"}, \code{"grouping_assignment"} and \code{"grouping_network"}
#' @param mc_cores The number of cores to use, allowing to parallellise the different tasks
#'
#' @importFrom utils capture.output
#' @importFrom readr read_file
#' @importFrom stringr str_length
#' @importFrom parallel mclapply
#' @export
execute_method <- function(
  tasks,
  method,
  parameters,
  give_priors = NULL,
  mc_cores = 1
) {
  # Run method on each task
  parallel::mclapply(seq_len(nrow(tasks)), mc.cores = mc_cores, function(i) {
    # start the timer
    time0 <- Sys.time()

    # get the task
    task <- dynutils::extract_row_to_list(tasks, i)

    # find out a method's required params (counts, expression, or any prior information)
    required_params <- formals(method$run_fun) %>%
      as.list %>%
      map_chr(class) %>%
      keep(~.=="name") %>%
      names

    # find out wether the method wants the counts, expression, or both
    task_args <- task[intersect(required_params, c("counts", "expression"))]

    # determine which prior information is strictly required by the method
    required_priors <- setdiff(required_params, c("counts", "expression"))

    # collect all prior information
    prior_information <- task$prior_information
    prior_information$task <- task

    # determine which priors to give and give it
    prior_names <- union(give_priors, required_priors)
    prior_type <- ifelse(prior_names %in% required_priors, "required", "optional")

    if (!all(prior_names %in% names(prior_information))) {
      stop("Prior information ", paste(setdiff(prior_names, names(prior_information)), collapse = ";"), " is missing from ", task$id)
    }

    prior_args <- as.list(prior_information[prior_names])
    prior_df <- data_frame(prior_type, prior_names)

    # create arglist. content:
    # * counts or expression
    # * parameters
    # * any prior information
    arglist <- c(task_args, parameters, prior_args)

    # create a temporary directory to set as working directory,
    # to avoid polluting the working directory if a method starts
    # producing files :angry_face:
    tmp_dir <- tempfile(pattern = method$short_name)
    dir.create(tmp_dir)
    old_wd <- getwd()
    setwd(tmp_dir)

    # disable seed setting
    # a method shouldn't set seeds during regular execution,
    # it should be left up to the user instead
    orig_setseed <- base::set.seed
    setseed_detection_file <- tempfile(pattern = "seedsetcheck")

    # run the method and catch the error, if necessary
    out <-
      tryCatch({
        # run method
        model <- execute_method_internal(method, arglist, setseed_detection_file)

        # add task id and method names to the model
        model$task_id <- task$id
        model$method_name <- method$name
        model$method_short_name <- method$short_name

        c(model, list(error = NULL))
      }, error = function(e) {
        time_new <- Sys.time()
        timings_list <- list(
          method_start = time0,
          method_afterpreproc = time0,
          method_aftermethod = time_new,
          method_afterpostproc = time_new,
          method_stop = time_new
        )
        list(model = NULL, timings_list = timings_list, error = e)
      })

    # retrieve the model, error message, and timings
    model <- out$model
    error <- out$error
    timings_list <- out$timings_list

    # check whether the method produced output files and
    # wd to previous state
    num_files_created <- length(list.files(tmp_dir, recursive = TRUE))
    setwd(old_wd)

    # Remove temporary folder
    unlink(tmp_dir, recursive = TRUE, force = TRUE)

    # read how many seeds were set and
    # restore environment to previous state
    num_setseed_calls <-
      if (file.exists(setseed_detection_file)) {
        stringr::str_length(readr::read_file(setseed_detection_file))
      } else {
        0
      }
    if (file.exists(setseed_detection_file)) {
      file.remove(setseed_detection_file)
    }
    dynutils::override_setseed(orig_setseed)

    # stop the timer
    time3 <- Sys.time()

    # create a summary tibble
    summary <- tibble(
      method_name = method$name,
      method_short_name = method$short_name,
      task_id = task$id,
      time_sessionsetup = as.numeric(difftime(timings_list$method_start, time0, units = "sec")),
      time_preprocessing = as.numeric(difftime(timings_list$method_afterpreproc, timings_list$method_start, units = "sec")),
      time_method = as.numeric(difftime(timings_list$method_aftermethod, timings_list$method_afterpreproc, units = "sec")),
      time_postprocessing = as.numeric(difftime(timings_list$method_afterpostproc, timings_list$method_aftermethod, units = "sec")),
      time_wrapping = as.numeric(difftime(timings_list$method_stop, timings_list$method_afterpostproc, units = "sec")),
      time_sessioncleanup = as.numeric(difftime(time3, timings_list$method_stop, units = "sec")),
      error = list(error),
      num_files_created = num_files_created,
      num_setseed_calls = num_setseed_calls,
      prior_df = list(prior_df)
    )

    # return output
    list(model = model, summary = summary)
  })
}

#' Internal method for executing a method
#'
#' If you're reading this, you're supposed to be using \code{execute_method} instead.
#'
#' @param method A TI method wrapper
#' @param arglist The arguments to apply to the method
#' @param setseed_detection_file A file to which will be written if a method
#'   uses the set.seed function.
#'
#' @export
#' @importFrom readr write_file
execute_method_internal <- function(method, arglist, setseed_detection_file) {
  # disable seed setting
  # a method shouldn't set seeds during regular execution,
  # it should be left up to the user instead
  new_setseed <- function(i) {
    readr::write_file("1", setseed_detection_file, append = TRUE)
  }
  dynutils::override_setseed(new_setseed)

  # Load required packages and namespaces
  for (pack in method$package_loaded) {
    suppressMessages(do.call(require, list(pack)))
  }
  for (pack in method$package_required) {
    suppressMessages(do.call(requireNamespace, list(pack)))
  }

  # measure second time point
  time_start <- Sys.time()

  # execute method and return model
  model <- do.call(method$run_fun, arglist)

  # measure third time point
  time_stop <- Sys.time()

  # fetch timings from within method (and place them in order of execution, just to make sure)
  timings_list <- c(
    list(method_start = time_start),
    get_timings_attribute(model),
    list(method_stop = time_stop)
  )

  # return output
  lst(timings_list, model)
}
