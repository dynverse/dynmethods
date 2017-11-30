#' Run a method on a set of tasks with a set of parameters
#' @param tasks The tasks on which to evaluate.
#' @param method The method to evaluate.
#' @param parameters The parameters to evaluate with.
#' @param give_priors All the priors a method is allowed to receive. Must be a subset of: \code{"start_milestones"},
#'  \code{"start_cells"}, \code{"end_milestones"}, \code{"end_cells"}, \code{"grouping_assignment"} and \code{"grouping_network"}
#' @param timeout Kill execution after a given amount of time.
#' @param debug_timeout Setting debug to \code{TRUE} will avoid running the method in a separate R session
#'   using \code{\link[dynutils]{eval_with_timeout}} and run the method directly. Note that the timeout functionality
#'   will not work when \code{debug} is \code{TRUE}.
#'
#' @importFrom utils capture.output
#' @importFrom readr read_file
#' @importFrom stringr str_length
#' @export
execute_method <- function(
  tasks,
  method,
  parameters,
  give_priors = NULL,
  timeout = 3600 * 24 * 365.25 * 10,
  debug_timeout = FALSE
) {
  # Run method on each task
  lapply(seq_len(nrow(tasks)), function(i) {
    # start the timer
    time0 <- Sys.time()

    # get the task
    task <- dynutils::extract_row_to_list(tasks, i)

    # required params
    required_params <- formals(method$run_fun) %>%
      as.list %>%
      map_chr(class) %>%
      keep(~.=="name") %>%
      names

    # counts or expression arguments
    task_params <-
      if ("counts" %in% required_params) {
        list(counts = task$counts)
      } else if("expression" %in% required_params) {
        list(expression = task$expression)
      }

    # determine which prior information is strictly required by the method
    required_priors <- required_params
      discard(~.%in%c("counts", "expression"))

    # collect all prior information
    prior_information <- task$prior_information
    prior_information$task <- task

    # determine which priors to give and give it
    prior_names <- union(give_priors, required_priors)
    prior_type <- ifelse(prior_names %in% required_priors, "required", "optional")
    prior_list <- as.list(prior_information[prior_names])
    prior_df <- data_frame(prior_type, prior_names)

    # create arglist. content:
    # * counts or expression
    # * parameters
    # * any priors
    arglist <- c(task_params, parameters, prior_list)

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
        if (debug_timeout) {
          cat("Running ", method$name, " on ", task$id, " in debug mode!\n", sep = "")

          # Execute method
          model <- execute_method_internal(method, arglist, setseed_detection_file)
        } else {
          # Run the method on each of the tasks
          model <- dynutils::eval_with_timeout(
            timeout = timeout,
            expr = {
              tryCatch({
                # create a temporary directory to set as working directory,
                # to avoid polluting the working directory if a method starts
                # producing files :angry_face:
                setwd(tmp_dir)

                # run method
                execute_method_internal(method, arglist, setseed_detection_file)
              }, finally = {
                # return to old_wd
                # (in the likely event that this is a different thread)
                setwd(old_wd)
              })
            }
          )
        }

        c(model, list(error = NULL))
      }, error = function(e) {
        list(model = NULL, time1 = time0, time2 = Sys.time(), error = e)
      })

    # retrieve the model and error message
    model <- out$model
    error <- out$error
    time1 <- out$time1
    time2 <- out$time2

    # check whether the method produced output files and
    # wd to previous state
    num_files_created <- length(list.files(tmp_dir, recursive = TRUE))
    setwd(old_wd)

    # Temporary fix: do not remove the whole tmp wd, as futures might still be in this wd.
    # TODO: solve it
    # unlink(tmp_dir, recursive = TRUE, force = TRUE)
    for (x in list.dirs(tmp_dir, recursive = FALSE)) {
      unlink(x, recursive = TRUE, force = TRUE)
    }
    file.remove(list.files(tmp_dir))

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
      time_setup = as.numeric(difftime(time1, time0, units = "sec")),
      time_method = as.numeric(difftime(time2, time1, units = "sec")),
      time_cleanup = as.numeric(difftime(time3, time2, units = "sec")),
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
execute_method_internal <- function(method, arglist,setseed_detection_file) {
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
  time1 <- Sys.time()

  # execute method and return model
  model <- do.call(method$run_fun, arglist)

  # measure third time point
  time2 <- Sys.time()

  list(time1 = time1, time2 = time2, model = model)
}
