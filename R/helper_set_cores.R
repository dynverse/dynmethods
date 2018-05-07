#' Helper function for nasty Python code
#'
#' @param n_codes the number of cores Python is allowed to use.
set_cores <- function(n_cores = 1) {
  Sys.setenv(
    MKL_NUM_THREADS = n_cores,
    NUMEXPR_NUM_THREADS = n_cores,
    OMP_NUM_THREADS = n_cores
  )
}
