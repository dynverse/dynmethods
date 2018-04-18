set_cores <- function(n_cores = 1) {
  Sys.setenv(MKL_NUM_THREADS=n_cores, NUMEXPR_NUM_THREADS=n_cores, OMP_NUM_THREADS=n_cores)
}
