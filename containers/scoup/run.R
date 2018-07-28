library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#' @examples
#' data <- dyntoy::generate_dataset(unique_id = "test", num_cells = 300, num_genes = 300, model = "linear") %>% c(., .$prior_information)
#' params <- yaml::read_yaml("containers/scoup/definition.yml")$parameters %>%
#'   {.[names(.) != "forbidden"]} %>%
#'   map(~ .$default)

expression <- data$expression
groups_id <- data$groups_id
start_id <- data$start_id
end_n <- data$end_n

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

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

# create distribution on starting population
vars <- apply(expression[start_ix,, drop = F], 2, stats::var)
vars[vars == 0] <- diff(range(expression)) * 1e-3
means <- apply(expression[start_ix,, drop=F], 2, mean)
distr_df <- data.frame(i = seq_along(vars) - 1, means, vars)

# write data to files
utils::write.table(t(expression), file = "data", sep = "\t", row.names = FALSE, col.names = FALSE)
utils::write.table(distr_df, file = "init", sep = "\t", row.names = FALSE, col.names = FALSE)

# TIMING: done with preproc
checkpoints <- list(method_afterpreproc = as.numeric(Sys.time()))

# execute sp
cmd <- paste0("SCOUP/sp data init time_sp dimred ", ncol(expression), " ", nrow(expression), " ", params$ndim)
cat(cmd, "\n", sep = "")
system(cmd)

# execute scoup
cmd <- paste0(
  "SCOUP/scoup data init time_sp gpara cpara ll ",
  ncol(expression), " ", nrow(expression),
  " -k ", end_n,
  " -m ", params$max_ite1,
  " -M ", params$max_ite2,
  " -a ", params$alpha_min,
  " -A ", params$alpha_max,
  " -t ", params$t_min,
  " -T ", params$t_max,
  " -s ", params$sigma_squared_min,
  " -e ", params$thresh
)
cat(cmd, "\n", sep = "")
system(cmd)

# read dimred
dimred <- utils::read.table("dimred", col.names = c("i", paste0("Comp", seq_len(params$ndim))))

# last line is root node
root <- dimred[nrow(dimred),-1,drop=F]
dimred <- dimred[-nrow(dimred),-1]
rownames(dimred) <- rownames(expression)

# read cell params
cpara <- utils::read.table("cpara", col.names = c("time", paste0("M", seq_len(end_n))))
rownames(cpara) <- rownames(expression)

# loglik
ll <- utils::read.table("ll")[[1]]

# TIMING: done with method
checkpoints$method_aftermethod <- as.numeric(Sys.time())

if (any(is.na(ll))) {
  stop("SCOUP returned NaNs", call. = FALSE)
}

pseudotime <- cpara %>% {set_names(.$time, rownames(.))}
esp <- cpara %>% select(-time) %>% tibble::rownames_to_column("cell_id")

# return output
output <- lst(
  end_state_probabilities = esp,
  pseudotime,
  do_scale_minmax = TRUE,
  dimred = dimred %>% as.matrix,
  timings = checkpoints
)

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(output, "/output/output.rds")
