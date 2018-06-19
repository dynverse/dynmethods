# here we generate an overview of all methods

library(tidyverse)
library(googlesheets)
devtools::load_all()

load("data/methods_containerised.rda")

# process all code of method wrappers
methods_processed <- list.files("R", pattern = "ti_", full.names = T) %>%
  map_df(function(file) {
    file_text <- readr::read_lines(file)

    line_numbers <- which(str_detect(file_text, "^ti_[ a-zA-Z_0-9]*<-"))
    descr_funs <- str_replace(file_text[line_numbers], "(ti_[a-zA-Z_0-9]*) <-.*", "\\1")

    method_ids <- str_replace(file_text[line_numbers], "ti_([a-zA-Z_0-9]*).*", "\\1")

    data_frame(
      r_wrapper_location = file,
      line_number = line_numbers,
      fun_name = descr_funs,
      method_id = method_ids
    )
  }) %>%
  group_by(method_id) %>%
  filter(n() == 1 | r_wrapper_location != "R/ti_container.R") %>%
  mutate(r_wrapped = r_wrapper_location != "R/ti_container.R") %>%
  ungroup()

# now create every method and extract information
methods_info <- map(methods_processed$fun_name, function(fun_name) {
  method <- get(fun_name, asNamespace("dynmethods"))()
  method %>% discard(is.function)
}) %>% dynutils::list_as_tibble()

methods <-
  full_join(
    methods_processed,
    methods_containerised,
    "method_id"
  ) %>%
  left_join(
    methods_info,
    "method_id"
  )

if (any(is.na(methods$fun_name))) {
  stop("Some methods were containerised but the code location was not found!")
}

# determine wrapper location depending on docker
methods <- methods %>% mutate(
  containerised = ifelse(is.na(containerised), FALSE, containerised),
  r_wrapped = ifelse(is.na(r_wrapped), FALSE, r_wrapped),
  wrapper_location = case_when(
    r_wrapped ~ paste0(r_wrapper_location, "#L", line_number),
    containerised ~ docker_wrapper_location
  )
)

if (any(is.na(methods$wrapper_location))) {
  stop("Some wrapper locations were not found!")
}


usethis::use_data(methods, overwrite = TRUE)
