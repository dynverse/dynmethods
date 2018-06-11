# here we generate an overview of all methods

library(tidyverse)
library(googlesheets)

load("data/methods_info.rda")
load("data/methods_containerised.rda")

# process all method wrappers
methods_processed <- list.files("R", pattern = "ti_", full.names = T) %>%
  map_df(function(file) {
    file_text <- readr::read_lines(file)

    line_numbers <- which(str_detect(file_text, "^ti_[ a-zA-Z]*<-"))
    descr_funs <- str_replace(file_text[line_numbers], "(ti_[a-zA-Z]*) <-.*", "\\1")

    method_ids <- str_replace(file_text[line_numbers], "ti_([a-zA-Z]*).*", "\\1")

    data_frame(
      r_wrapper_location = file,
      line_number = line_numbers,
      fun_name = descr_funs,
      method_id = method_ids
    )
  })

if (!all(methods_processed$method_id %in% methods_info$method_id)) {
  stop("Not all methods found in sheet: \n", setdiff(methods_processed$method_id, methods_info$method_id))
}

methods <-
  left_join(
    methods_processed,
    methods_info %>% select(method_id, method_name, type, DOI, code_location),
    "method_id"
  ) %>%
  left_join(
    methods_containerised,
    "method_id"
  )

# determine wrapper location depending on docker
methods <- methods %>% mutate(
  containerised = ifelse(is.na(containerised), FALSE, containerised),
  wrapper_location = case_when(
    !is.na(docker_wrapper_location) ~ docker_wrapper_location,
    !is.na(r_wrapper_location) ~ paste0(r_wrapper_location, "#", line_number)
  )
)


usethis::use_data(methods, overwrite = TRUE)
