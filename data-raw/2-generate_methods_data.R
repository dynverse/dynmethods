# here we generate an overview of all methods

library(tidyverse)
library(googlesheets)

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

# load from googlesheets some extra information on the methods
method_infos <- gs_key("1Mug0yz8BebzWt8cmEW306ie645SBh_tDHwjVw4OFhlE") %>%
  gs_read(ws = "Methods", skip = 1)

implementation_infos <- gs_key("1Mug0yz8BebzWt8cmEW306ie645SBh_tDHwjVw4OFhlE") %>%
  gs_read(ws = "Implementations", skip = 1)

method_infos <- left_join(method_infos, implementation_infos, "implementation_id")

if (!all(methods$method_id %in% method_infos$method_id)) {
  stop("Not all methods found in sheet: \n", setdiff(methods$method_id, method_infos$method_id))
}

methods <-
  left_join(
    methods_processed,
    method_infos %>% select(method_id, method_name, type, DOI, code_location),
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
