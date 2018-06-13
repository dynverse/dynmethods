# here we do some metaprogramming to generate the ti_{method} functions inside ti_{method}.R

library(tidyverse)
library(dynwrap)
library(googlesheets)

load("data/methods_info.rda")

#   ____________________________________________________________________________
#   Run & push dockers                                                      ####
# get all methods
method_ids <- list.dirs("./containers/", full.names = FALSE)[-1]

if (!all(method_ids %in% methods_info$method_id)) {
  stop("Some methods not in methods_info! \n", setdiff(method_ids, methods_info$method_id))
}

# rebuild & push all containers
rebuild <- FALSE
if (rebuild) {
  walk(method_ids, function(method_id) {
    system(str_glue("docker build containers/{method_id} -t dynverse/{method_id}"))
    system(str_glue("docker push dynverse/{method_id}"))
  })
}

#   ____________________________________________________________________________
#   Functions to generate documentation and parameters given a definition   ####
generate_documentation_from_definition <- function(method_id, definition) {
  # first generate some more complex parts
  # code
  code_text <-
    if (!is.na(definition$code_location)) {
      glue::glue('The original code of this method is available [here]({definition$code_location}).')
    } else {
      ""
    }

  # citation
  citation_text <-
    if (!is.na(definition$DOI)) {
      paste0("The method is described in: [", rcrossref::cr_cn(dois = definition$DOI[[1]], format = "text", style="elsevier-harvard"), "](https://doi.org/", definition$DOI, ")")
    } else {
      ""
    }

  # url within name
  url_name <- if(!is.na(definition$DOI)) {
    glue::glue("[{definition$name}](https://doi.org/{definition$DOI[[1]]})")
  } else if (!is.na(definition$code_location)) {
    glue::glue("[{definition$name}](https://doi.org/{definition$code_location})")
  } else {
    definition$name
  }

  # parameters
  params_text <- map2(names(definition$parameters), definition$parameters, function(parameter_id, parameter) {
    if (parameter_id != "forbidden") {
      values_text <- if (!is.null(parameter$values)) {
        paste0("possible values: ", glue::collapse(parameter$values, ", "))
      } else if (!is.null(parameter$lower) && !is.null(parameter$upper)){
        paste0("possible values between ", parameter$lower, " and ", parameter$upper, "")
      }

      if (is.null(parameter$description)) {
        parameter$description <- ""
      }

      # escape {} for glue
      parameter$description <- parameter$description %>%
        str_replace_all("\\{", "{{") %>%
        str_replace_all("\\}", "}}")

      # remove newlines
      parameter$description <- parameter$description %>%
        str_replace_all("\n", "")

      c(
        glue::glue("@param {parameter_id} {Hmisc::capitalize(parameter$description)} \\cr "),
        glue::glue("    {parameter$type}; default: {deparse(parameter$default)}; {values_text}")
      )
    }
  }) %>% unlist()

  # now combine everything
  documentation <- c(
    # title
    "Inferring a trajectory inference using {url_name}",
    "",
    # description, including citations & code link
    "Will generate a trajectory using {url_name}. This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/{method_id}).",
    "",
    code_text,
    "",
    citation_text,
    "",
    params_text,
    "",
    "@return The trajectory model",
    "@export"
  ) %>%
    paste0("#' ", .) %>%
    glue::collapse("\n")

  glue::glue_data(definition, documentation)
}

# generate parameters
generate_parameters <- function(definition) {
  map2_chr(
    names(definition$parameters),
    definition$parameters,
    function(parameter_id, parameter) {
      if (parameter_id != "forbidden") {
        glue::glue("{parameter_id} = {deparse(parameter$default, width.cutoff = 500)}")
      } else {
        ""
      }
    }
  ) %>% paste0("    ", .) %>% glue::collapse(",\n")
}


generate_func <- function(method_id, definition) {
  func <- "
ti_{method_id} <- function(
{generate_parameters(definition)}
) {{
  args <- as.list(environment())
  method <- create_docker_ti_method('dynverse/{method_id}')
  do.call(method, args)
}}
" %>% glue::glue()
}


#   ____________________________________________________________________________
#   Generate documentation and function for each method                     ####
get_method_definition <- function(method_id) {
  definition <- extract_definition_from_docker_image(paste0("dynverse/", method_id), docker_client = stevedore::docker_client())

  if (!method_id %in% methods_info$method_id) {stop(method_id, " not found in google sheet!")}

  method_info <- methods_info %>% slice(match(!!method_id, method_id)) %>% dynutils::extract_row_to_list(1)

  definition <- c(definition, method_info)

  definition
}


file_location <- "R/ti_container.R"
write_file("
################################### DO NOT EDIT #####################################
#### This file is automatically generated from 2_build_and_generate_containers.R ####
#####################################################################################
\n\n\n
", file_location)

for (method_id in method_ids) {
  print(method_id)

  definition <- get_method_definition(method_id)
  documentation <- generate_documentation_from_definition(method_id, definition)
  func <- generate_func(method_id, definition)

  # now save to file
  file <- paste(
    documentation,
    func,
    sep="\n"
  )

  write_file(file, file_location, append = TRUE)
  write_file("\n\n\n\n", file_location, append = TRUE)
}

#   ____________________________________________________________________________
#   Save data on containerised methods                                      ####
methods_containerised <- tibble(
  method_id = method_ids,
  docker_container = paste0("dynverse/", method_ids),
  containerised = TRUE
) %>%
  mutate(
    docker_wrapper_location = paste0("containers/", method_id),
    dockerhub_location = paste0("https://hub.docker.com/r/", docker_container)
  )
usethis::use_data(methods_containerised, overwrite = TRUE)
