# here we do some metaprogramming to generate the ti_{method} functions inside ti_containers.R

library(tidyverse)
library(dynwrap)
library(furrr)
plan(multiprocess)

load("data/methods_info.rda")

#   ____________________________________________________________________________
#   Run & push dockers                                                      ####
# get all methods
method_ids <- list.dirs("./containers/", full.names = FALSE)[-1]

if (!all(method_ids %in% methods_info$method_id)) {
  stop("Some methods not in methods_info! \n", setdiff(method_ids, methods_info$method_id))
}

# rebuild & push all containers
rebuild <- TRUE
if (rebuild) {
  future_map(method_ids, function(method_id) {
    system(str_glue("docker build containers/{method_id} -t dynverse/{method_id}"))
    system(str_glue("docker push dynverse/{method_id}"))
  })
}

#   ____________________________________________________________________________
#   Functions to generate documentation and parameters given a definition   ####
generate_documentation_from_definition <- function(method_id, definition, r_wrapped = FALSE) {
  # first generate some more complex parts
  # r wrapped
  r_wrapped_text <- if(r_wrapped) {
    glue::glue('This methods was first wrapped inside R, see [ti_{method_id}]')
  } else {
    ""
  }


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
  # do not add parameters documentation if r_wrapped, as the documentation is than already contained at the original function location
  params_text <- if (!r_wrapped) {
    map2(names(definition$parameters), definition$parameters, function(parameter_id, parameter) {
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
  } else {
    params_text <- "@param docker Whether to use the docker container or the R wrapper"
  }

  # now combine everything
  documentation <- c(
    # title
    "Inferring a trajectory inference using {url_name}",
    "",
    # description, including citations & code link
    "Will generate a trajectory using {url_name}. This method was wrapped inside a [container](https://github.com/dynverse/dynmethods/tree/master/containers/{method_id}).",
    "",
    r_wrapped_text,
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
      # do not include forbidden
      if (parameter_id != "forbidden") {
        glue::glue("{parameter_id} = {deparse(parameter$default, width.cutoff = 500)}")
      } else {
        NA
      }
    }
  ) %>% discard(is.na) %>% paste0("    ", .) %>% glue::collapse(",\n")
}


generate_func <- function(method_id, definition, r_wrapped) {
  if (r_wrapped) {
    func <- "
ti_{method_id} <- create_ti_method_chooser(ti_{method_id}, 'dynverse/{method_id}')
" %>% glue::glue()
  } else {
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
  func
}


#   ____________________________________________________________________________
#   Generate documentation and function for each method                     ####
get_method_definition <- function(method_id) {
  definition <- extract_definition_from_docker_image(paste0("dynverse/", method_id))

  if (!method_id %in% methods_info$method_id) {stop(method_id, " not found in google sheet!")}

  method_info <- methods_info %>% slice(match(!!method_id, method_id)) %>% dynutils::extract_row_to_list(1)

  definition <- c(definition, method_info)

  definition
}

# beginning of the file
file_location <- "R/ti_container.R"
write_file("
################################### DO NOT EDIT #####################################
#### This file is automatically generated from 2_build_and_generate_containers.R ####
#####################################################################################
\n
", file_location)

# make sure all other ti_* files are loaded first, so that the ti_* functions can be used for the containers
write_file(paste0(
  "#' @include ",
  list.files("R") %>% str_subset("ti_.*") %>% str_subset("^(?!ti_container.R)") %>% glue::collapse(" "),
  "\n",

  "#' @include wrapper_method_chooser.R\n"
), file_location, append = TRUE)

# get methods which were wrapped inside R
lines <- list.files("R", full.names = TRUE) %>% map(read_lines) %>% unlist()
method_ids_r <- str_subset(lines, "^ti_[A-Za-z0-9_]* <-.*") %>%
  str_replace_all("ti_([A-Za-z0-9_]*) <-.*", "\\1")

# now generate the ti_* functions
for (method_id in method_ids) {
  print(method_id)

  r_wrapped <- method_id %in% method_ids_r

  definition <- get_method_definition(method_id)
  documentation <- generate_documentation_from_definition(
    method_id,
    definition,
    r_wrapped = r_wrapped
    )
  func <- generate_func(method_id, definition, r_wrapped)

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
