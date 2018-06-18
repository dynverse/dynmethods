# here we generate containers for all methods which are directly implemented inside the dynmethods package
# this will use the wrappers contained in for example R/ti_scorpius.R and will save new files inside containers/scorpius

library(tidyverse)
library(dynwrap)
library(googlesheets)
library(desc)

load("data/methods_info.rda")

write_file("", "R/ti_container.R")

method_id <- "monocle_ddrtree"
devtools::load_all()
method <- get(paste0("ti_", method_id), asNamespace("dynmethods"))()

# generate the definition file
get_definition <- function(method) {
  if (is.null(method$parameters)) {
    stop(method$short_name, " does not have a list of parameters!")
  }

  definition <- list(
    name = method$name,
    short_name = method$short_name,
    parameters = method$parameters,
    input = list(
      format = "rds"
    ),
    output = list(
      format = "dynwrap"
    )
  )

  optional_inputs <- method$inputs %>% filter(!required, type != "parameter") %>% pull(input_id)
  if (length(optional_inputs) > 0) {definition$input$optional <- optional_inputs}

  required_inputs <- method$inputs %>% filter(required, type != "parameter") %>% pull(input_id)
  if (length(required_inputs) > 0) {definition$input$required <- required_inputs}

  definition
}

# generate the docker file
get_dockerfile <- function(method) {
  remotes <- desc::desc_get_remotes()
  remotes <- set_names(remotes, remotes %>% str_replace(".*/([^@]*).*", "\\1"))
  dependencies <- c(method$package_required, method$package_loaded)

  if (length(dependencies) > 0) {
    install_dependencies <- map(dependencies, function(dependency) {
      if(dependency %in% names(remotes)) {
        glue::glue("devtools::install_github('{remotes[dependency]}')")
      } else {
        glue::glue("install.packages('{dependency}')")
      }
    }) %>%
      paste0("RUN R -e \"", ., "\"") %>%
      glue::collapse("\n")
  } else {
    install_dependencies <- ""
  }

  install_apt <- if(!is.null(method$apt_dependencies)) {
    glue::glue("RUN apt-get install -y {glue::collapse(method$apt_dependencies, ' ')}")
  } else {
    ""
  }

glue::glue("
FROM rocker/tidyverse

RUN echo 'utils::setRepositories(ind=1:4)' > ~/.Rprofile

RUN R -e 'devtools::install_github(\"dynverse/dynwrap\")'

{install_apt}

{install_dependencies}

ADD . /code

ENTRYPOINT Rscript code/run.R
")
}

# generate the run.R file
get_runr <- function(method) {
  glue::glue("
library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

{glue::collapse(paste0('library(', method$package_required, ')'), '\\n')}

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- {method$run_fun %>% deparse() %>% glue::collapse('\\n')}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')
  ")
}

folder <- glue::glue("containers/{method$short_name}")
dir.create(folder, showWarnings = FALSE)
get_definition(method) %>% yaml::write_yaml(file.path(folder, "definition.yml"))
get_dockerfile(method) %>% write_file(file.path(folder, "Dockerfile"))
get_runr(method) %>% write_file(file.path(folder, "run.R"))

# system(glue::glue("docker build containers/{method_id} -t dynverse/{method_id}"))
