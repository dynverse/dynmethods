# here we generate containers for all methods which are directly implemented inside the dynmethods package

library(tidyverse)
library(dynwrap)
library(googlesheets)
library(desc)

load("data/methods_info.rda")

get_ti_methods()

method_id <- "stemid2"

# extract par_set, process to parameters list
file <- paste0("man/ti_", method_id, ".Rd")
rd <- Rd2roxygen::parse_file(file)
parameter_descriptions = rd$params %>% {set_names(gsub("[^ ]* (.*)", "\\1", .), gsub("([^ ]*) .*", "\\1", .))}

parameters <- map2(names(method$par_set$pars), method$par_set$pars, function(id, par) {
  list(
    type = par$type,
    default = par$default,
    upper = par$upper,
    lower = par$lower,
    values = as.character(par$values),
    description = parameter_descriptions[[id]]
  ) %>% discard(is.null) %>% discard(~length(.) == 0)
}) %>% set_names(names(method$par_set$pars))
if (!is.null(method$par_set$forbidden)) {
  parameters$forbidden <- deparse(method$par_set$forbidden)
}

# copy parameters
parameters %>%
  deparse(width.cutoff=500L) %>%
  gsub("([A-Za-z_\\.]* = )", "\n\\1", .) %>%
  clipr::write_clip()

##
devtools::load_all()

method <- get(paste0("ti_", method_id), asNamespace("dynmethods"))()

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

get_dockerfile <- function(method) {
  remotes <- desc::desc_get_remotes()
  remotes <- set_names(remotes, remotes %>% str_replace(".*/(.*)", "\\1"))
  dependencies <- c(method$package_required, method$package_loaded)

  install_dependencies <- map(dependencies, function(dependency) {
    if(dependency %in% names(remotes)) {
      glue::glue("setRepositories(ind=1:2);devtools::install_github('{remotes[dependency]}')")
    } else {
      glue::glue("setRepositories(ind=1:2);install.packages('{dependency}')")
    }
  }) %>%
    paste0("RUN R -e \"", ., "\"") %>%
    glue::collapse("\n")

glue::glue("
FROM rocker/tidyverse

RUN R -e 'devtools::install_github(\"dynverse/dynwrap\")'

{install_dependencies}

ADD . /code

ENTRYPOINT Rscript code/run.R
")
}

get_runr <- function(method) {
  glue::glue("
library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)

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

system(glue::glue("docker build containers/{method_id} -t dynverse/{method_id}"))
