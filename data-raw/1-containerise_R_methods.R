# here we generate containers for all methods which are directly implemented inside the dynmethods package
# this will use the wrappers contained in for example R/ti_scorpius.R and will save new files inside containers/scorpius

library(tidyverse)
library(dynwrap)
library(desc)

write_file("", "R/ti_container.R")

devtools::load_all()

method_ids <- dynwrap::get_ti_methods(packages = "dynmethods")$method_id
method_ids <- "mpath"

walk(method_ids, function(method_id) {
  cat("Running ", method_id, "\n", sep = "")
  method <- get(paste0("ti_", method_id))()

  if (is.null(method$run_fun_name) || !grepl("dynmethods::", method$run_fun_name)) {
    return(NULL)
  }

  # generate the definition file
  get_definition <- function(method) {
    if (is.null(method$parameters)) {
      stop(method$short_name, " does not have a list of parameters!")
    }

    definition <- method %>% keep(~is.numeric(.) || is.character(.))

    definition <- c(
      definition,
      list(
        parameters = method$parameters,
        input = list(
          format = "rds"
        ),
        output = list(
          format = "dynwrap"
        )
      )
    )

    if (!is.null(method$authors)) {
      definition$authors <- method$authors
    }

    optional_inputs <- method$inputs %>% filter(!required, type != "parameter") %>% pull(input_id)
    if (length(optional_inputs) > 0) {definition$input$optional <- optional_inputs}

    required_inputs <- method$inputs %>% filter(required, type != "parameter") %>% pull(input_id)
    if (length(required_inputs) > 0) {definition$input$required <- required_inputs}

    definition[!names(definition) %in% c("method_id", "run_fun_name")]
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
    if (is.null(method$run_fun_name) || !grepl("dynmethods::", method$run_fun_name)) {
      return(NA)
    }

    run_fun_name <- gsub("dynmethods::", "", method$run_fun_name)

    df <- list.files("R", pattern = "ti_", full.names = T) %>%
      map_df(function(file) {
        file_text <- readr::read_lines(file)
        line_numbers_start <- which(str_detect(file_text, paste0("^", run_fun_name, " <- function")))

        if (length(line_numbers_start) == 0) {
          NULL
        } else {
          line_numbers_middle <- which(str_detect(file_text, "^\\) +\\{ *$")) %>% keep(. >= line_numbers_start) %>% head(1)
          line_numbers_end <- which(str_detect(file_text, "^\\} *$")) %>% keep(. >= line_numbers_middle) %>% head(1)

          formals <- formals(method$run_fun)
          default_params <- paste0("  ", mapply(
            names(formals),
            formals,
            FUN = function(n, v) {
              if (is.symbol(v)) n else paste0(n, " = ", deparse(v))
            }
          ), collapse = ",\n")


          # deparse_text <- deparse(method$run_fun)
          # line_numbers <- which(deparse_text == "{")
          # paste0(deparse_text[seq_len(line_numbers-1)], collapse = " ")

          run_fun_code <- paste0(c(
            "run_fun <- function(",
            default_params,
            file_text[seq(line_numbers_middle, line_numbers_end)]
          ), collapse = "\n")

          data_frame(run_fun_code)
        }
      })

    if (nrow(df) == 0) {
      return(NA)
    }

    run_fun_code <- df$run_fun_code[[1]]

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

{run_fun_code}

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
})
