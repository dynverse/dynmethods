#' Install remotes using git tags and R versions
#'
#' @param remotes A list of remotes in the format of `parse_github_repo_spec`
#' @param is_interactive Whether running in an interactive session
#'
#' @importFrom purrr map_chr map2_dbl map_chr map remotes
#' @importFrom remotes parse_github_repo_spec
install_github_tagged_version <- function(remotes, is_interactive = interactive()) {
  parsed <- purrr::map(remotes, parse_github_repo_spec) %>% purrr::set_names(remotes)

  current_versions <- parsed %>% purrr::map_chr("ref")
  requested_versions <- purrr::map_chr(parsed, function(x) devtools::package_info(x$package, dependencies = NULL)$ondiskversion)

  version_comparisons <- purrr::map2_dbl(current_versions, requested_versions, compareVersion)

  to_install <- remotes[version_comparisons != 0]

  # only install when we need to install
  if (length(to_install) > 0) {

    # only prompt user if interactive
    if (is_interactive) {
      message(paste0(
        "Following packages have to be installed: ",
        glue::glue_collapse(crayon::bold(parsed[[to_install]]$package), ", ", last = " and "),
        "\n",
        "Do you want to install these packages? \n",
        "1: Yes [default]\n",
        "2: No"
      ))

      answer <- readline()

      if (answer %in% c("2", "n", "no", "No")) {
        stop("Installation was interrupted.")
      }
    }

    walk(to_install, remotes::install_github)
  }

  TRUE
}

parse_github_repo_spec <- function(x) {
  parsed <- remotes::parse_github_repo_spec(x)
  if(parsed$package == "")
    parsed$package <- parsed$repo

  parsed
}
