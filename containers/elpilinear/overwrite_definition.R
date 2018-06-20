library(yaml)

definition <- read_yaml("/code/definition.yml")
definition$name = "ElPiGraph linear"
definition$short_name = "elpilinear"
definition$parameters$topology$default <- "linear"
definition$parameters$topology$values <- "linear"
definition$input$optional <- NULL

write_yaml(definition, "/code/definition.yml")
