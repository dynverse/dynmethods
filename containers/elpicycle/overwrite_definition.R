library(yaml)

definition <- read_yaml("/code/definition.yml")
definition$name = "ElPiGraph cyclic"
definition$short_name = "elpicycle"
definition$parameters$topology$default <- "cycle"
definition$parameters$topology$values <- "cycle"
definition$input$optional <- NULL

write_yaml(definition, "/code/definition.yml")
