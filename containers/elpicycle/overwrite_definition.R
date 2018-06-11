library(yaml)

definition <- read_yaml("/code/definition.yml")
definition$name = "elpilinear"
definition$parameters$topology$default <- "cycle"
definition$parameters$topology$values <- "cycle"
definition$input$optional <- NULL

write_yaml(definition, "/code/definition.yml")
