# # This can be removed by adding ParamHelpers to DESCRIPTION
# .onLoad <- function(libname, pkgname){
#   packageStartupMessage("Loading discreteNameToValue in your global environment -- this is a dirty fix.")
#   assign(
#     "discreteNameToValue",
#     ParamHelpers::discreteNameToValue,
#     envir = .GlobalEnv
#   )
# }
