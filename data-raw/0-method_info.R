library(googlesheets)
library(tidyverse)

#   ____________________________________________________________________________
#   Load metadata of every method from google sheets                        ####
methods_info <- gs_key("1Mug0yz8BebzWt8cmEW306ie645SBh_tDHwjVw4OFhlE") %>%
  gs_read(ws = "Methods", skip = 1)

implementations_info <- gs_key("1Mug0yz8BebzWt8cmEW306ie645SBh_tDHwjVw4OFhlE") %>%
  gs_read(ws = "Implementations", skip = 1)

methods_info <- left_join(methods_info, implementations_info, "implementation_id")

usethis::use_data(methods_info, overwrite=TRUE)
