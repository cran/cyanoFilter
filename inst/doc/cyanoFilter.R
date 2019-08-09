## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----data_and_preprocessing----------------------------------------------
metadata <- system.file("extdata", "2019-03-25_Rstarted.csv", package = "cyanoFilter", 
                        mustWork = TRUE)

