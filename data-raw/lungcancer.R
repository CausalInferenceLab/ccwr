## code to prepare `lungcancer` dataset goes here
library(readr)
lungcancer <- read_csv("data-raw/lungcancer.csv")
attr(lungcancer, "spec") <- NULL
attr(lungcancer, "problems") <- NULL
usethis::use_data(lungcancer, overwrite = TRUE)
