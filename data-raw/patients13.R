## code to prepare `patients13` dataset goes here
patients13 <- 
  tibble::tribble(
    ~id, ~surgery, ~time_to_surgery, ~death, ~followup,
    "A", 1, 61, 1, 300,
    "B", 1, 150, 0, 365,
    "C", 1, 20, 1, 150,
    "D", 1, 140, 0, 160,
    "E", 1, 10, 0, 300,
    "F", 1, 220, 1, 320,
    "G", 1, 340, 0, 365,
    "H", 1, 200, 0, 300,
    "I", 0, NA, 0, 365,
    "J", 0, NA, 1, 260,
    "K", 0, NA, 1, 40,
    "L", 0, NA, 0, 140,
    "M", 0, NA, 0, 320
  )

usethis::use_data(patients13, overwrite = TRUE)
