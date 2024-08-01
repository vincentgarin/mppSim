##########################
# mppSim pkg development #
##########################

library(usethis)
library(devtools)
library(golem)

# setwd(dir = "D:/Mes Donnees/WD/R/packages/mppSim")

## Dependencies ----
usethis::use_package("simcross")
usethis::use_package("LDcorSV")
usethis::use_package("ggplot2")

## Add helper functions ----
add_fct("LD_computation", pkg = ".", with_test = FALSE)
add_fct("LD_plot", pkg = ".", with_test = FALSE)
add_fct("LD_decay", pkg = ".", with_test = FALSE)


