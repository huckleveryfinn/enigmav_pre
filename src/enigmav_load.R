#
# enigmav_load.R
#
# created on Fri Nov  4 17:38:56 2022
# Philipp Homan, <philipp dot homan at bli dot uzh dot ch>
#-----------------------------------------------------------------------
#

# load functions
source("enigmav_func.R")


# load main data file
evp <- read_csv("../data/enigmav_data.csv")

# load covariates file
cov <- read_csv("../data/enigmav_covariates.csv")
