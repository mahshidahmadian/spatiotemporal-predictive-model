# Install dependencies
install.packages(c("Rcpp", "mvtnorm", "VGAM", "sf", "dplyr", "ggplot2", 
                   "leaflet", "rnaturalearth"))

# Load code
library(Rcpp)
Rcpp::sourceCpp("c++/CPP_Functions.cpp")
source("R/R_Functions.R")

# Load data
fish_data <- read.csv("datasets/fish_data.csv")
receiver_data <- read.csv("datasets/receiver_data.csv")

# Run imputation (for fish 1)
fish1 <- fish_data[fish_data$tagname == "A69-9001-18453", ]
observed_locs <- as.matrix(fish1[fish1$miss_id == 0, c("longitude", "latitude")])
receiver_coords <- as.matrix(receiver_data[, c("longitude", "latitude")])

results <- impute_multiple_trajectories(
  NS = 1000,
  percent1 = 0.9,
  observed_cors = observed_locs,
  reciever_cors = receiver_coords,
  gamma_mean = 3,
  gamma_sd = 0.1,
  sigma_R2 = 0.00148,
  gap_vec = gap_count(fish1$miss_id),
  r = 0.004491556,
  beta_param = 3,
  phi_mean = 0.5,
  phi_sd = 0.01,
  ocean = coastline  # Load with rnaturalearth::ne_download()
)

# Visualize
cobia_map(results$imputed_trajectories, receiver_coords)


