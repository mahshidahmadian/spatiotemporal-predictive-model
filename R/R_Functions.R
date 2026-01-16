

################################################################################
# FISH TRAJECTORY IMPUTATION - R FUNCTIONS
################################################################################
#' R functions for paper:
#' *A New Perspective to Fish Trajectory Imputation:*
#' *A Spatiotemporal Probability Model for Simulating Acoustically Tagged Fish Movement*
#
#' Link to the Paper: *https://arxiv.org/abs/2408.13220*
#'
# Required packages: mvtnorm, VGAM, sf, ggplot2, ggforce, reshape2, leaflet, leaflet.extras
################################################################################

#------------------------------------------------------------------------------#
# MOVEMENT MODEL FUNCTIONS
#------------------------------------------------------------------------------#

#' generate One-Step ahead fish location
#' 
#' @about: 
#' A function that generates a new movement coordinate
#' based on step length (*D*) and turning angle (*psi*), constrained by 
#' receiver proximity and land boundaries.
#'
#' @param prev_step Previous location [lon, lat]
#' @param receiver_cors Matrix of receiver coordinates
#' @param theta1 Target bearing (radians)
#' @param sigma_psi Standard deviation of turning angle
#' @param sigma_d Standard deviation of step length
#' @param r Receiver detection radius
#' @param mu_d Mean step length
#' @param ocean Coastline shapefile (sf object)
#' @param line_buffer Buffer distance for line geometries (meters)
#' 
#'
generate_valid_move <- function(prev_step, reciever_cors, theta1, sigma_psi, sigma_d, 
                                r, mu_d, ocean, line_buffer) {
  
  while (TRUE) {
    # 1. movement parameters
    psi1   <- rnorm(1, 0, sigma_psi)
    theta_star <- reduce_radian(theta1 + psi1)
    D1     <- abs(rnorm(1, mu_d, sigma_d))
    
    # 2. propose new candidate coordinate
    # prev_step[1] is Longitude/X, prev_step[2] is Latitude/Y
    new_cor <- c(prev_step[1] + cos(theta_star) * D1, 
                 prev_step[2] + sin(theta_star) * D1)
    
    # 3. spatial constraints Check: 
    # recievers constraint
    dist_sq  <- apply(reciever_cors,1,
                      function(coord){euc_distance(new_cor[1], new_cor[2], coord[1], coord[2])})
    is_near_rec <- any(dist_sq < r)
    
    if (is_near_rec) {
      next # if not don't do land check
    }
    
    # land boundaries
    land_check <- check_land_constraint(new_cor, ocean, distance_to_check = D1, line_buffer)
    
    if (land_check == 1) {
      next
    }
    # Output:
    return(list(
      Imputed_cor    = new_cor, 
      psi            = psi1, 
      Djt            = D1
    ))
  }
}

#------------------------------------------------------------------------------#

#' Check if a Coordinate is on Land
#'
#' @about: 
#' Determines if a point is on *land* by use of a global shapefile using a 
#' bounding box and a spatial intersection check.
#'
#' @param location_coord Coordinate to check [lon, lat]
#' @param coast_shape Coastline shapefile (sf object)
#' @param distance_to_check Optional search radius (degrees)
#' @param line_buffer Buffer distance for line geometries (meters)
#' 
#' @return Integer: 1 if on land, 0 if in ocean
#'
#'
check_land_constraint <- function(location_coord, coast_shape, 
                                  distance_to_check = NULL, line_buffer = 100) {
  
  # check inputs are valid
  if (!inherits(coast_shape, "sf")) stop("coast_shape must be an sf object")
  
  # make it Coordinate Reference System (CRS)
  point <- st_sfc(st_point(location_coord), crs = st_crs(coast_shape))
  
  # bounding box calculation
  # we consider a ~30km buffer in degrees to subset the global shapefile
  # this speed up the intersection check
  lat_rad <- abs(location_coord[2]) * pi / 180
  deg_per_30km <- 30 / (111.32 * cos(lat_rad))
  
  # consider search radius
  search_radius <- if (is.null(distance_to_check)){deg_per_30km
    }else{
    max(distance_to_check, deg_per_30km)}
  
  lon <- location_coord[1]
  lat <- location_coord[2]
  
  bbox <- st_bbox(c(
    xmin = lon - search_radius, xmax = lon + search_radius,
    ymin = lat - search_radius, ymax = lat + search_radius
  ), crs = st_crs(coast_shape))
  
  # spatial sub
  # we do not check the whole world, we only check land near the point
  coast_subset <- st_filter(coast_shape, st_as_sfc(bbox))
  
  if (nrow(coast_subset) == 0){
    return(0L) # Return integer 0 for ocean (No land found)
  }
  geom_types <- unique(as.character(st_geometry_type(coast_subset)))
  
  if (any(grepl("POLYGON", geom_types))) {
    # check direct intersection
    is_land <- any(st_intersects(point, coast_subset, sparse = FALSE))
    
  } else if (any(grepl("LINE", geom_types))) {
    # Check if point is within line_buffer meters
    point_buf <- st_buffer(point, dist = line_buffer) 
    is_land   <- any(st_intersects(point_buf, coast_subset, sparse = FALSE))
    
  } else {
    is_land <- FALSE
  }
  
  # Return 1 for land, 0 for ocean
  return(as.integer(is_land)) # TRUE -> 1 (land), FALSE -> 0 (ocean)
}

#------------------------------------------------------------------------------#

#' Impute Position Around Receiver
#' 
#' @description
#' this generates a position near a receiver using rejection sampling from
#' a bivariate normal distribution constrained by detection radius.
#' 
#' @param R_lon Receiver longitude
#' @param R_lat Receiver latitude
#' @param Sigma_R Covariance matrix for detection error
#' @param radius Detection radius constraint
#' 
#' @return Vector [lon, lat] within radius of receiver

impute_receiver_observed <- function(R_lon, R_lat, Sigma_R, radius) {
  R_cor <- c(R_lon, R_lat)

  repeat {
    # generate from bivariate normal
    s1 <- rmvnorm(1, mean = R_cor, sigma = Sigma_R)
    
    # distance check
    d1 <- euc_distance(s1[1], s1[2], R_lon, R_lat)
    if (d1 < radius) {
      return(s1)
    }
  }
}
#--------------------------------------------------------------------#
# Log-likelihood for Step Length (Folded Normal)
log_like_dist <- function(dist_val, mu_d, sigma_d) {
  # Using VGAM::dfoldnorm directly on the value
  return(dfoldnorm(dist_val, mean = mu_d, sd = sigma_d, log = TRUE))
}
#--------------------------------------------------------------------#
# Log-likelihood for Position (N)
log_like_pos <- function(N_matrix, R_matrix, sigma_R) {
  # Vectorized calculation for two points
  Sigma_R <- diag(sigma_R^2, 2)
  # Compute log-likelihood for each point separately, then sum
  loglik1 <- dmvnorm(N_matrix[1, ], mean = R_matrix[1, ], sigma = Sigma_R, log = TRUE)
  loglik2 <- dmvnorm(N_matrix[2, ], mean = R_matrix[2, ], sigma = Sigma_R, log = TRUE)
  
  return(loglik1 + loglik2)
}
#--------------------------------------------------------------------#
get_z_params <- function(theta1, sigma_psi) {
  # Mean Vector Approximation
  mu_z <- c(cos(theta1) * (1 - 0.5 * sigma_psi),
            sin(theta1) * (1 - 0.5 * sigma_psi))
  
  # Variance-Covariance Matrix Approximation
  s2 <- sigma_psi
  s4 <- sigma_psi^2
  c1 <- cos(theta1)
  sn <- sin(theta1)
  
  v11 <- (sn^2 * s2) + (0.5 * c1^2 * s4)
  v12 <- (-c1 * sn * s2) + (0.5 * c1 * sn * s4)
  v22 <- (c1^2 * s2) + (0.5 * sn^2 * s4)
  
  return(list(mu = mu_z, sigma = matrix(c(v11, v12, v12, v22), 2, 2)))
}
#--------------------------------------------------------------------#
log_like_z <- function(X_curr, X_prev, D, theta1, sigma_psi) {
  params <- get_z_params(theta1, sigma_psi)
  mean1 <- X_prev + abs(D) * params$mu
  sigma1 <- (D^2 * params$sigma) + diag(1e-9, 2)
  log_like <- dmvnorm(X_curr, mean = mean1, sigma = sigma1, log = TRUE)
  return(log_like)
}

#--------------------------------------------------------------------#
gap_count <- function(miss_id) {
  # Find indices of zeros
  zeros1 <- which(miss_id == 0)
  gap1 <- c() # Initialize an empty vector to store gaps
  # Loop through zero positions and count 1s between them
  for (j in 1:(length(zeros1) - 1)) {
    # Count the number of 1s between consecutive zeros
    gap1 <- c(gap1, sum(miss_id[(zeros1[j] + 1):(zeros1[j + 1] - 1)] == 1))
  }
  return(gap1)
}

#--------------------------------------------------------------------#
#' Reduce angle to [0, 2*pi) range
reduce_radian <- function(theta) {
  reduced_angle <- theta %% (2 * pi)
  return(reduced_angle)
}
#--------------------------------------------------------------------#

#' Calculate Euclidean distance between two points
euc_distance <- function(x1, y1, x2, y2) {
  d <- sqrt(((x2 - x1)^2) + ((y2 - y1)^2))
  return(d)
}
#--------------------------------------------------------------------#

#' Find angle of a direction vector
#' @param dir1 A 1x2 matrix with direction components [dx, dy]
#' @return Angle in radians [0, 2*pi)
angle_finder <- function(dir1){
  if((dir1[1,2]==0) & (dir1[1,1]<0)){
    dir1r <- ((3*pi)/2)
  }else if((dir1[1,2]==0) & (dir1[1,1]>0)){
    dir1r <- (pi/2)
  }else if((dir1[1,1]==0) & (dir1[1,2]>0)){
    dir1r <- 0
  }else if((dir1[1,1]==0) & (dir1[1,2]<0)){
    dir1r <- (pi)
  }else if((dir1[1,2]>0) & (dir1[1,1]<0)){         #2nd Quarter
    dir1r <- (atan( dir1[1,2]/dir1[1,1])) + pi
  }else if((dir1[1,2]<0) & (dir1[1,1]<0)){         #3rd Quarter
    dir1r <- (atan( dir1[1,2]/dir1[1,1])) + pi
  }else if((dir1[1,2]<0) & (dir1[1,1]>0)){         #4th Quarter
    dir1r <- (atan( dir1[1,2]/dir1[1,1])) + (2*pi)
  }else{                                           #1st Quarter
    dir1r <- (atan( dir1[1,2]/dir1[1,1]))
  }
  return(dir1r)
}
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
library(ggplot2)
library(ggforce)
library(sf)
library(reshape2)

plot_Cobia <- function(seen_locs, Receivers_cors, Res1, coastline_sf, zoom1, r, Receiver_numbers) {
  
  # Receiver and gap circles
  receiver_circles <- data.frame(
    x0 = Receivers_cors[, 1],
    y0 = Receivers_cors[, 2],
    r = r
  )
  
  gap_circles <- data.frame(
    x0 = seen_locs[, 1],
    y0 = seen_locs[, 2],
    r = r
  )
  
  # Crop coastline to bounding box with padding
  bbox <- st_bbox(c(
    xmin = min(seen_locs[, 1]) - zoom1,
    xmax = max(seen_locs[, 1]) + zoom1,
    ymin = min(seen_locs[, 2]) - zoom1,
    ymax = max(seen_locs[, 2]) + zoom1
  ), crs = st_crs(coastline_sf))
  coastline_crop <- st_crop(coastline_sf, bbox)
  
  # FIXED: More robust reshaping
  Res1_long <- melt(Res1)
  colnames(Res1_long) <- c("step", "coord", "imp", "value")
  
  # Use dcast instead of reshape for reliability
  Res1_wide <- dcast(Res1_long, step + imp ~ coord, value.var = "value")
  colnames(Res1_wide) <- c("step", "imp", "x", "y")
  
  # Plot
  ggplot() +
    geom_sf(data = coastline_crop, fill = "palegreen4", color = "coral4", linewidth = 0.3) +
    geom_point(data = Res1_wide, aes(x = x, y = y), 
               color = alpha("burlywood3", 0.6), size = 1) +
    geom_point(data = as.data.frame(gap_circles), aes(x = x0, y = y0), 
               color = "lightblue", size = 2) +
    geom_point(data = as.data.frame(receiver_circles), aes(x = x0, y = y0),
               shape = 21, fill = "red", color = "black", size = 0.75) +
    geom_circle(data = receiver_circles, aes(x0 = x0, y0 = y0, r = r), 
                color = "blue", linewidth = 0.5) +
    geom_circle(data = gap_circles, aes(x0 = x0, y0 = y0, r = r), 
                color = "navy", fill = alpha("maroon", 0.4), linewidth = 0.5) +
    geom_text(data = gap_circles, aes(x = x0, y = y0, label = Receiver_numbers),
              vjust = -1, color = 'black', size = 5) +
    coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
             ylim = c(bbox["ymin"], bbox["ymax"]),
             expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "lightblue1", color = NA),
      panel.grid = element_blank(),
      axis.title = element_text(size = 25, face = "bold"),
      axis.text = element_text(size = 15, face = "bold")
    )
}
#--------------------------------------------------------------------#
cobia_map <- function(Imputed_points, Receivers_cors, receiver_labels = NULL) {
  # Extract imputed points into vectors
  Point_vector1 <- as.vector(Imputed_points[, 1, ])
  Point_vector2 <- as.vector(Imputed_points[, 2, ])
  
  # Create data frame for imputed points
  dff1 <- data.frame(x = Point_vector1, y = Point_vector2)
  
  # Create a base map with tiles
  map_imp <- leaflet() %>%
    addProviderTiles(providers$CartoDB.Positron) %>%
    
    # Add markers for Receivers_cors 
    addCircleMarkers(
      lng = Receivers_cors[, 1],
      lat = Receivers_cors[, 2],
      radius = 4,
      color = "black",
      fillColor = "red",
      fillOpacity = 0.8
    ) %>%
    
    # Add a heatmap layer for density representation
    addHeatmap(
      lng = ~x,
      lat = ~y,
      intensity = ~1,
      data = dff1,
      radius = 15,
      blur = 20,
      max = 0.4
    ) %>%
    
    # Add individual imputed point markers
    addCircleMarkers(
      lng = dff1$x,
      lat = dff1$y,
      radius = 0.1,
      color = "darkgreen",
      fillColor = "lightgreen",
      fillOpacity = 1
    )
  
  # FIXED: Only add labels if provided
  if (!is.null(receiver_labels)) {
    map_imp <- map_imp %>%
      addLabelOnlyMarkers(
        lng = Receivers_cors[, 1],
        lat = Receivers_cors[, 2],
        label = as.character(receiver_labels),
        labelOptions = labelOptions(
          noHide = TRUE, 
          textOnly = TRUE, 
          direction = "top", 
          offset = c(0, -10), 
          style = list('font-size' = '15px', 'font-weight' = 'bold', 'color' = 'black')
        )
      )
  }
  
  return(map_imp)
}
#--------------------------------------------------------------------#