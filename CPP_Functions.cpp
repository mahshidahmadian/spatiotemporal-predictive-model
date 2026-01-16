////////////////////////////////////////////////////////////////////////////////
//                     C++ FUNCTIONS
////////////////////////////////////////////////////////////////////////////////
// C++ functions for paper:
// *A New Perspective to Fish Trajectory Imputation:*
// *A Spatiotemporal Probability Model for Simulating Acoustically Tagged Fish Movement*
// 
// Link to the Paper: *https://arxiv.org/abs/2408.13220*
//
////////////////////////////////////////////////////////////////////////////////

#define RCPP_NO_MODULES
#include <Rcpp.h>
#include <cmath>
#include <numeric>
using namespace Rcpp;

//------------------------------------------------------------------------------
// HELPER FUNCTIONS
//------------------------------------------------------------------------------

/**
 * Calculate standard deviation for (psi)
 */
double sd_psi_calculator2(double gamma1, int nj, double t1, double beta1, double u1) {
  return (beta1 <= nj) ? gamma1 * exp((nj - (t1 - 1))) : gamma1 * u1;
}


// [[Rcpp::export]]
List generate_one_segment(NumericVector start, NumericVector end,
                                       NumericMatrix reciever_cors, double gamma1, 
                                       double sigma_R2, 
                                       int nj, double r, int beta1, 
                                       double phi, SEXP ocean) {
  double Log_likelihood1 = 0.0;
  int T_tot = nj + 2;
  NumericMatrix rp1(T_tot, 2);
  
  // trajectory parameters:
  NumericVector mu_d(nj), theta1(nj), sigma_d(nj), sigma_psi(nj);
  
  // start and end of segment
  // we call'impute_receiver_observed'
  Function impute_R = Environment::global_env()["impute_receiver_observed"];
  NumericMatrix Sigma_R(2, 2);
  Sigma_R(0, 0) = sigma_R2;
  Sigma_R(1, 1) = sigma_R2;
  
  NumericVector ss_start = impute_R(start[0], start[1], Sigma_R, r);
  rp1(0, _) = ss_start;
  
  NumericVector ss_end = impute_R(end[0], end[1], Sigma_R, r);
  rp1(T_tot - 1, _) = ss_end;
  
  // likelihood (above)
  Function like_pos = Environment::global_env()["log_like_pos"];
  NumericMatrix endpoints(2, 2);
  endpoints(0, _) = ss_start; endpoints(1, _) = ss_end;
  
  NumericMatrix rec_pts(2, 2);
  rec_pts(0, _) = start; rec_pts(1, _) = end;
  
  Log_likelihood1 += as<double>(like_pos(endpoints, rec_pts, std::sqrt(sigma_R2)));
  
  // loop for segment generation
  NumericVector current_pos = ss_start;
  for (int j = 1; j <= nj; ++j) {
    double u1 = R::runif(0, beta1);
    
    // distance to target logic
    double dist_to_end = std::sqrt(std::pow(ss_end[0] - current_pos[0], 2) + 
                                   std::pow(ss_end[1] - current_pos[1], 2));

    mu_d[j-1] = dist_to_end / (nj - j + 1);
    sigma_d[j - 1] = std::pow(dist_to_end, phi);
    
    // Angle to endpoint
    NumericMatrix dir_vec(1, 2);
    dir_vec(0, 0) = ss_end[0] - current_pos[0];  // dx
    dir_vec(0, 1) = ss_end[1] - current_pos[1];  // dy
    Function angle_finder_R = Environment::global_env()["angle_finder"];
    theta1[j - 1] = as<double>(angle_finder_R(dir_vec));
    
    //  Calculate psi
    sigma_psi[j-1] = sd_psi_calculator2(gamma1, nj, j, beta1, u1);
    
    Function generate_valid_move = Environment::global_env()["generate_valid_move"];
    List step_res = generate_valid_move(current_pos, reciever_cors, theta1[j-1], 
                             sigma_psi[j-1], sigma_d[j - 1], r, mu_d[j-1], ocean, 1000);
    
    current_pos = as<NumericVector>(step_res["Imputed_cor"]);
    double psi_val = as<double>(step_res["psi"]);
    double d_val = as<double>(step_res["Djt"]);
    rp1(j, _) = current_pos;
    
    // Likelihood
    Function like_dist = Environment::global_env()["log_like_dist"];
    Function like_z = Environment::global_env()["log_like_z"];
    
    //NumericVector z_vec = {std::cos(theta1[j-1] + psi), std::sin(theta1[j-1] + psi)};
    Log_likelihood1 += as<double>(like_dist(d_val, mu_d[j-1], sigma_d[j - 1]));
    NumericVector row_prev = rp1(j-1, _);
    NumericVector row_curr = rp1(j, _);
    Log_likelihood1 += as<double>(like_z(row_curr, row_prev, d_val, theta1[j-1], sigma_psi[j-1]));
  }
  
  return List::create(Named("path_points") = rp1, Named("likelihood") = Log_likelihood1);
}


// [[Rcpp::export]]
List impute_multiple_trajectories(int NS, double percent1, NumericMatrix observed_cors,
                                  NumericMatrix reciever_cors, double gamma_mean, double gamma_sd,
                                  double sigma_R2, IntegerVector gap_vec,
                                  double r, int beta_param,
                                  double phi_mean, double phi_sd, SEXP ocean) {
  GetRNGstate();
  
  int path_segments = observed_cors.nrow() - 1;
  int Ntot = observed_cors.nrow() + sum(gap_vec);
  int NKeep = std::max(1, (int)std::round(NS * percent1)); // Ensure at least 1 path is kept
  
  // Storage for all generated paths and their Log-Likelihoods
  // Dimension: [Ntot rows, 2 columns (lon/lat), NS simulations]
  NumericVector all_paths(Ntot * 2 * NS);
  all_paths.attr("dim") = Dimension(Ntot, 2, NS);
  NumericVector log_likes(NS);
  
  // Access the segment generator (using the refined name from earlier)
  Environment env = Environment::global_env();
  Function generate_segment = env["generate_one_segment"];
  
  for (int s = 0; s < NS; ++s) {
    // 1. Stochastic Parameter Sampling
    double phi   = R::rlnorm(phi_mean, phi_sd);
    int beta1    = R::rpois(beta_param); 
    double gamma = R::rlnorm(gamma_mean, gamma_sd);
    
    NumericMatrix full_trajectory(Ntot, 2);
    int current_row = 0;
    double total_log_like = 0.0;
    
    // 2. Stitch Segments Together
    for (int i = 0; i < path_segments; ++i) {
      NumericVector start = observed_cors(i, _);
      NumericVector end   = observed_cors(i + 1, _);
      int nj = gap_vec[i];
      
      List segment_res = generate_one_segment(start, end, reciever_cors, gamma,
                                              sigma_R2, nj, r, beta1, phi, ocean);
      
      NumericMatrix segment_pts = as<NumericMatrix>(segment_res["path_points"]);
      total_log_like += as<double>(segment_res["likelihood"]);
      
      // Copy points to the full trajectory, avoiding the duplicate endpoint
      int start_idx = (i == 0) ? 0 : 1;
      int num_pts = segment_pts.nrow();
      
      for (int k = start_idx; k < num_pts; ++k) {
        full_trajectory(current_row, 0) = segment_pts(k, 0);
        full_trajectory(current_row, 1) = segment_pts(k, 1);
        current_row++;
      }
    }
    
    // 3. Store Results (Flattened indexing for 3D array)
    for (int row = 0; row < Ntot; ++row) {
      all_paths[row + 0 * Ntot + s * Ntot * 2] = full_trajectory(row, 0);
      all_paths[row + 1 * Ntot + s * Ntot * 2] = full_trajectory(row, 1);
    }
    log_likes[s] = total_log_like;
  }
  
  // 4. Selection Logic (Importance Sampling)
  // Sort indices based on log-likelihood descending
  std::vector<int> indices(NS);
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
            [&](int a, int b) { return log_likes[a] > log_likes[b]; });
  
  // Final Output Storage
  NumericVector kept_paths(Ntot * 2 * NKeep);
  kept_paths.attr("dim") = Dimension(Ntot, 2, NKeep);
  
  for (int k = 0; k < NKeep; ++k) {
    int top_idx = indices[k];
    for (int row = 0; row < Ntot; ++row) {
      kept_paths[row + 0 * Ntot + k * Ntot * 2] = all_paths[row + 0 * Ntot + top_idx * Ntot * 2];
      kept_paths[row + 1 * Ntot + k * Ntot * 2] = all_paths[row + 1 * Ntot + top_idx * Ntot * 2];
    }
  }
  
  PutRNGstate();
  return List::create(
    Named("imputed_trajectories") = kept_paths,
    Named("log_likelihood")  = log_likes[indices[0]]
  );
}


