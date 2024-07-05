#include <Rcpp.h>
#include <unordered_set>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Function to get unique values from a NumericVector
// [[Rcpp::export]]
NumericVector unique_cpp(NumericVector x) {
  std::unordered_set<double> seen;
  std::vector<double> unique_vals;
  
  for (double val : x) {
    if (seen.insert(val).second) {
      unique_vals.push_back(val);
    }
  }
  
  return wrap(unique_vals);
}

// Function to find the best split
double find_best_split(NumericVector x, double split_point) {
  double right = R_PosInf;
  double left = R_NegInf;
  
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] > split_point && x[i] < right) {
      right = x[i];
    }
    if (x[i] <= split_point && x[i] > left) {
      left = x[i];
    }
  }
  
  return (left + right) / 2;
}

// [[Rcpp::export]]
NumericVector search_split_cpp(NumericVector x, NumericMatrix y, NumericVector splits) {
  int N = y.nrow();
  int G = y.ncol();
  
  // Sort x and get sorted indices
  IntegerVector ord = Rcpp::seq(0, N - 1);
  std::sort(ord.begin(), ord.end(), [&](int i, int j) { return x[i] < x[j]; });
  
  NumericVector x_sorted(N);
  NumericMatrix y_sorted(N, G);
  for (int i = 0; i < N; ++i) {
    x_sorted[i] = x[ord[i]];
    for (int j = 0; j < G; ++j) {
      y_sorted(i, j) = y(ord[i], j);
    }
  }
  
  // Cumulative sums
  NumericMatrix cum_sum_y(N, G);
  NumericMatrix cum_sum_y2(N, G);
  for (int j = 0; j < G; ++j) {
    cum_sum_y(0, j) = y_sorted(0, j);
    cum_sum_y2(0, j) = y_sorted(0, j) * y_sorted(0, j);
    for (int i = 1; i < N; ++i) {
      cum_sum_y(i, j) = cum_sum_y(i - 1, j) + y_sorted(i, j);
      cum_sum_y2(i, j) = cum_sum_y2(i - 1, j) + y_sorted(i, j) * y_sorted(i, j);
    }
  }
  
  // Calculate the objective values for all splits
  int num_splits = splits.size();
  NumericVector objectives(num_splits, R_PosInf);
  for (int k = 0; k < num_splits; ++k) {
    double split = splits[k];
    int N_L = 0;
    while (N_L < N && x_sorted[N_L] <= split) {
      N_L++;
    }
    if (N_L == 0 || N_L == N) {
      continue; // Invalid split
    }
    
    double var_L = 0.0, var_R = 0.0;
    for (int j = 0; j < G; ++j) {
      double S_L = cum_sum_y(N_L - 1, j);
      double S_R = cum_sum_y(N - 1, j) - S_L;
      double SS_L = cum_sum_y2(N_L - 1, j);
      double SS_R = cum_sum_y2(N - 1, j) - SS_L;
      
      var_L += SS_L - (S_L * S_L) / N_L;
      var_R += SS_R - (S_R * S_R) / (N - N_L);
    }
    
    objectives[k] = var_L + var_R;
  }
  
  return objectives;
}

// [[Rcpp::export]]
DataFrame best_split_cpp(DataFrame X, NumericMatrix y, Nullable<List> splits_list = R_NilValue) {
  int n_features = X.size();
  CharacterVector feature_names = X.names();
  
  NumericVector split_points(n_features);
  NumericVector objective_values(n_features);
  LogicalVector best_split(n_features);
  CharacterVector features(n_features);
  
  for (int i = 0; i < n_features; ++i) {
    NumericVector x = X[i];
    
    NumericVector split_vector;
    if (splits_list.isNotNull()) {
      List splits = as<List>(splits_list);
      split_vector = splits[i];
    } else {
      split_vector = unique_cpp(x);
    }
    
    NumericVector objectives = search_split_cpp(x, y, split_vector);
    int best = which_min(objectives);
    double best_split_point = split_vector[best];
    
    split_points[i] = find_best_split(x, best_split_point);
    
    objective_values[i] = objectives[best];
    features[i] = feature_names[i];
  }
  
  best_split = (objective_values == min(objective_values));
  
  return DataFrame::create(
    Named("feature") = features,
    Named("split.points") = split_points,
    Named("objective.value") = objective_values,
    Named("best.split") = best_split
  );
}