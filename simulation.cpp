// #include <Rcpp.h>
#include <vector>
#include <random>
// using namespace Rcpp;

std::random_device rd;
std::mt19937 generator(rd());

// [[Rcpp::export]]
std::vector<double> sim_variation_cpp (int N,
                                       double mu,
                                       double s,
                                       double h,
                                       int gen,
                                       double q0) {
  std::vector<double> q (gen, 0);
  q[0] = q0;


  for (int i = 1; i < gen; i++) {

    std::binomial_distribution<int> binom(2, q[i - 1]);
    std::vector<int> geno (N, 0);
    for (int j = 0; j < N; j++) {
      geno[j] = binom(generator);
    }

    std::vector<double> fitness (N, 1);
    for (int j = 0; j < N; j++) {
      switch (geno[j]) {
        case 1: fitness[j] = 1 - h * s; break;
        case 2: fitness[j] = 1 - s; break;
      }
    }

    std::uniform_real_distribution<> unif(0, 1);
    std::vector<double> survival_roll (N, 0);
    for (int j = 0; j < N; j++) {
      survival_roll[j] = unif(generator);
    }

    std::vector<bool> survival (N, false);
    for (int j = 0; j < N; j++) {
      if (survival_roll[j] < fitness[j]) {
        survival[j] = true;
      }
    }

    std::vector<int> surviving;
    for (int j = 0; j < N; j++) {
      if (survival[j]) {
        surviving.push_back(geno[j]);
      }
    }

    std::vector<double> mutation_roll (N, 0);
    for (int j = 0; j < surviving.size(); j++) {
      mutation_roll[j] = unif(generator);
    }

    std::vector<int> mutated (surviving);
    for (int j = 0; j < surviving.size(); j++) {
      if (mutation_roll[j] < mu && surviving[j] == 0) {
        mutated[j] = 1;
      }
      if (mutation_roll[j] < mu && surviving[j] == 1) {
        mutated[j] = 2;
      }

    }

    double n_P = 0;
    double n_H = 0;
    double n_Q = 0;
    for (int j = 0; j < mutated.size(); j++) {
      switch(mutated[j]) {
        case 0: n_P++; break;
        case 1: n_H++; break;
        case 2: n_Q++; break;
      }
    }
    
    q[i] = (n_H + n_Q * 2) / (2 * (n_P + n_H + n_Q));
  }
  return q;
}
