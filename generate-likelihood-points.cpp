#include <algorithm>
#include "src/gaussian-interpolator/GaussianInterpolator.hpp"
#include <chrono>
#include <cmath>
#include <fstream>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <iostream>
#include <sstream>
#include <limits>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <vector>

int main(int argc, char *argv[]) {
  if (argc < 3 || argc > 3) {
    printf("You must provide input\n");
    printf("The input is: \n\nnumber data points;\noutput file name;\n");
    exit(0);
  }

  unsigned N = std::stoi(argv[1]);
  std::string output_file_name = argv[2];
  unsigned i;
  std::vector<likelihood_point> data_points = std::vector<likelihood_point> (N);

  static int counter = 0;
  static gsl_rng* r_ptr_threadprivate;
#pragma omp threadprivate(counter, r_ptr_threadprivate)
  omp_set_dynamic(0);
  omp_set_num_threads(20);

//   BivariateGaussianKernelBasis basis_positive =
//     BivariateGaussianKernelBasis(dx,
// 				 rho_basis,
// 				 sigma_x_basis,
// 				 sigma_y_basis,
// 				 1.0,
// 				 1.0);

  long unsigned seed_init = 10;
   gsl_rng * r_ptr_local;
   const gsl_rng_type * Type;
   gsl_rng_env_setup();
   Type = gsl_rng_default;
   r_ptr_local = gsl_rng_alloc(Type);
   gsl_rng_set(r_ptr_local, seed_init + N);

  int tid=0;
#pragma omp parallel default(none) private(tid, i) shared(r_ptr_local)
  {
    tid = omp_get_thread_num();

    r_ptr_threadprivate = gsl_rng_clone(r_ptr_local);
    gsl_rng_set(r_ptr_threadprivate, tid);

    printf("Thread %d: counter %d\n", tid, counter);
  }


//   std::vector<likelihood_point> points_for_kriging (N);
//   std::ifstream input_file(input_file_name);

//   if (input_file.is_open()) {
//     for (i=0; i<N; ++i) {
//       input_file >> points_for_kriging[i];
//       points_for_kriging[i].likelihood = 0.0;
//     }
//   }

//   std::vector<likelihood_point> points_for_integration (1);

  auto t1 = std::chrono::high_resolution_clock::now();
#pragma omp parallel default(none) private(i) shared(data_points, N, seed_init)
    {
#pragma omp for
      for (i=0; i<N; ++i) {
	long unsigned seed = seed_init + i;
	gsl_rng_set(r_ptr_threadprivate, seed);

	double rho = gsl_rng_uniform(r_ptr_threadprivate)*1.90 - 0.95;
	double log_sigma_x = 1 + gsl_ran_gaussian(r_ptr_threadprivate, 1.0);
	double log_sigma_y = 1 + gsl_ran_gaussian(r_ptr_threadprivate, 1.0);

	double sigma_x = exp(log_sigma_x);
	double sigma_y = exp(log_sigma_y);

	BrownianMotion BM_current = BrownianMotion(seed,
						   1e6,
						   rho,
						   sigma_x,
						   sigma_y,
						   0.0,
						   0.0,
						   1.0);
	data_points[i] = BM_current;
	data_points[i].print_point();

	printf("Thread %d with address -- produces likelihood --\n",
	       omp_get_thread_num());
      }
    }
    auto t2 = std::chrono::high_resolution_clock::now();

//     std::string output_file_name = file_prefix +
//       "-number-points-" + argv[1] +
//       "-rho_basis-" + argv[2] +
//       "-sigma_x_basis-" + argv[3] +
//       "-sigma_y_basis-" + argv[4] +
//       "-dx_likelihood-" + argv[5] +
//       ".csv";

    std::ofstream output_file;
    output_file.open(output_file_name);

    for (i=0; i<N; ++i) {
      output_file << data_points[i];
    }
    output_file.close();

//     parameters_nominal params = parameters_nominal();
//     GaussianInterpolator GP_prior = GaussianInterpolator(points_for_integration,
// 							 points_for_kriging,
// 							 params);
//     // for (unsigned i=0; i<points_for_kriging.size(); ++i) {
//     //   for (unsigned j=0; j<points_for_kriging.size(); ++j) {
//     // 	std::cout << gsl_matrix_get(GP_prior.Cinv, i,j) << " ";
//     //   }
//     //   std::cout << std::endl;
//     // }

//     double integral = 0;
//     for (unsigned i=0; i<points_for_integration.size(); ++i) {
//       double add = GP_prior(points_for_integration[i]);
//       integral = integral +
// 	add;
//       if (add < 0) {
// 	std::cout << points_for_integration[i] << std::endl;
//       }
//     }
//     std::cout << "Integral = " << integral << std::endl;

//     // MultivariateNormal mvtnorm = MultivariateNormal();
//     // std::cout << mvtnorm.dmvnorm(N,y,mu,C) << std::endl;

//     nlopt::opt opt(nlopt::LN_NELDERMEAD, 32);
//     //    std::vector<double> lb =

//     // std::vector<double> x = params.as_vector();
//     // std::cout << optimization_wrapper(x, NULL, &points_for_kriging) << std::endl;

//     gsl_rng_free(r_ptr_local);


    return 0;
}
