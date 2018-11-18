#include "2DBrownianMotionPath.hpp"
#include <algorithm>
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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <limits>
#include <omp.h>
#include <stdio.h>
#include <string.h>
#include <vector>

// Generates N OCHL data points to be consumed by MLE code. Assumes
// \mu_x = \mu_y = 0, \sigma_x = \sigma_y = 1, \rho is specified by
// input. Output is formatted for consuption by 2DMLEFiniteElement object.
//
// The random seed is also initialized by the user.

int main(int argc, char *argv[]) {
  if (argc < 6 || argc > 6) {
    printf("You must provide input\n");
    printf("The input is: \n\nnumber data points;\nrho;\nseed;\noutput file name; \noutput file name as Brownian Motion\n");
    exit(0);
  }

  unsigned N = std::stoi(argv[1]);
  double rho = std::stod(argv[2]);
  long unsigned seed_init = std::stoi(argv[3]);
  std::string output_file_name = argv[4];
  std::string output_file_name_as_BM = argv[5];
  unsigned i;
  std::vector<BrownianMotion> data_points = std::vector<BrownianMotion> (N);
  double sigma_x = 1.0;
  double sigma_y = 1.0;

  static int counter = 0;
  static gsl_rng* r_ptr_threadprivate;
#pragma omp threadprivate(counter, r_ptr_threadprivate)
  omp_set_dynamic(0);
  omp_set_num_threads(10);


  gsl_rng * r_ptr_local;
  const gsl_rng_type * Type;
  gsl_rng_env_setup();
  Type = gsl_rng_default;
  r_ptr_local = gsl_rng_alloc(Type);
  gsl_rng_set(r_ptr_local, seed_init);
  
  int tid=0;
#pragma omp parallel default(none) private(tid, i) shared(r_ptr_local)
  {
    tid = omp_get_thread_num();

    r_ptr_threadprivate = gsl_rng_clone(r_ptr_local);
    gsl_rng_set(r_ptr_threadprivate, tid);

    printf("Thread %d: counter %d\n", tid, counter);
  }

  std::vector<long unsigned int> seeds_for_points (N);
  for (i=0; i<N; ++i) {
    seeds_for_points[i] = gsl_rng_get(r_ptr_local);
  }

  auto t1 = std::chrono::high_resolution_clock::now();
#pragma omp parallel default(none) private(i) shared(data_points, N, seeds_for_points, rho, sigma_x,sigma_y)
    {
#pragma omp for
      for (i=0; i<N; ++i) {
	long unsigned seed = seeds_for_points[i];
	gsl_rng_set(r_ptr_threadprivate, seed);

	BrownianMotion BM_current = BrownianMotion(seed,
						   1e6,
						   rho,
						   sigma_x,
						   sigma_y,
						   0.0,
						   0.0,
						   1.0);
	data_points[i] = BM_current;
      }
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "data generation took " 
	      << std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() 
	      << " milliseconds"
	      << std::endl;

    std::ofstream output_file;
    std::ofstream output_file_BM;
    output_file.open(output_file_name);
    output_file_BM.open(output_file_name_as_BM);

    // header 
    output_file << "sigma_x, sigma_y, rho, x_0, y_0, t, a, x_T, b, c, y_T, d\n";

    for (i=0; i<N; ++i) {
      output_file_BM << data_points[i];

      double x_0 = data_points[i].get_x_0();
      double y_0 = data_points[i].get_y_0();
      // 
      double t = data_points[i].get_t();
      //
      double a = data_points[i].get_a();
      double x_T = data_points[i].get_x_T();
      double b = data_points[i].get_b();
      //
      double c = data_points[i].get_c();
      double y_T = data_points[i].get_y_T();
      double d = data_points[i].get_d();

      output_file << std::scientific << std::setprecision(14) << sigma_x << ","
		  << std::scientific << std::setprecision(14) << sigma_y << ","
		  << std::scientific << std::setprecision(14) << rho << ","
		  << std::scientific << std::setprecision(14) << x_0 << ","
		  << std::scientific << std::setprecision(14) << y_0 << ","
		  << std::scientific << std::setprecision(14) << t << ","
		  << std::scientific << std::setprecision(14) << a << ","
		  << std::scientific << std::setprecision(14) << x_T << ","
		  << std::scientific << std::setprecision(14) << b << ","
		  << std::scientific << std::setprecision(14) << c << ","
		  << std::scientific << std::setprecision(14) << y_T << ","
		  << std::scientific << std::setprecision(14) << d;
      if (i<N-1) {
	output_file << "\n";
      }
    }
    output_file.close();
    output_file_BM.close();

    return 0;
}
