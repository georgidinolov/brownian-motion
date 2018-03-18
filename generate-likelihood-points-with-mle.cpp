#include <algorithm>
#include "src/gaussian-interpolator/GaussianInterpolator.hpp"
#include "src/multivariate-normal/MultivariateNormal.hpp"
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

double log_likelihood(const std::vector<double> &x,
		      void * data) {
  // // //
  BrownianMotion * lp =
    reinterpret_cast<BrownianMotion*> (data);
  // // //

  double sigma_x = x[0];
  double sigma_y = x[1];
  double rho = x[2];

  double t = lp->get_t();
  // // //
  double out = -log(2*M_PI*1.0*sigma_x*sigma_y*std::sqrt(1-rho*rho)*t) +
    -1/(2*t*(1-rho*rho)) * (
  			    std::pow(lp->get_x_T() - lp->get_x_0(), 2)/std::pow(sigma_x,2.0) +
  			    std::pow(lp->get_y_T() - lp->get_y_0(), 2)/std::pow(sigma_y,2.0) -
  			    2*rho*(lp->get_x_T() - lp->get_x_0())*(lp->get_y_T() - lp->get_y_0())/
			    (sigma_x * sigma_y)
  			    ) ;

  return out;
}

double log_likelihood_transformed(const std::vector<double> &x,
				  void * data) {
  // // //
  double logit_sigma_y_tilde = x[0];
  double log_t = x[1];
  double logit_rho = x[2];

  double sigma_y_tilde = gp_logit_inv(logit_sigma_y_tilde);
  double t = exp(log_t);
  double rho = gp_logit_inv(logit_rho)*2.0 - 1.0;

  std::vector<double> nominal_input = std::vector<double> {sigma_y_tilde,
							   t,
							   rho};
  // // //
  double out = log_likelihood(nominal_input,
			      data);

  // ADDING DERIV OF TRANSFORMATION
  out = out +
    log(sigma_y_tilde) + log(1-sigma_y_tilde) +
    log_t +
    log(2.0) + log(std::abs( (rho+1.0)/2.0 )) + log(std::abs( 1.0 - (rho+1.0)/2.0 ));

  return out;
}

double numeric_log_deriv_rho(const std::vector<double> &x,
			     void * data,
			     double drho) {
  std::vector<double> my_x = x;
  my_x[2] = x[2] + drho;
  double out_plus = log_likelihood(my_x, data);

  my_x[2] = x[2] - drho;
  double out_minus = log_likelihood(my_x, data);

  return (out_plus - out_minus)/(2*drho);
}

double log_deriv_rho(const std::vector<double> &x,
		     void * data) {

  likelihood_point * lp =
    reinterpret_cast<likelihood_point*> (data);

  double sigma_y_tilde = x[0];
  double t = x[1];
  double rho = x[2];

  double S = (lp->x_t_tilde - lp->x_0_tilde)*(lp->x_t_tilde - lp->x_0_tilde)/std::pow(1.0,2.0) +
    (lp->y_t_tilde - lp->y_0_tilde)*(lp->y_t_tilde - lp->y_0_tilde)/std::pow(sigma_y_tilde,2.0) -
    2*rho*(lp->x_t_tilde - lp->x_0_tilde)*(lp->y_t_tilde - lp->y_0_tilde)/(1.0 * sigma_y_tilde);

  double S_rho =
    -2*(lp->x_t_tilde - lp->x_0_tilde)*(lp->y_t_tilde - lp->y_0_tilde)/(1.0 * sigma_y_tilde);

  double out = rho/(1-std::pow(rho,2.0)) -
    rho/(t*std::pow((1-rho*rho),2.0))*S -
    1.0/(2.0*t*(1-rho*rho))*S_rho;

  return out;
}


double numeric_log_deriv_t(const std::vector<double> &x,
			     void * data,
			     double dt) {
  std::vector<double> my_x = x;
  my_x[1] = x[1] + dt;
  double out_plus = log_likelihood(my_x, data);

  my_x[1] = x[1] - dt;
  double out_minus = log_likelihood(my_x, data);

  return (out_plus - out_minus)/(2*dt);
}

double log_deriv_t(const std::vector<double> &x,
		     void * data) {

  likelihood_point * lp =
    reinterpret_cast<likelihood_point*> (data);

  double sigma_y_tilde = x[0];
  double t = x[1];
  double rho = x[2];

  double S = (lp->x_t_tilde - lp->x_0_tilde)*(lp->x_t_tilde - lp->x_0_tilde)/std::pow(1.0,2.0) +
    (lp->y_t_tilde - lp->y_0_tilde)*(lp->y_t_tilde - lp->y_0_tilde)/std::pow(sigma_y_tilde,2.0) -
    2*rho*(lp->x_t_tilde - lp->x_0_tilde)*(lp->y_t_tilde - lp->y_0_tilde)/(1.0 * sigma_y_tilde);

  double out = -1.0/t + S/(2.0*std::pow(t,2.0)*(1.0-std::pow(rho,2)));

  return out;
}

double numeric_log_deriv_sigma(const std::vector<double> &x,
			       void * data,
			       double dsigma) {
  std::vector<double> my_x = x;
  my_x[0] = x[0] + dsigma;
  double out_plus = log_likelihood(my_x, data);

  my_x[0] = x[0] - dsigma;
  double out_minus = log_likelihood(my_x, data);

  return (out_plus - out_minus)/(2*dsigma);
}

double log_deriv_sigma(const std::vector<double> &x,
		       void * data) {

  likelihood_point * lp =
    reinterpret_cast<likelihood_point*> (data);

  double sigma_y_tilde = x[0];
  double t = x[1];
  double rho = x[2];

  double S = (lp->x_t_tilde - lp->x_0_tilde)*(lp->x_t_tilde - lp->x_0_tilde)/std::pow(1.0,2.0) +
    (lp->y_t_tilde - lp->y_0_tilde)*(lp->y_t_tilde - lp->y_0_tilde)/std::pow(sigma_y_tilde,2.0) -
    2*rho*(lp->x_t_tilde - lp->x_0_tilde)*(lp->y_t_tilde - lp->y_0_tilde)/(1.0 * sigma_y_tilde);

  double S_sigma =
    -2.0*(lp->y_t_tilde - lp->y_0_tilde)*(lp->y_t_tilde - lp->y_0_tilde)/std::pow(sigma_y_tilde,3.0) -
    -1.0*2*rho*(lp->x_t_tilde - lp->x_0_tilde)*(lp->y_t_tilde - lp->y_0_tilde)/(1.0 * std::pow(sigma_y_tilde,2.0));

  double out = -1.0/sigma_y_tilde - 1.0/(2.0*t*(1-std::pow(rho,2.0))) * S_sigma;

  return out;
}

double numeric_log_deriv_sigma_sigma(const std::vector<double> &x,
				     void * data,
				     double dsigma) {
  std::vector<double> my_x = x;
  my_x[0] = x[0] + dsigma;
  double out_plus = log_deriv_sigma(my_x, data);

  my_x = x;
  my_x[0] = x[0] - dsigma;
  double out_minus = log_deriv_sigma(my_x, data);

  return (out_plus - out_minus)/(2*dsigma);
}

double numeric_log_deriv_sigma_t(const std::vector<double> &x,
				   void * data,
				   double dt) {
  std::vector<double> my_x = x;
  my_x[1] = x[1] + dt;
  double out_plus = log_deriv_sigma(my_x, data);

  my_x = x;
  my_x[1] = x[1] - dt;
  double out_minus = log_deriv_sigma(my_x, data);

  return (out_plus - out_minus)/(2*dt);
}

double numeric_log_deriv_sigma_rho(const std::vector<double> &x,
				   void * data,
				   double drho) {
  std::vector<double> my_x = x;
  my_x[2] = x[2] + drho;
  double out_plus = log_deriv_sigma(my_x, data);

  my_x = x;
  my_x[2] = x[2] - drho;
  double out_minus = log_deriv_sigma(my_x, data);

  return (out_plus - out_minus)/(2*drho);
}

double numeric_log_deriv_t_t(const std::vector<double> &x,
				   void * data,
				   double dt) {
  std::vector<double> my_x = x;
  my_x[1] = x[1] + dt;
  double out_plus = log_deriv_t(my_x, data);

  my_x = x;
  my_x[1] = x[1] - dt;
  double out_minus = log_deriv_t(my_x, data);

  return (out_plus - out_minus)/(2*dt);
}

double numeric_log_deriv_t_rho(const std::vector<double> &x,
				   void * data,
				   double drho) {
  std::vector<double> my_x = x;
  my_x[2] = x[2] + drho;
  double out_plus = log_deriv_t(my_x, data);

  my_x = x;
  my_x[2] = x[2] - drho;
  double out_minus = log_deriv_t(my_x, data);

  return (out_plus - out_minus)/(2*drho);
}

double numeric_log_deriv_rho_rho(const std::vector<double> &x,
				   void * data,
				   double drho) {
  std::vector<double> my_x = x;
  my_x[2] = x[2] + drho;
  double out_plus = log_deriv_rho(my_x, data);

  my_x = x;
  my_x[2] = x[2] - drho;
  double out_minus = log_deriv_rho(my_x, data);

  return (out_plus - out_minus)/(2*drho);
}

double myfunc(const std::vector<double> &x,
	      std::vector<double> &grad,
	      void * data) {

  // if (!grad.empty()) {
  //   grad[0] = log_deriv_sigma(x, data);
  //   grad[1] = log_deriv_t(x, data);
  //   grad[2] = log_deriv_rho(x, data);
  // }

  std::vector<BrownianMotion> * lps =
    reinterpret_cast<std::vector<BrownianMotion>*> (data);

  double out = 0.0;
  for (unsigned i=0; i<lps->size(); ++i) {
    BrownianMotion lp = (*lps)[i];
    out = out + log_likelihood(x, &lp);
  }
  return out;
}

int main(int argc, char *argv[]) {
  if (argc < 4 || argc > 4) {
    printf("You must provide input\n");
    printf("The input is: \n\ninput number data points;\ninput file name (as BrownianMotions);\noutput file name;\n");
    exit(0);
  }

  unsigned N = std::stoi(argv[1]);
  unsigned m = 80;
  std::ifstream data_file(argv[2]);
  std::string output_file_name = argv[3];

  unsigned i;
  std::vector<BrownianMotion> data_points = std::vector<BrownianMotion> (N);

  static int counter = 0;
  static gsl_rng* r_ptr_threadprivate;
#pragma omp threadprivate(counter, r_ptr_threadprivate)
  omp_set_dynamic(0);
  omp_set_num_threads(1);

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

  for (unsigned i=0; i<N; ++i) {
    data_file >>  data_points[i];
  }
  std::vector<likelihood_point> emulator_points (N*m);

  std::vector<double> lb = {0.001, 0.001, -0.95};
  std::vector<double> ub = {HUGE_VAL, HUGE_VAL,   0.95};
  std::vector<double> MSE = std::vector<double> (N);

  auto t1 = std::chrono::high_resolution_clock::now();
#pragma omp parallel default(none) private(i) shared(data_points, N, seed_init, lb, ub, MSE, m, emulator_points)
    {
#pragma omp for
      for (i=0; i<N; ++i) {
	unsigned ss = 4;
	std::vector<BrownianMotion> current_lps = std::vector<BrownianMotion> (ss);
	if (i < N- ss-1) {
	  for (unsigned j=0; j<ss; ++j) {
	    current_lps[j] = data_points[i+j];
	  }
	} else {
	  for (unsigned j=0; j<ss; ++j) {
	    current_lps[j] = data_points[i-j];
	  }
	}
	  // for (likelihood_point lp : current_lps) {
	  //   lp.print_point();
	  // }

	  nlopt::opt opt(nlopt::LN_NELDERMEAD, 3);
	  opt.set_lower_bounds(lb);
	  opt.set_upper_bounds(ub);
	  opt.set_ftol_rel(0.0001);
	  double maxf;

	  // CONFIGURATION 1
	  opt.set_max_objective(myfunc,
				&current_lps);
	  std::vector<double> optimal_params = std::vector<double> (3);
	  optimal_params[0] = 1.0;
	  optimal_params[1] = 1.0;
	  optimal_params[2] = 0.0;

	  opt.optimize(optimal_params, maxf);

	  gsl_matrix * MLE_cov = gsl_matrix_alloc(2,2);
	  gsl_matrix * sample_cov = gsl_matrix_alloc(2,2);
	  gsl_matrix_set(MLE_cov, 0, 0, (ss + 2 + 1)*std::pow(optimal_params[0],2));
	  gsl_matrix_set(MLE_cov, 0, 1, (ss + 2 + 1)*optimal_params[0]*optimal_params[1]*optimal_params[2]);
	  gsl_matrix_set(MLE_cov, 1, 0, (ss + 2 + 1)*optimal_params[0]*optimal_params[1]*optimal_params[2]);
	  gsl_matrix_set(MLE_cov, 1, 1, (ss + 2 + 1)*std::pow(optimal_params[1],2));
	//   double dsigma = 0.001;
	//   double dt = 0.001;
	//   double drho = 0.001;

	//   double sd_sigma = 0.5;
	//   if (numeric_log_deriv_sigma_sigma(optimal_params, &current_lps[0], dsigma) < 0) {
	//     sd_sigma = sqrt(1/-numeric_log_deriv_sigma_sigma(optimal_params, &current_lps[0], dsigma));
	//   }
	//   double sd_t = sqrt(1/-numeric_log_deriv_t_t(optimal_params, &current_lps[0], dt))*10;
	//   double sd_rho = 1.0;
	//   if (numeric_log_deriv_rho_rho(optimal_params, &current_lps[0], drho) < 0) {
	//     sd_rho = sqrt(1/-numeric_log_deriv_rho_rho(optimal_params, &current_lps[0], drho));
	//   }
	  MultivariateNormal mvtnorm = MultivariateNormal();
	  for (unsigned j=0; j<m; ++j) {
	    mvtnorm.rinvwishart(r_ptr_threadprivate,
				2,
				ss,
				MLE_cov,
				sample_cov);
	    // //
	    double sigma_x = sqrt( gsl_matrix_get(sample_cov, 0,0) );
	    double sigma_y = sqrt( gsl_matrix_get(sample_cov, 1,1) );
	    double rho = gsl_matrix_get(sample_cov, 0, 1)/(sigma_x*sigma_y);

	    // sigma_x = optimal_params[0];
	    // sigma_y = optimal_params[1];

	    likelihood_point current_lp = likelihood_point(rho,sigma_x,sigma_y,
							   data_points[i].get_x_0(),
							   data_points[i].get_y_0(),
							   data_points[i].get_x_T(),
							   data_points[i].get_y_T(),
							   data_points[i].get_t(),
							   data_points[i].get_a(),
							   data_points[i].get_b(),
							   data_points[i].get_c(),
							   data_points[i].get_d());
	    while ((current_lp.sigma_y_tilde < 0.40) ||
		   (current_lp.t_tilde < 0.30) ||
		   (std::abs(current_lp.rho) > 0.95)) {
	      mvtnorm.rinvwishart(r_ptr_threadprivate,
				  2,
				  ss,
				  MLE_cov,
				  sample_cov);
	      // //
	      sigma_x = sqrt( gsl_matrix_get(sample_cov, 0,0) );
	      sigma_y = sqrt( gsl_matrix_get(sample_cov, 1,1) );
	      rho = gsl_matrix_get(sample_cov, 0, 1)/(sigma_x*sigma_y);

	      current_lp = likelihood_point(rho,sigma_x,sigma_y,
					    data_points[i].get_x_0(),
					    data_points[i].get_y_0(),
					    data_points[i].get_x_T(),
					    data_points[i].get_y_T(),
					    data_points[i].get_t(),
					    data_points[i].get_a(),
					    data_points[i].get_b(),
					    data_points[i].get_c(),
					    data_points[i].get_d());
	    }

	    current_lp.log_likelihood = 
	      log_likelihood(std::vector<double> {sigma_x,sigma_y,rho},
			     &data_points[i]) +
	      log(data_points[i].get_b() - data_points[i].get_a()) +
	      log(data_points[i].get_d() - data_points[i].get_c());
	    
	    emulator_points[i*m + j] = current_lp;
	  }

	  gsl_matrix_free(MLE_cov);
	  gsl_matrix_free(sample_cov);

	//   gsl_matrix_set(repeated_sampling_cov, 0,0, -numeric_log_deriv_sigma_sigma(optimal_params, &current_lps[0], dsigma));
	//   gsl_matrix_set(repeated_sampling_cov, 0,1, -numeric_log_deriv_sigma_t(optimal_params, &current_lps[0], dt));
	//   gsl_matrix_set(repeated_sampling_cov, 0,2, -numeric_log_deriv_sigma_rho(optimal_params, &current_lps[0], drho));
	//   //
	//   gsl_matrix_set(repeated_sampling_cov, 1,1, -numeric_log_deriv_t_t(optimal_params, &current_lps[0], dt));
	//   gsl_matrix_set(repeated_sampling_cov, 1,2, -numeric_log_deriv_t_rho(optimal_params, &current_lps[0], drho));
	//   //
	//   gsl_matrix_set(repeated_sampling_cov, 2,2, -numeric_log_deriv_rho_rho(optimal_params, &current_lps[0], drho));

	//   for (unsigned j=0; j<3; ++j) {
	//     for (unsigned k=j; k<3; ++k) {
	//       gsl_matrix_set(repeated_sampling_cov, k,j,
	// 		     gsl_matrix_get(repeated_sampling_cov, j,k)) ;
	//     }
	//   }

	  printf("Thread %d with address and data point %d -- produces likelihood %f: (%f, %f, %f)\n",
		 omp_get_thread_num(),
		 i,
		 maxf,
		 optimal_params[0],
		 optimal_params[1],
		 optimal_params[2]);

	// // // CONFIGURATION 2
	// // current_lps[0] = data_points[i];
	// // current_lps[0].x_0_tilde = data_points[i].y_0_tilde;
	// // current_lps[0].y_0_tilde = data_points[i].x_0_tilde;

	// // current_lps[0].x_t_tilde = data_points[i].y_t_tilde;
	// // current_lps[0].y_t_tilde = data_points[i].x_t_tilde;

	// // opt.set_max_objective(myfunc,
	// // 		      &current_lps[0]);

	// // optimal_params = std::vector<double> (3);
	// // optimal_params[0] = 0.5;
	// // optimal_params[1] = 1.0;
	// // optimal_params[2] = 0.0;
	// // opt.optimize(optimal_params, maxf);

	// // dsigma = 0.001;
	// // dt = 0.001;
	// // drho = 0.001;

	// // sd_sigma = 0.5;
	// // if (numeric_log_deriv_sigma_sigma(optimal_params, &current_lps[0], dsigma) < 0) {
	// //   sd_sigma = sqrt(1/-numeric_log_deriv_sigma_sigma(optimal_params, &current_lps[0], dsigma));
	// // }
	// // sd_t = sqrt(1/-numeric_log_deriv_t_t(optimal_params, &current_lps[0], dt))*10;
	// // sd_rho = 1.0;
	// // if (numeric_log_deriv_rho_rho(optimal_params, &current_lps[0], drho) < 0) {
	// //   sd_rho = sqrt(1/-numeric_log_deriv_rho_rho(optimal_params, &current_lps[0], drho));
	// // }

	// // for (unsigned j=m; j<2*m; ++j) {
	// //   double sigma_sample = gsl_ran_gaussian(r_ptr_threadprivate, sd_sigma) + optimal_params[0];
	// //   while (sigma_sample < 0.0 || sigma_sample > 1.0) {
	// //     sigma_sample = gsl_ran_gaussian(r_ptr_threadprivate, sd_sigma) + optimal_params[0];
	// //   }

	// //   double t_sample = gsl_ran_gaussian(r_ptr_threadprivate, sd_t) + optimal_params[1];
	// //   while (t_sample <= 0.0) {
	// //     t_sample = gsl_ran_gaussian(r_ptr_threadprivate, sd_t) + optimal_params[1];
	// //   }

	// //   double rho_sample = gsl_ran_gaussian(r_ptr_threadprivate, sd_rho) + optimal_params[2];
	// //   while (rho_sample <= -0.95 || rho_sample >= 0.95) {
	// //     rho_sample = gsl_ran_gaussian(r_ptr_threadprivate, sd_rho) + optimal_params[2];
	// //   }
	// //   current_lps[0].sigma_y_tilde = sigma_sample;
	// //   current_lps[0].t_tilde = t_sample;
	// //   current_lps[0].rho = rho_sample;
	// //   emulator_points[i*2*m + j] = current_lps[0];
	// // }

	// // gsl_matrix_set(repeated_sampling_cov, 0,0, -numeric_log_deriv_sigma_sigma(optimal_params, &current_lps[0], dsigma));
	// // gsl_matrix_set(repeated_sampling_cov, 0,1, -numeric_log_deriv_sigma_t(optimal_params, &current_lps[0], dt));
	// // gsl_matrix_set(repeated_sampling_cov, 0,2, -numeric_log_deriv_sigma_rho(optimal_params, &current_lps[0], drho));
	// // //
	// // gsl_matrix_set(repeated_sampling_cov, 1,1, -numeric_log_deriv_t_t(optimal_params, &current_lps[0], dt));
	// // gsl_matrix_set(repeated_sampling_cov, 1,2, -numeric_log_deriv_t_rho(optimal_params, &current_lps[0], drho));
	// // //
	// // gsl_matrix_set(repeated_sampling_cov, 2,2, -numeric_log_deriv_rho_rho(optimal_params, &current_lps[0], drho));

	// // for (unsigned j=0; j<3; ++j) {
	// //   for (unsigned k=j; k<3; ++k) {
	// //     gsl_matrix_set(repeated_sampling_cov, k,j,
	// // 		   gsl_matrix_get(repeated_sampling_cov, j,k)) ;
	// //   }
	// // }


	// // // gsl_linalg_cholesky_decomp(repeated_sampling_cov);
	// // // gsl_linalg_cholesky_invert(repeated_sampling_cov);



	// // MSE[i] = sqrt( std::pow((current_lps[0].sigma_y_tilde - optimal_params[0]), 2) +
	// // 	       std::pow((current_lps[0].t_tilde - optimal_params[1]), 2) +
	// // 	       std::pow((current_lps[0].rho - optimal_params[3]), 2) );

      }
    }
    auto t2 = std::chrono::high_resolution_clock::now();

    double out = 0;
    for (double mse : MSE) {
      out = out + mse;
    }
    std::cout << "MSE = " << out / N << std::endl;

    std::ofstream output_file;
    output_file.open(output_file_name);

    for (likelihood_point lp : emulator_points) {
      output_file << lp;
    }
    output_file.close();

    return 0;
}
