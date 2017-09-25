#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "2DBrownianMotionPath.hpp"

using namespace std;

int main() 
{
  unsigned order = 5e6;
  double sigma_x = 1.0;
  double sigma_y = 1.0;
  double rho_base = 0.6;
  double rho = rho_base;
  double x_0 = 0;
  double y_0 = 0;
  double t = 1;
  unsigned number_data_sets = 500;
  unsigned number_observations_per_data_set = 128;

#pragma omp parallel for default(shared)
  for (unsigned i=0; i<number_data_sets; ++i) {

    std::ofstream path_file;
    std::string file_name = "data/data-set-" + std::to_string(i+1) + ".csv";
    path_file.open(file_name);
    path_file << std::fixed << std::setprecision(32);
    
      // header
      path_file << "sigma_x, sigma_y, rho, x_0, y_0, t, a, x_T, b, c, y_T, d\n";
      
      for (unsigned j=0; j<number_observations_per_data_set; ++j) {

  	// long unsigned seed = (i-1)*number_observations_per_data_set + j + 1 + 10;
	long unsigned seed = i*number_observations_per_data_set + j;

  	BrownianMotion BM = BrownianMotion(seed,
					   order,
  					   rho,
  					   sigma_x,
  					   sigma_y,
  					   x_0,
  					   y_0,
  					   t);
  	path_file << sigma_x << ","
  		  << sigma_y << ","
  		  << rho << ","
  		  << x_0 << ","
  		  << y_0 << ","
  		  << t << ","
  		  << BM.get_a() << ","
  		  << BM.get_x_T() << ","
  		  << BM.get_b() << ","
  		  << BM.get_c() << ","
  		  << BM.get_y_T() << ","
		  << BM.get_d() << "\n";

	std::cout << "seed = " << BM.get_seed() << "\n";
      }
      std::cout << std::endl;
      path_file.close();
  }
  return 0;
}
