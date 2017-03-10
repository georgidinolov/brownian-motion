#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include "2DBrownianMotionPath.hpp"

using namespace std;

int main() 
{
  unsigned order = 10000;
  double sigma_x = 1.0;
  double sigma_y = 1.0;
  double rho_base = 0.6;
  double rho = rho_base;
  double x_0 = 0;
  double y_0 = 0;
  double t = 1;
  unsigned number_data_sets = 1;
  unsigned number_observations_per_data_set = 1000;

  for (unsigned i=0+1; i<number_data_sets+1; ++i) {
    // if (i > 500 + 1000) {
    //   rho = -1.0 * rho_base;
    // } else {
    //   rho = rho_base;
    // }

    std::ofstream path_file;
    std::string file_name = "data/data-set-" +
      std::to_string(i+1) +
      ".csv";
    path_file.open(file_name);
    path_file << std::fixed << std::setprecision(32);
    
      // header
      path_file << "sigma_x, sigma_y, rho, x_0, y_0, t, a, x_T, b, c, y_T, d \n";
      
      for (unsigned j=0; j<number_observations_per_data_set; ++j) {

  	long unsigned seed = (i-1)*number_observations_per_data_set + j + 1;

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
  		  << -BM.get_b() << ","
  		  << -BM.get_x_T() << ","
  		  << -BM.get_a() << ","
  		  << BM.get_c() << ","
  		  << BM.get_y_T() << ","
		  << BM.get_d() << "\n";
      }
      
      path_file.close();
  }
  return 0;
}
