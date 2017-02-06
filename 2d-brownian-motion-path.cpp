#include <iostream>
#include <fstream>
#include <vector>
#include "2DBrownianMotionPath.hpp"

using namespace std;

int main() 
{
  BrownianMotion* BM = new BrownianMotion();
  delete BM;

  unsigned order = 10000;
  double sigma_x = 1.0;
  double sigma_y = 1.0;
  double rho = 0.4;
  double x_0 = 0;
  double y_0 = 0;
  double t = 1;

  BM = new BrownianMotion(order,
			  rho,
			  sigma_x,
			  sigma_y,
			  x_0,
			  y_0,
			  t);

  cout << "a=" << BM->get_a() << endl;
  cout << "b=" << BM->get_b() << endl;
  cout << "c=" << BM->get_c() << endl;
  cout << "d=" << BM->get_d() << endl;

  std::ofstream path_file ("2d-path.csv");
  if (path_file.is_open()) {
    std::cout << "file is open" << std::endl;
  }
  // header
  path_file << "xpath, ypath, sigma_x, sigma_y, x_0, y_0, t, a, b, c, d \n";
  std::vector<std::vector<double>> path = BM->get_path();
  std::vector<double> xpath = path[0];
  std::vector<double> ypath = path[1];
  unsigned size = path[0].size();
  
  // filling it up
  for (unsigned i=0; i<size; ++i) {
    if (i==0) {
      path_file << xpath[i] 
		<< "," << ypath[i]
		<< "," << BM->get_sigma_x() 
		<< "," << BM->get_sigma_y() 
		<< "," << BM->get_x_0() 
		<< "," << BM->get_y_0() 
		<< "," << BM->get_t()
		<< "," << BM->get_a()
		<< "," << BM->get_b() 
		<< "," << BM->get_c() 
		<< "," << BM->get_d() 
		<< "\n";
    } else {
      path_file << xpath[i] << "," << ypath[i] << "\n";
    }
  }
  path_file.close();
  return 0;
}
