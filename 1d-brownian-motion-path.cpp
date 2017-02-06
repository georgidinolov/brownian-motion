#include <iostream>
#include <fstream>
#include <vector>
#include "1DBrownianMotionPath.hpp"

using namespace std;

int main() 
{
  BrownianMotion* BM = new BrownianMotion();
  delete BM;

  unsigned order = 2048;
  double sigma_x = 0.1;
  double x_0 = 0;
  double t = 10;

  BM = new BrownianMotion(order,sigma_x,x_0,t);

  cout << "a=" << BM->get_a() << endl;
  cout << "b=" << BM->get_b() << endl;

  ofstream path_file;
  path_file.open("/home/gdinolov/Research/PDE-solvers/src/brownian-motion/path.csv");
  // header
  path_file << "path, sigma, x_0, t, a, b \n";
  std::vector<double> path = BM->get_path();
  unsigned size = path.size();
  
  // filling it up
  for (unsigned i=0; i<size; ++i) {
    if (i==0) {
      path_file << path[i] 
		<< "," << BM->get_sigma() 
		<< "," << BM->get_x_0() 
		<< "," << BM->get_t()
		<< "," << BM->get_a()
		<< "," << BM->get_b() << "\n";
    } else {
      path_file << path[i] << "\n";
    }
  }
  path_file.close();
  return 0;
}
