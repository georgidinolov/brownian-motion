#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include <iostream>
#include <vector>
#include "1DBrownianMotionPath.hpp"

BrownianMotion::BrownianMotion()
  : order_(64),
    sigma_x_(1.0),
    x_0_(0.0),
    t_(1.0)
{
  generate_path();
}

BrownianMotion::BrownianMotion(unsigned order,
			       double sigma_x,
			       double x_0,
			       double t)
  : order_(order),
    sigma_x_(sigma_x),
    x_0_(x_0),
    t_(t)
{
  generate_path();
}

double BrownianMotion::get_a() const
{
  return a_;
}

double BrownianMotion::get_b() const
{
  return b_;
}

double BrownianMotion::get_sigma() const
{
  return sigma_x_;
}
double BrownianMotion::get_x_0() const
{
  return x_0_;
}
double BrownianMotion::get_x_T() const
{
  return path_[order_];
}
double BrownianMotion::get_t() const
{
  return t_;
 }
double BrownianMotion::get_order() const
{
  return order_;
}

const std::vector<double>& BrownianMotion::get_path() const 
{
  return path_;
}

void BrownianMotion::generate_path()
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);

  // std::default_random_engine generator;
  
  std::normal_distribution<double> rnorm(0.0, 1.0);

  double dt = t_/order_;
  path_ = std::vector<double>(order_+1);
  a_ = x_0_;
  b_ = x_0_;

  path_[0] = x_0_;
  for (unsigned i=1; i<order_+1; ++i) {
    double new_element = 
      path_[i-1] + std::sqrt(dt) * sigma_x_ * rnorm(generator);
    
    if (new_element >= b_) {
      b_ = new_element;
    }

    if (new_element <= a_) {
      a_ = new_element;
    }

    path_[i] = new_element;
  }
}
