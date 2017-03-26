#include <algorithm>
#include <chrono>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <random>
#include <iostream>
#include <vector>
#include "2DBrownianMotionPath.hpp"

BrownianMotion::BrownianMotion()
  : order_(64),
    rho_(0.0),
    sigma_x_(1.0),
    sigma_y_(1.0),
    x_0_(0.0),
    y_0_(0.0),
    t_(1.0)
{
  generate_path();
}

BrownianMotion::BrownianMotion(unsigned order,
			       double rho,
			       double sigma_x,
			       double sigma_y,
			       double x_0,
			       double y_0,
			       double t)
  : order_(order),
    rho_(rho),
    sigma_x_(sigma_x),
    sigma_y_(sigma_y),
    x_0_(x_0),
    y_0_(y_0),
    t_(t)
{
  generate_path();
}

BrownianMotion::BrownianMotion(long unsigned seed,
			       unsigned order,
			       double rho,
			       double sigma_x,
			       double sigma_y,
			       double x_0,
			       double y_0,
			       double t)
  : order_(order),
    rho_(rho),
    sigma_x_(sigma_x),
    sigma_y_(sigma_y),
    x_0_(x_0),
    y_0_(y_0),
    t_(t),
    seed_(seed)
{
  generate_path(seed_);
}

double BrownianMotion::get_a() const
{
  return a_;
}

double BrownianMotion::get_b() const
{
  return b_;
}

double BrownianMotion::get_c() const
{
  return c_;
}

double BrownianMotion::get_d() const
{
  return d_;
}

double BrownianMotion::get_sigma_x() const
{
  return sigma_x_;
}

double BrownianMotion::get_sigma_y() const
{
  return sigma_y_;
}

double BrownianMotion::get_rho() const
{
  return rho_;
}

double BrownianMotion::get_x_0() const
{
  return x_0_;
}

double BrownianMotion::get_y_0() const
{
  return y_0_;
}

double BrownianMotion::get_x_T() const
{
  std::vector<double> xmoves = path_[0];
  return xmoves[order_];
}

double BrownianMotion::get_y_T() const
{
  std::vector<double> ymoves = path_[1];
  return ymoves[order_];
}

double BrownianMotion::get_t() const
{
  return t_;
 }
double BrownianMotion::get_order() const
{
  return order_;
}
long unsigned BrownianMotion::get_seed() const
{
  return seed_;
}

const std::vector<std::vector<double>>& BrownianMotion::get_path() const 
{
  return path_;
}

void BrownianMotion::generate_path()
{
  seed_ = std::chrono::system_clock::now().time_since_epoch().count();
  generate_path(seed_);
}

void BrownianMotion::generate_path(unsigned long int seed)
{
  const gsl_rng_type * T;
  gsl_rng * r;
  T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);

  double dt = t_/order_;
  std::vector<double> x_path = std::vector<double>(order_+1);
  std::vector<double> y_path = std::vector<double>(order_+1);
  a_ = x_0_;
  b_ = x_0_;
  c_ = y_0_;
  d_ = y_0_;

  x_path[0] = x_0_;
  y_path[0] = y_0_;
  for (unsigned i=1; i<order_+1; ++i) {
    double new_y = 
      y_path[i-1] + std::sqrt(dt) * sigma_y_ * gsl_ran_gaussian(r, 1.0);

    double new_x = 
      x_path[i-1] + sigma_x_/sigma_y_*rho_*(new_y-y_path[i-1]) + 
      std::sqrt((1-std::pow(rho_,2))*dt) * sigma_x_ * gsl_ran_gaussian(r, 1.0);

    if (new_x >= b_) {
      b_ = new_x;
    }

    if (new_x <= a_) {
      a_ = new_x;
    }

    if (new_y >= d_) {
      d_ = new_y;
    }

    if (new_y <= c_) {
      c_ = new_y;
    }

    x_path[i] = new_x;
    y_path[i] = new_y;
  }

  path_ = std::vector<std::vector<double>> (2);
  path_[0] = x_path;
  path_[1] = y_path;

  gsl_rng_free(r);
}
