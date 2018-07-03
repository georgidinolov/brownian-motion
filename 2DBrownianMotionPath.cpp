#include <algorithm>
#include <chrono>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <random>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include "2DBrownianMotionPath.hpp"


BrownianMotion::BrownianMotion()
  : order_(8),
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

BrownianMotion::BrownianMotion(const BrownianMotion& rhs)
  : order_(rhs.order_),
    rho_(rhs.rho_),
    sigma_x_(rhs.sigma_x_),
    sigma_y_(rhs.sigma_y_),
    x_0_(rhs.x_0_),
    y_0_(rhs.y_0_),
    t_(rhs.t_),
    seed_(rhs.seed_),
    path_(rhs.path_),
    a_(rhs.a_),
    b_(rhs.b_),
    c_(rhs.c_),
    d_(rhs.d_)
{}
    
BrownianMotion& BrownianMotion::operator=(const BrownianMotion& rhs)
{
  if (this==&rhs) {
    return *this;
  } else {
    order_ = rhs.get_order();
    rho_ = rhs.get_rho();
    sigma_x_ = rhs.get_sigma_x();
    sigma_y_ = rhs.get_sigma_y();
    x_0_ = rhs.get_x_0();
    y_0_ = rhs.get_y_0();
    t_ = rhs.get_t();
    seed_ = rhs.get_seed();
    path_ = rhs.get_path();
    a_ = rhs.get_a();
    b_ = rhs.get_b();
    c_ = rhs.get_c();
    d_ = rhs.get_d();

    return *this;
  }
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

std::vector<std::vector<double>> BrownianMotion::get_path() const 
{
  return path_;
}

void BrownianMotion::set_path(const std::vector<std::vector<double>>& new_path)
{
  path_ = new_path;
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

void BrownianMotion::print_point() const 
{
  printf("rho=%f\nsigma_x=%f\nsigma_y=%f\nx_0=%f\ny_0=%f\nx_T=%f\ny_T=%f\nt=%f\na=%f\nb=%f\nc=%f\nd=%f\n",
	 rho_,
	 sigma_x_,
	 sigma_y_,
	 x_0_,
	 y_0_,
	 get_x_T(),
	 get_y_T(), //path_[1][path_[1].size()],
	 t_,
	 a_,
	 b_,
	 c_,
	 d_);
}

std::ostream& operator<<(std::ostream& os,
			 const BrownianMotion& current_BM) 
{
  os << "rho=" << std::scientific << std::setprecision(14) << current_BM.rho_ << ";\n";
  os << "sigma_x=" << std::scientific << std::setprecision(14) << current_BM.sigma_x_ << ";\n";
  os << "sigma_y=" << std::scientific << std::setprecision(14) << current_BM.sigma_y_ << ";\n";
  os << "x_0=" << std::scientific << std::setprecision(14) << current_BM.x_0_ << ";\n";
  os << "y_0=" << std::scientific << std::setprecision(14) << current_BM.y_0_ << ";\n";
  os << "x_T=" << std::scientific << std::setprecision(14) << current_BM.get_x_T() << ";\n";
  os << "y_T=" << std::scientific << std::setprecision(14) << current_BM.get_y_T() << ";\n";
  os << "t=" <<std::scientific << std::setprecision(14) << current_BM.t_ << ";\n";
  os << "a=" << std::scientific << std::setprecision(14) << current_BM.a_ << ";\n";
  os << "b=" << std::scientific << std::setprecision(14) << current_BM.b_ << ";\n";
  os << "c=" << std::scientific << std::setprecision(14) << current_BM.c_ << ";\n";
  os << "d=" << std::scientific << std::setprecision(14) << current_BM.d_ << ";\n";
  return os;
}

std::istream& operator>>(std::istream& is,
			 BrownianMotion& current_BM) 
{
  std::string name;
  std::string value;
  // rho_
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  current_BM.rho_ = std::stod(value);
  // sigma_x_
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  current_BM.sigma_x_ = std::stod(value);
  // sigma_y_
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  current_BM.sigma_y_ = std::stod(value);
  // x_0_
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  current_BM.x_0_ = std::stod(value);
  // y_0_
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  current_BM.y_0_ = std::stod(value);
  // x_T
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  std::vector<double>& xmoves = current_BM.path_[0];
  xmoves[current_BM.order_] = std::stod(value);
  // y_T
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  std::vector<double>& ymoves = current_BM.path_[1];
  ymoves[current_BM.order_] = std::stod(value);
  // t_
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  current_BM.t_ = std::stod(value);
  // a_
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  current_BM.a_ = std::stod(value);
  // b_
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  current_BM.b_ = std::stod(value);
  // c_
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  current_BM.c_ = std::stod(value);
  // d_
  std::getline(is, name, '=');
  std::getline(is, value, ';');
  current_BM.d_ = std::stod(value);
  return is;
}

/////////////////////////////////////////////////////

BrownianMotionWithDrift::BrownianMotionWithDrift()
  : BrownianMotion(),
    mu_x_(0.0),
    mu_y_(0.0)
{}
