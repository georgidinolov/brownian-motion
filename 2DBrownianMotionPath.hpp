#include <vector>

class BrownianMotion
{
public:
  BrownianMotion();
  BrownianMotion(unsigned order,
		 double rho,
		 double sigma_x,
		 double sigma_y,
		 double x_0,
		 double y_0,
		 double t);

  BrownianMotion(long unsigned seed,
		 unsigned order,
		 double rho,
		 double sigma_x,
		 double sigma_y,
		 double x_0,
		 double y_0,
		 double t);

  double get_a() const;
  double get_b() const;
  double get_c() const;
  double get_d() const;
  double get_sigma_x() const;
  double get_sigma_y() const;
  double get_rho() const;
  double get_x_0() const;
  double get_x_T() const;
  double get_y_0() const;
  double get_y_T() const;
  double get_t() const;
  double get_order() const;
  long unsigned get_seed() const;
  const std::vector<std::vector<double>>& get_path() const;
  void generate_path();
  void generate_path(unsigned long seed);
  
private:
  unsigned order_;
  double rho_;
  double sigma_x_;
  double sigma_y_;
  double x_0_;
  double y_0_;
  double t_;
  long unsigned seed_;
  std::vector<std::vector<double>> path_;
  double a_; // min obtained x-dir
  double b_; // max obtained x-dir
  double c_; // min obtained y-dir
  double d_; // max obtained y-dir
};
