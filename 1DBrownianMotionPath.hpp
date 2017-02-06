#include <vector>

class BrownianMotion
{
public:
  BrownianMotion();
  BrownianMotion(unsigned order,
		 double sigma_x,
		 double x_0,
		 double t);

  double get_a() const;
  double get_b() const;
  double get_sigma() const;
  double get_x_0() const;
  double get_x_T() const;
  double get_t() const;
  double get_order() const;
  const std::vector<double>& get_path() const;
  void generate_path();
  
private:
  unsigned order_;
  double sigma_x_;
  double x_0_;
  double t_;
  std::vector<double> path_;
  double a_; // min obtained
  double b_; // max obtained
};
