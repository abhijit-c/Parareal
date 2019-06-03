#include "backward_euler.h"

inline int
backward_euler_newton_method(
    const ode_system &sys, 
    double t,
    double dt, 
    const Eigen::VectorXd &y0, 
    Eigen::VectorXd &y)
{
  int D = sys.dimension, k = 0;
  Eigen::VectorXd dydt = y0; y = y0;
  Eigen::MatrixXd J(D,D);
  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(D,D);

  sys.f(t+dt, y, dydt);
  double residual = 0.0;
  do
  { //TODO: Code motion b out of loop to avoid recomputation in while.
    sys.J(t+dt, y, J);
    Eigen::MatrixXd A = I - dt*J;
    Eigen::VectorXd b = dt*dydt + y0 - y;
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);
    y += x; k += 1;
    
    sys.f(t+dt, y, dydt);
  } while ( (dt*dydt+y0-y).norm() > 1e-3 && k < 100);
  return 0;
}

inline int 
backward_euler(ode_system &sys, double dt, Eigen::VectorXd &yf)
{
  yf = sys.y0;
  Eigen::VectorXd y_temp(sys.dimension);
  double t = sys.t_init;
  while (t + dt/2 < sys.t_final)
  {
    y_temp = yf;
    backward_euler_newton_method(sys, t, dt, y_temp, yf);

    t = t + dt;
  }
  return 0;
}

inline int 
backward_euler_allt(ode_system &sys, double dt, Eigen::MatrixXd &steps)
{
  Eigen::VectorXd yf = sys.y0; steps.row(0) = sys.y0;
  Eigen::VectorXd y_temp(sys.dimension);
  double t = sys.t_init; int k = 1;
  while (t + dt/2 < sys.t_final)
  {
    y_temp = yf;
    backward_euler_newton_method(sys, t, dt, y_temp, yf);
    steps.row(k) = yf;

    t = t + dt; k = k + 1;
  }
  return 0;
}
