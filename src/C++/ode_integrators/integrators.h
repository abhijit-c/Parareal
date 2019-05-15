#ifndef INTEGRATORS_H_INCLUDED
#define INTEGRATORS_H_INCLUDED

#include <Eigen/Dense>
#include <functional>

class ode_system; class time_stepper;

// int ode(double t, Eigen::VectorXd y, Eigen::VectorXd &dydt)
typedef std::function<int(double, Eigen::VectorXd &, Eigen::VectorXd &)> ode_rhs;
// int integrator(ode_system sys, double dt, Eigen::VectorXd &yf)
typedef std::function<int(ode_system, double, Eigen::VectorXd &)> ode_method;

class ode_system
{
  public:
    int dimension;
    double t_init, t_final;
    Eigen::VectorXd y0;
    ode_rhs f;
};

class time_stepper
{
  public:
    double dt;
    ode_method F;
    int integrate(ode_system sys, Eigen::VectorXd &yf)
    {
      return F(sys, dt, yf);
    }
    int sstep_integrate(ode_system sys, Eigen::VectorXd &yf)
    {
      ode_system step_sys = sys;
      step_sys.t_final = sys.t_init + dt;
      return F(step_sys, dt, yf);
    }
};

// Explicit Integrators.
int forward_euler(ode_system, double, Eigen::VectorXd &);
int rk4(ode_system, double, Eigen::VectorXd &);

// Implicit Integrators.
int backward_euler(ode_system, double, Eigen::VectorXd &);

// Parareal Method.
int parareal(ode_system, time_stepper, time_stepper, Eigen::VectorXd &);

#endif
