#ifndef INTEGRATORS_H_INCLUDED
#define INTEGRATORS_H_INCLUDED

#include <functional>
#include "../Eigen/Dense"

class ode_system; class time_stepper;

// int ode(double t, Eigen::VectorXd y, Eigen::VectorXd &dydt)
typedef std::function<int(double,Eigen::VectorXd &,Eigen::VectorXd &)> ode_rhs;
// int integrator(ode_system sys, double dt, Eigen::VectorXd &yf)
typedef std::function<int(ode_system&,double,Eigen::VectorXd &)> ode_intg;
// int step_integrator(ode_system sys, double dt, Eigen::MatrixXd &yf)
typedef std::function<int(ode_system&,double,Eigen::MatrixXd &)> ode_intg_allt;

class ode_system
{
  public:
    int dimension;
    double t_init, t_final;
    Eigen::VectorXd y0;
    ode_rhs f;
    int num_steps(double dt)
    {
      return 1 + (int) ceil( (t_final-t_init)/dt - 1/2);
    }
};

class time_stepper
{
  public:
    double dt;
    ode_intg F;
    ode_intg_allt F_allt;
    int integrate(ode_system &sys, Eigen::VectorXd &yf)
    {
      return F(sys, dt, yf);
    }
    int integrate_allt(ode_system &sys, Eigen::MatrixXd &steps)
    {
      F_allt(sys, dt, steps);
      return 0;
    }
};

// Explicit Integrators.

// Euler Methods
int forward_euler(ode_system&, double, Eigen::VectorXd &);
int forward_euler_allt(ode_system&, double, Eigen::MatrixXd &);

// Runge-Kutta Methods.
int rk4(ode_system&, double, Eigen::VectorXd &);
int rk4_allt(ode_system&, double, Eigen::MatrixXd &);

// Parareal Method.
int parareal(ode_system&, time_stepper, time_stepper, int, Eigen::MatrixXd &);

#endif
