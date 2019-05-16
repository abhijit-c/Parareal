#ifndef INTEGRATORS_H_INCLUDED
#define INTEGRATORS_H_INCLUDED

#include <Eigen/Dense>
#include <functional>

class ode_system; class time_stepper;

// int ode(double t, Eigen::VectorXd y, Eigen::VectorXd &dydt)
typedef std::function<int(double,Eigen::VectorXd &,Eigen::VectorXd &)> ode_rhs;
// int integrator(ode_system sys, double dt, Eigen::VectorXd &yf)
typedef std::function<int(ode_system,double,Eigen::VectorXd &)> ode_intg;
// int step_integrator(ode_system sys, double dt, Eigen::MatrixXd &yf)
typedef std::function<int(ode_system,double,Eigen::MatrixXd &)> sode_intg;

class ode_system
{
  public:
    int dimension;
    double t_init, t_final;
    Eigen::VectorXd y0;
    ode_rhs f;
    int num_steps(double dt)
    {
      return 1 + (int) ceil( (t_final+dt/2-t_init)/dt );
    }
};

class time_stepper
{
  public:
    double dt;
    ode_intg F;
    sode_intg F_steps;
    int integrate(ode_system sys, Eigen::VectorXd &yf)
    {
      return F(sys, dt, yf);
    }
    int s_integrate(ode_system sys, Eigen::MatrixXd &steps)
    {
      return F_steps(sys, dt, steps);
    }
};

// Explicit Integrators.

// Euler Methods
int forward_euler(ode_system, double, Eigen::VectorXd &);
int forward_euler_steps(ode_system, double, Eigen::MatrixXd &);

// Runge-Kutta Methods.
int rk4(ode_system, double, Eigen::VectorXd &);

// Implicit Integrators.
int backward_euler(ode_system, double, Eigen::VectorXd &);

// Parareal Method.
int parareal(ode_system, time_stepper, time_stepper, Eigen::MatrixXd &);

#endif
