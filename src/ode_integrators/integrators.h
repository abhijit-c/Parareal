#ifndef INTEGRATORS_H_INCLUDED
#define INTEGRATORS_H_INCLUDED

#include <functional>
#include <iostream>
#include "../Eigen/Dense"

class ode_system; class time_stepper;

// Functional typedefs

// Functions describing ode system.
typedef std::function<int(double,Eigen::VectorXd &,Eigen::VectorXd &)> ode_rhs;
typedef std::function<int(double,Eigen::VectorXd &,Eigen::MatrixXd &)> ode_jac;

// Functions describing temporal integrator.
typedef std::function<int(ode_system&,double,Eigen::VectorXd &)> ode_intg;
typedef std::function<int(ode_system&,double,Eigen::MatrixXd &)> ode_intg_allt;

class ode_system
{
  public:
    int dimension;
    double t_init, t_final;
    Eigen::VectorXd y0;
    ode_rhs f;
    ode_jac J;
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

// Euler Methods
int forward_euler(ode_system&, double, Eigen::VectorXd &);
int forward_euler_allt(ode_system&, double, Eigen::MatrixXd &);

// Runge-Kutta Methods.
int rk2(ode_system&, double, Eigen::VectorXd &);
int rk2_allt(ode_system&, double, Eigen::MatrixXd &);
int rk4(ode_system&, double, Eigen::VectorXd &);
int rk4_allt(ode_system&, double, Eigen::MatrixXd &);

// Parareal Method.
int parareal(ode_system&, time_stepper, time_stepper, int, Eigen::MatrixXd &);
int pipelined_parareal(ode_system&, time_stepper, time_stepper, int, Eigen::MatrixXd &);


#endif
