#ifndef INTEGRATORS_H_INCLUDED
#define INTEGRATORS_H_INCLUDED

#include <functional>

/* integrator function declaration syntax.
 * first argument should be ode system, which should be a function object:
 *    double *ode(int dimension, double t, double *y, double *dydt)
 */

class ode_system; class time_stepper;

typedef std::function<int(double, double *, double *)> ode_rhs;
typedef std::function<int(ode_system, double, double *)> ode_method;

class ode_system
{
  public:
    int dimension;
    double t_init, t_final;
    double *y0;
    ode_rhs f;
};

class time_stepper
{
  public:
    double dt;
    ode_method F;
    int integrate(ode_system sys, double *yf)
    {
      return F(sys, dt, yf);
    }
};

// Explicit Integrators.
int forward_euler(ode_system, double, double *);
int forward_euler_step(ode_system, double, double *);
int rk4(ode_system, double, double *);
int rk4_step(ode_system, double, double *);

// Implicit Integrators.
int backward_euler(ode_system, double, double *);
int trapezoidal_rule(ode_system, double, double *);

// Parareal Method.
int parareal(ode_system, time_stepper, time_stepper, double *);

#endif
