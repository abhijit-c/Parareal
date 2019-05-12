#ifndef INTEGRATORS_H_INCLUDED
#define INTEGRATORS_H_INCLUDED

#include <functional>

/* integrator function declaration syntax.
 * first argument should be ode system, which should be a function object:
 *    double *ode(int dimension, double t, double *y, double *dydt)
 */

typedef std::function<int(double, double *, double *)> ode_rhs;
typedef std::function<int(ode_system, double, double *)> integrator;

class ode_system
{
  public:
    int dimension;
    double t_init, t_final;
    double *y0;
    ode_rhs f;
};

// Explicit Integrators.
int forward_euler(ode_system, double, double *);
int rk4(ode_system, double, double *);

// Implicit Integrators.
int backward_euler(ode_system, double, double *);
int trapezoidal_rule(ode_system, double, double *);

// Parareal Method.
int parareal(ode_system, integrator, integrator, double *);

#endif
