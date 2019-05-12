#ifndef INTEGRATORS_H_INCLUDED
#define INTEGRATORS_H_INCLUDED

#include <functional>

/* integrator function declaration syntax.
 * first argument should be ode system, which should be a function object:
 *    double *ode(int dimension, double t, double *y, double *dydt)
 */

typedef std::function<int(double, double *, double *)> ode_rhs;

class ode_system
{
  public:
    int dimension;
    double t_init, t_final;
    double *y0;
    ode_rhs f;
};

int forward_euler(ode_system, double, double *);

#endif
