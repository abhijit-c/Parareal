#ifndef RHS_LORENZ_H
#define RHS_LORENZ_H

#include <Eigen/Core>

#define SIGMA 10
#define RHO 28
#define BETA 2.6667

/* Lorenz Attractor:
 * dx/dt = SIGMA*(y-x)
 * dy/dt = x*(RHO-z)-y
 * dz/dt = x*y-BETA*z
 */
inline int 
lorenz(double t, Eigen::VectorXd &yf, Eigen::VectorXd &dydt)
{
  dydt(0) = SIGMA*(yf(1) - yf(0));
  dydt(1) = yf(0)*(RHO - yf(2)) - yf(1);
  dydt(2) = yf(0)*yf(1) - BETA*yf(2);
  return 0;
}

// Jacobian of Lorenz attractor
inline int 
lorenz_jac(double t, Eigen::VectorXd &yf, Eigen::MatrixXd &J)
{
  J(0,0) = -SIGMA; J(1,0) = SIGMA; J(2,0) = 0    ;
  J(0,1) = RHO   ; J(1,1) = -1   ; J(2,1) = 0    ;
  J(0,2) = 0     ; J(1,2) = 0    ; J(2,2) = -BETA;
  return 0;
}

#endif
