#ifndef RHS_LINEARD1D_H
#define RHS_LINEARD1D_H

#include <Eigen/Core>

#define LAMBDA 1

// u' = (LAMBDA)*u
inline int 
linear1d(double t, Eigen::VectorXd &yf, Eigen::VectorXd &dydt)
{
  dydt(0) = LAMBDA*yf(0);
  return 0;
}

#endif
