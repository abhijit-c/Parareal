#ifndef RK4_H_INCLUDED
#define RK4_H_INCLUDED

#include "core.h"

inline int 
rk4(ode_system&, double, Eigen::VectorXd &);

inline int 
rk4(ode_system&, double, Eigen::MatrixXd &);

#include "rk4.cpp"

#endif
