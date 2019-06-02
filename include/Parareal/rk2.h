#ifndef RK2_H_INCLUDED
#define RK2_H_INCLUDED

#include "core.h"

inline int 
rk2(ode_system&, double, Eigen::VectorXd &);

inline int 
rk2(ode_system&, double, Eigen::MatrixXd &);

#include "rk2.cpp"

#endif
