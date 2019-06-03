#ifndef BACKWARD_EULER_H_INCLUDED
#define BACKWARD_EULER_H_INCLUDED

#include "core.h"

inline int 
backward_euler_newton_method(ode_system&, double, Eigen::VectorXd &);

inline int 
backward_euler(ode_system&, double, Eigen::VectorXd &);
inline int 
backward_euler_allt(ode_system&, double, Eigen::MatrixXd &);

#include "backward_euler.cpp"

#endif
