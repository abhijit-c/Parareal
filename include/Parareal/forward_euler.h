#ifndef FORWARD_EULER_H_INCLUDED
#define FORWARD_EULER_H_INCLUDED

#include "core.h"

inline int 
forward_euler(ode_system&, double, Eigen::VectorXd &);
inline int 
forward_euler_allt(ode_system&, double, Eigen::MatrixXd &);


#include "forward_euler.cpp"

#endif
