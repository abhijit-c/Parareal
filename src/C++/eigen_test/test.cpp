#include <iostream>
#include <Eigen/Dense>
#include <omp.h>

void wow(Eigen::VectorXd &m)
{
      m(1) = 22;
}

int main()
{
  Eigen::VectorXd v(5);
  Eigen::VectorXd b(5);
  b = v;
  std::cout << v << std::endl;
  wow(v);
  std::cout << v << std::endl;
  std::cout << b << std::endl;
  return 0;
}
