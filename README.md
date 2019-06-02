# Parareal

## What is this

The final project for both my Numerical Methods II and High Performance
Computing graduate classes. In general, the objective was to explore the
parareal algorithm: 

- Implementing it and write an ODE ecosystem around it.
- Efficiency goals:
    - Implement in OpenMP
    - Optimize with pipelining.
    - Observe how choices with course and fine propagators affect speedup.
    - Understand and derive theoretical efficiency limits.
    - Examine weak and strong scaling of the algorithm.
- Numerical Analysis Goals:
    - Understand how Parareal's stability changes under propagator choices.
    - Derive convergence results under specific cases.

## How to use

If you would like to compile and run this code, there's a few requirements.
First, you must have a gcc compiler compliant with the C++11 standard, and
OpenMP installed. In addition, you must have access to the make unix tool. 

After this, Parareal is a header-only library and therefore all you need to do
is include the relevant .h files after adding the `include` directory to your
include path. For example, a sample compilation call would be:

```
g++ -std=c++11 -I /path/to/Parareal/include sample.c
```

In your program, you should always include Eigen and the core file. The core
file contains the typedefs and classes needed to interact with the library, and
the library itself uses Eigen's datastructures internally.

```c
#include <Eigen/Core>
#include <Parareal/core.h>
```

There is currently no way to statically link to this library.
