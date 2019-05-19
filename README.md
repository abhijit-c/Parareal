# Parareal

## What is this

The final project for both my Numerical Methods II and High Performance
Computing graduate classes. In general, the objective was to explore the
parareal algorithm: 

- Implementing it and write an ODE ecosystem around it.
- Efficiency goals:
    - Implement in OpenMP
    - Optimize with pipelining.
    - Observe how choices with course and fine propagators affect runtime.
    - Understand and derive theoretical efficiency limits.
- Numerical Analysis Goals:
    - Understand how Parareal's stability changes under propagator choices.
    - Derive convergence results under specific cases.

## How to use

If you would like to compile and run this code, there's a few requirements.
First, you must have a gcc compiler compliant with the C++11 standard, and
OpenMP installed. In addition, you must have access to the make and ar unix
tools. 

First, you must compile the parareal library, to do so cd into
`src/ode_integrators/` and then call `make`. This will construct the object file
`lib_parareal.a`. Now anytime you want to use the library, include the
`integrators.h` file, and link `lib_parareal.a`.

## Todo

- Confirm convergence analysis with plots.
- Write implicit euler solver, and compare results numerically.
- Test Scalability (strong and weak).
- Write up parallel efficiency analysis.
- Detail Parareal algorithm with pseudocode.
- Write up quick slides for presentation.
