# Parallel Heat Simulation

## Video Demonstration

[View the simulation in action](https://testusjedu-my.sharepoint.com/:v:/g/personal/rodolph_khoury_net_usj_edu_lb/IQDT1FTt_VKHRrFFpnb6u5FdAQHGeqv7cVF4cx-vQFr2D0A?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=xcbZa1)

## Overview

This project solves the 2D transient heat-conduction equation using the explicit finite-difference FTCS (Forward Time Centered Space) method. The implementation is based on *Numerical Methods for Engineers* by Chapra & Canale (Chapter 30).

![Heat Conduction PDE](images/heat_conduction_pde.png)

## The Heat Equation

The heat-conduction equation is a parabolic PDE—meaning it has a first derivative in time and second derivatives in space. This characterizes diffusion problems where heat gradually spreads through a domain over time. (PT8.1 — PDE Classification)

![PDE Classification](images/pde_classification.png)

## Numerical Method

### Time Discretization

We use forward differences (Forward Euler) for the time derivative (Ch. 30.2):

![Forward Euler Method](images/forward_euler.png)

### Space Discretization

Central differences approximate the spatial derivatives (Ch. 30.2):

![Central Differences](images/central_differences.png)

### The Update Rule

Combining time and space discretization gives us the FTCS update (Eq. 30.2):

![FTCS Update Equation](images/ftcs_update_equation.png)

The stability coefficient is:

![Stability Coefficient](images/stability_coefficient.png)

Which leads to the discrete update:

```cpp
unew[i][j] =
    u[i][j] + (DT / (H * H)) *
    (u[i+1][j] + u[i-1][j] +
     u[i][j+1] + u[i][j-1] - 4 * u[i][j]);
```

## Boundary Conditions

We apply Dirichlet conditions (Ch. 29.3):
- Top edge: 100°C (heat source)
- All other edges: 0°C

```cpp
if (y == LY) T = 100;
if (x == 0 || x == LX || y == 0) T = 0;
```

## Stability

For stability in 2D, the timestep must satisfy (Ch. 30.2.1):

```
Δt ≤ (Δx² * Δy²) / (2 * (Δx² + Δy²))
```

Violating this causes temperature spikes.

## Computational Stencil

The 5-point stencil for the 2D heat equation (PT8.2):

```
    (i, j+1)
(i-1,j) (i,j) (i+1,j)
    (i, j-1)
```

## Parallelization

The domain is decomposed into horizontal slices, with each process handling a portion of rows. Processes exchange boundary rows at each timestep to update ghost cells.

## Implementation Details

| Concept | Source | Implementation |
|---|---|---|
| PDE | Ch. 30.1 | Base model |
| FTCS Method | Ch. 30.2 | `heatEquation()` |
| Stability | Ch. 30.2.1 | Timestep selection |
| Boundary Conditions | Ch. 29.3 | Edge handling |
| Stencil | PT8.2 | 5-point pattern |

## Prerequisites

- C++ compiler with C++17 support
- CMake ≥ 3.16

## Build Instructions

1. Open a terminal and navigate to the project root:
```bash
cd parallel_computing_project
```

Create and enter the build directory:

```bash
mkdir -p build
cd build
```

Configure the project using CMake:

```bash
cmake ..
```

Build the project:

```bash
cmake --build .
```

After building, the executable `heat_solver.exe` will be located in the build directory.

## Run Instructions

From the build directory, run the program:

```bash
export PATH="/c/Program Files/Microsoft MPI/Bin:$PATH"
mpiexec -n x ./heat_solver.exe   # x is the number of processes
```

On Linux:

```bash
mpiexec --oversubscribe -n 4 ./heat_solver
```
