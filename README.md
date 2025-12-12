# Parallel Heat Simulation

## Video Demonstration

Watch the simulation in action:

[View full video on SharePoint](https://testusjedu-my.sharepoint.com/:v:/g/personal/rodolph_khoury_net_usj_edu_lb/IQALHFH7zMkuT6lbS4ABE7PkAXb0Y7St7a_HTc6kH2pqmcA?nav=eyJyZWZlcnJhbEluZm8iOnsicmVmZXJyYWxBcHAiOiJPbmVEcml2ZUZvckJ1c2luZXNzIiwicmVmZXJyYWxBcHBQbGF0Zm9ybSI6IldlYiIsInJlZmVycmFsTW9kZSI6InZpZXciLCJyZWZlcnJhbFZpZXciOiJNeUZpbGVzTGlua0NvcHkifX0&e=vAS8Gd)

## 1. Introduction

In this project, we analyze and solve the 2D transient heat-conduction equation using the explicit finite-difference FTCS method.
The formulation follows *Numerical Methods for Engineers* by Chapra & Canale.

The mathematical model comes from Chapter 30.1 — Heat-Conduction Equation, where the authors derive the PDE governing heat diffusion in a plate.

## 2. PDE Classification — Why It Is Parabolic

In PT8.1 Motivation (Table PT8.1), PDEs are classified as:

- Elliptic – steady-state (Laplace equation)
- Parabolic – time-dependent diffusion (heat equation)
- Hyperbolic – wave propagation

The heat-conduction equation is explicitly categorized as a parabolic PDE.

**Why it is parabolic**

A parabolic PDE contains:

- a first derivative in time
- second derivatives in space
- no second derivative in time

The heat equation has exactly this structure.

**Example**

A plate heated along one edge warms gradually over time — a typical parabolic diffusion behavior.

## 3. Finite Difference Discretization (Explicit Method)

The explicit method is described in Chapter 30.2 — Explicit Methods.

### 3.1 Time Approximation (Forward Euler)

Uses a forward difference for the time derivative.

### 3.2 Space Approximation (Central Differences)

Uses central differences for the second spatial derivatives.

### 3.3 Combined Explicit FTCS Update

From Equation (30.2), we derive the update used in the code:

```cpp
unew[i][j] =
    u[i][j] + (DT / (H * H)) *
    (u[i+1][j] + u[i-1][j] +
     u[i][j+1] + u[i][j-1] - 4 * u[i][j]);
```

**Example**

If a cell is 60°C and neighbors are 68, 55, 62, and 50°C, the update computes the new temperature using surrounding diffusion.

## 4. Boundary Conditions

Based on Chapter 29.3 — Boundary Conditions.

We apply Dirichlet conditions:

- Top boundary: 100°C
- Left, right, bottom: 0°C

In code:

```cpp
if (y == LY) T = 100;
if (x == 0 || x == LX || y == 0) T = 0;
```

✔ **Example**

The fixed 100°C top boundary continuously injects heat into the domain.

## 5. Stability Condition

Explicit FTCS is conditionally stable.

From Chapter 30.2.1 — Convergence and Stability, the 2D condition is:

```
Δt ≤ (Δx² * Δy²) / (2 * (Δx² + Δy²))
```

**Why it matters**

Large Δt causes unstable, non-physical temperature explosions.

**Example**

Too large Δt → values jump to thousands of degrees.

## 6. Computational Molecule (Stencil)

Defined in PT8.2.

Our 5-point stencil:

```
    (i, j+1)
(i-1,j) (i,j) (i+1,j)
    (i, j-1)
```

**Example**

The new temperature depends on the four neighbors around each grid point.

## 7. Parallel Implementation (MPI)

Although the book does not cover MPI, the numerical scheme extends naturally:

- The domain is split into horizontal slices
- Neighboring processes exchange halo rows
- Each process runs the same FTCS update
- Boundary conditions are enforced locally

**Example**

If the full domain has 100 rows:

- Process 0 → rows 0–49
- Process 1 → rows 50–99

They exchange one boundary row per timestep.

## 8. Summary of Theory → Code Mapping

| Book Section | Concept Used | Code Location |
|---|---|---|
| PT8.1 | PDE classification | Choosing the heat equation |
| 30.1 | Heat-conduction PDE | Base model |
| 30.2 | Explicit FTCS | heatEquation() |
| 30.2.1 | Stability | Choosing Δt |
| 29.3 | Boundary conditions | 0°C / 100°C edges |
| PT8.2 | Computational molecules | 5-point stencil |

This confirms that the implementation aligns with the textbook.

## 9. Conclusion

This project:

- identifies the PDE as parabolic (PT8.1)
- uses the heat-conduction equation (30.1)
- discretizes using explicit FTCS (30.2)
- respects the stability condition (30.2.1)
- applies Dirichlet boundaries (29.3)
- uses the 5-point stencil (PT8.2)
- extends the model to a parallel MPI simulation

The approach follows the recommended methodology from *Numerical Methods for Engineers*.

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
