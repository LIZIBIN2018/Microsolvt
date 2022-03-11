# Microsolvt
This project is a mini C++ solver of harmonic equation on some "regular" domain.

## Dependency
Our project will depend on two third-party libraries:
1. jsoncpp
2. eigen

## What we need to do
1. Define a json file to specify the input.
2. Define the grid about the domain.
3. Find out the linear equations about the discrete grid of the domain.
4. Solve the linear equations and store the solution.
5. Transform the solution to a proper format.

## Some details about our workflow

### json file
1. Specify the step length of grid, the center and radius of circular disk, the boundary condition type and the boundary value in `.json` file as the following properties:
`grid_h`,`circ_c`,`circ_r`,`bdry_type`.`bdry_val`.
2. Maybe we can generate the json file by another program.....

### Grid
1. If we need a **adaptive grid**, maybe we need not predefine the step length....

### Linear equations
nothing

### Output format
**To be determined**
