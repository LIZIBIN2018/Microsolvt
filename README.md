# Microsolvt
This project is a mini C++ solver for Poisson equation on some "regular" domain.

## Dependency
Our project will depend on two third-party libraries:
1. jsoncpp
2. eigen


## The code framework
1. Most important class is **Multigrid<int>**, which contains two special tool members and **CycleSolver<int>**;
2. **CycleSolver<int>** is used to indicate the type of cycle, and solve the linear system.\
   **CycleSolver<int>** contains two special tool member **InterpolationOperator<int>** and **RestrictionOperator**, to indicate which type of interpolation and restriction.
  
## The pipeline
1. Generate the grid according to the json file.\
    1.1 Initialize the grid properties.\
    1.2 Generate the tools objects.
2. Solve the grid.\
    2.1 Call the solve function of the grid. (feed f)\
    2.2 Write the result into the grid.
3. Output the result.          
4. Post processing.

## The workflow
1. Build the framwork of multigrid, transfer the main part of work to the solver.
2. Implement
