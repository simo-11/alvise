# rbfcub
Beam Cross Area based on https://github.com/alvisesommariva/RBFCUB_2021

# README from source

Cubature by RBF over polygonal regions RBFCUB: Matlab package for near-optimal meshless cubature on general polygons (2021)

Matlab codes for computing moments and cubature rules on polygonal regions (via triangulation or alternatively by Gauss-Green theorem in polar coordinates). The regions can be general, i.e. polygons that are not simply connected or even disconnected. It requires Matlab built-in environment polyshape. With respect to a previous version, depending on the values of the integrand at the centers, we optimize the shape parameter so to achieve better cubature results. The authors of the following work are R. Cavoretto, A. De Rossi, A. Sommariva and M. Vianello.

Object: RBF cubature on polygonal regions (with almost optimal shape parameter).

The main routines are: RBF_cub_polygon_OPT.m that determines the RBF based cubature rule on the polygonal region; RBF_moms_polar.m and RBF_moms_tri.m that determines the RBF moments on the polygonal region; RBF_optimal_scale.m that chooses the optimal RBF scale parameter;

The software includes several demos: demo_RBF_cub_OPT.m that performs a battery of tests, allowing different integrands, polygonal regions, RBF, number of centers and choices of the RBF shape parameter. demo_RBF_cub_OPT_basic that performs a simple example, on a certain integrand, polygonal region, RBF, number of centers and with almost optimal RBF shape parameter (useful to adapt the code to user experiments).

For the demos of the previous version published in RBF moment computation and meshless cubature on general polygonal regions, see: RBF_cub_polygon.m that determines the RBF based cubature rule on the polygonal region (not optimized shape parameter);

Source:

Papers:
    R. Cavoretto, A. De Rossi, A. Sommariva and M. Vianello, RBFCUB: a numerical package for near-optimal meshless cubature on general polygons;
    A. Sommariva and M. Vianello, RBF moment computation and meshless cubature on general polygonal regions 

# Additions

## Radial Basis Functions (RBF)
https://en.wikipedia.org/wiki/Radial_basis_function
