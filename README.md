# alvise
While studying various methods on Beam Cross Area analysis getting familiar with work of Alvise Sommariva
https://www.math.unipd.it/~alvise/software.html

Short notes are collected on this page, see more details in each directory

## polcub
uses matlab polyshape to define geometry and subdivides triangles until rquired precision is met.
[cub_polygon_adaptive.m](polcub/cub_polygon_adaptive.m) is modified slightly and used in [cross-section](https://github.com/simo-11/ancf-rhs/tree/master/cross-section). 
See also [cross-section in ancf main entry](https://github.com/simo-11/ancf-rhs)

## splcub
Software computes algebraic cubature rules with fixed degree of precision over spline-curvilinear polygons. It includes an in-domain routine incurvpolygon.m that determines if a point is inside one of such regions.
The cubature formula has all points inside the domain as well as positive weights.
[paper](https://www.math.unipd.it/~marcov/pdf/splinecatch.pdf)

## polygauss
Software computes algebraic cubature rules with fixed degree of precision over general polygons (convex or not convex, connected or disconnected, simply or not simply connected). It needs Matlab-built in polyshape environment.
[paper](https://www.math.unipd.it/~marcov/pdf/polygons.pdf)


## rbfcub
based on https://github.com/alvisesommariva/RBFCUB_2021
Work is focused on [Radial Basis Functions - RBF](https://en.wikipedia.org/wiki/Radial_basis_function) which are not directly usable for Beam Cross Area as
r^2 which is tricky part is not RBF.
