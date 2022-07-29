# alvise
While studying various methods on Beam Cross Area analysis getting familiar with work of Alvise Sommariva
https://www.math.unipd.it/~alvise/software.html

Short notes are collected on this page, see more details in each directory

## polcub
uses matlab polyshape to define geometry and subdivides triangles until rquired precision is met.
This is used in [cross-section](https://github.com/simo-11/ancf-rhs/tree/master/cross-section). 
See also [cross-section in ancf main entry](https://github.com/simo-11/ancf-rhs)

## rbfcub
based on https://github.com/alvisesommariva/RBFCUB_2021
Work is focused on [Radial Basis Functions - RBF](https://en.wikipedia.org/wiki/Radial_basis_function) which are not directly usable for Beam Cross Area as
r^2 which is tricky part is not RBF.
