# LowerBounds

Compute lower bounds of eigenvalues of Hermitian operators.
This is a simple code to compute and test lower bounds of hermitian operators,
typically the Hamiltonian operator of some interesting physical system.
It reads the tridiagonal representation of the operator (of whatever dimensions) 
as given by the Lanczos algorithm, please provide the diagonal and first-diagonal elements 
in alpha.dat and beta.dat, respectively. 
The largest dimensional calculation provides reference values to compute the error 
of both the upper and lower bound estimates. Please ignore these information if you know that 
the Lanczos calculation is not yet converged.

Requirements 

Linux/MacOS operating system
gfortran    compiler
LAPACK/BLAS linear package libraries

Compile

"make lower" creates the exec file in your $HOME/bin directory.

Usage

Type "lower" and answer the questions to specify the precise lower bound estimate you'd like to have/test.
The latter include some improved versions of the approximate bounds described in 

Pollak E (2019) "An improved lower bound to the ground-state energy", Journal of Chemical Theory and Computation15(3):1498–1502 

Pollak E (2019) "A tight lower bound to the ground-state energy", Journal of Chemical Theory and Computation15(7):4079–4087

and the rigorous theory described in

Martinazzo R. and Pollak E. (2020) "Lower bounds to eigenvalues of the Schrödinger equation by solution of a ninety year challenge", submitted
