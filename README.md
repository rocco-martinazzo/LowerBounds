# LowerBounds

Compute lower bounds of eigenvalues of Hermitian operators
This is a simple code to compute and test lower bounds of hermitian operators,
typically the Hamiltonian operator of some interesting physical system.
It reads the tridiagonal representation of the operator (of whatever dimensions) 
as given by the Lanczos algorithm, please provide the diagonal and first-diagonal elements 
in alpha.dat and beta.dat, respectively. 
The largest dimensional calculation provides reference values to compute the error 
of both the upper and lower bound estimates. Please ignore these information if you know that 
the Lanczos calculation is not yet converged.

Requirements 

Lapack/Blas linear package libraries


Compile

"make lower" creates the exec file in your $HOME/bin directory.
