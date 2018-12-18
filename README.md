# determinant
Calculating determinant of matrices using MPI (Message Passing Interface) and C.
1. Calculating determinant of matrices using G with O(n^3)
2. Calculating determinant of matrices using Laplace expansion (recursive) with O(n!)!
   https://en.wikipedia.org/wiki/Laplace_expansion

# Installing mpich2
For using mpi library you need to install mpich. you can find more information for installing mpich2 at here :
   http://mpitutorial.com/tutorials/installing-mpich2
   
#Compiling and running using mpi
for compiling mpi based c file you need to run :
    mpicc filename.c –lm –o filename

for running files:
    mpicc –f machinefile –n nm ./filename
    nm is number of machins which you want to run in parallel.
   

