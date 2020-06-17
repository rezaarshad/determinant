# determinant
Calculating determinant of matrices using MPI (Message Passing Interface) and C.
1. using Gaussian algorithm with time complexity of O(n^3)
2. using Laplace expansion (recursive) with time complexity of O(n!)!
  You can find more information about this algorithm at here: https://en.wikipedia.org/wiki/Laplace_expansion

# Installing mpich2
For using mpi library you need to install mpich. you can find more information for installing mpich2 at here :
   http://mpitutorial.com/tutorials/installing-mpich2
   
# How to compile and run
mpicc filename.c –lm –o filename  
mpicc –f machinefile –n 4 ./filename  

Notes:
*-n 4 means four-process parallel application, running four copies of the executable named filename.  
*machinefile is the file which indicate which hosts to start the processes on. so the file contains contributed machine    names.
   

