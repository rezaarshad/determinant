# Determinant Calculation
Calculating the determinant of matrices using MPI (Message Passing Interface) and C.

1. **Using Gaussian Algorithm**  
   Time complexity: O(n³)

2. **Using Laplace Expansion (Recursive)**  
   Time complexity: O(n!).  
   You can find more information about this algorithm here:  
   https://en.wikipedia.org/wiki/Laplace_expansion

# Installing MPICH2

To use the MPI library, you need to install MPICH. You can find more information on how to install MPICH2 here:  
http://mpitutorial.com/tutorials/installing-mpich2

# How to Compile and Run

```bash
mpicc filename.c -lm -o filename  
mpirun -f machinefile -n 4 ./filename  
