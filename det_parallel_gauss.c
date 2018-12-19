#include <math.h>
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define MASTER 0            /* rank of first task */
#define FROM_MASTER 1            /* setting a message type */
#define FROM_WORKER 2

double Partition(double **a, int s, int end, int n);

int main(argc, argv)
        int argc;
        char *argv[];
{
    int NumberProcesses,    /* number of tasks in partition */
            rank,        /* a task identifier */
            NumberWorkers,    /* number of worker tasks */
            source,        /* task id of message source */
            destination,        /* task id of message destinationination */
            messagetype,        /* message type */
            extra, offset,
            i, j, k, rc, len;    /* misc */
    double det, StartTime, EndTime, read_StartTime, read_EndTime, print_StartTime, print_EndTime, **matrix, *buffer, determinant_of_matrix;
    int n;    /*number of rows and columns in matrix */

    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Status status;
    rc = MPI_Init(&argc, &argv);
    rc |= MPI_Comm_size(MPI_COMM_WORLD, &NumberProcesses);
    rc |= MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(hostname, &len);
    StartTime = MPI_Wtime();
    NumberWorkers = NumberProcesses - 1;


    if (!rank) {
        read_StartTime = MPI_Wtime();
        printf("Enter the n for n*n matrix : ");
        fflush(stdout);
        scanf("%d", &n);
        read_EndTime = MPI_Wtime();

        /**************** Create random number for matrix **************************/
        buffer = (double *) malloc(sizeof(double) * n * n);
        printf("Number of  tasks = %d\n", NumberProcesses);
        print_StartTime = MPI_Wtime();
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                buffer[i * n + j] = (double) (rand() % 20) + 1;;
                //printf("%3.2lf    ",buffer[i*n+j]);
            }
            //printf("\n");
        }
        print_EndTime = MPI_Wtime();

        /* send matrix data to the worker tasks */
        float temp;
        temp = n / NumberProcesses;
        temp += 0.5;
        offset = temp;

        extra = n % NumberProcesses;
        messagetype = FROM_MASTER;
        //MPI_Bcast (&offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // MPI_Bcast (buffer, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (destination = 1; destination <= NumberWorkers; destination++) {

            //printf("   sending %d task to task %d\n",offset,destination);
            MPI_Send(&n, 1, MPI_INT, destination, messagetype, MPI_COMM_WORLD);
            //MPI_Send(&offset, 1, MPI_INT, destination, messagetype, MPI_COMM_WORLD);
            MPI_Send(buffer, n * n, MPI_DOUBLE, destination, messagetype, MPI_COMM_WORLD);

        }
        matrix = (double **) malloc((n) * sizeof(double[n]));
        for (k = 0; k < n; ++k)
            matrix[k] = (double *) malloc((n) * sizeof(double));

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                matrix[i][j] = buffer[i * n + j];
            }


        determinant_of_matrix = Partition(matrix, 0, offset, n);
        printf("%s calculate it's part with determinant=%3.2lf\n", hostname, determinant_of_matrix);

        free(buffer);

        /* wait for results from all worker tasks and computer determinant of matrix */
        messagetype = FROM_WORKER;

        for (i = 1; i <= NumberWorkers; i++) {
            source = i;
            MPI_Recv(&det, 1, MPI_DOUBLE, source, messagetype, MPI_COMM_WORLD, &status);
            determinant_of_matrix += det;
        }
        //end time
        EndTime = MPI_Wtime();
        printf("Elapsed time is %f\n",
               ((EndTime - StartTime) - (print_EndTime - print_StartTime) - (read_EndTime - read_StartTime)));
        printf("Determinant of matrix is :%3.2lf\n", determinant_of_matrix);

    }

/**************************** worker task ************************************/
    if (rank) {
        messagetype = FROM_MASTER;
        MPI_Recv(&n, 1, MPI_INT, MASTER, messagetype, MPI_COMM_WORLD, &status);
        buffer = (double *) malloc(sizeof(double) * n * n);
        //MPI_Recv(&offset, 1, MPI_INT, MASTER, messagetype, MPI_COMM_WORLD, &status);
        MPI_Recv(buffer, n * n, MPI_DOUBLE, MASTER, messagetype, MPI_COMM_WORLD, &status);

        float temp;
        temp = n / NumberProcesses;
        temp += 0.5;
        offset = temp;

        int end;
        int start = (rank) * offset;
        if ((rank) == NumberWorkers)
            end = n;
        else
            end = (start + offset);

        matrix = (double **) malloc((n) * sizeof(double[n]));
        for (k = 0; k < n; ++k)
            matrix[k] = (double *) malloc((n) * sizeof(double));

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                matrix[i][j] = buffer[i * n + j];
            }


        det = Partition(matrix, start, end, n);
        int h = 0;
        for (h = 0; h < n; ++h)
            free(matrix[h]);
        free(matrix);
        free(buffer);
        messagetype = FROM_WORKER;
        /* 	Send answer to master node. */
        MPI_Send(&det, 1, MPI_DOUBLE, MASTER, messagetype, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}


double MatrixDeterminant(int nDim, double *pfMatr) {
    double fDet = 1.;
    double fMaxElem;
    double fAcc;
    int i, j, k, m;

    for (k = 0; k < (nDim - 1); k++) // base row of matrix
    {
        // search of line with max element
        fMaxElem = fabs(pfMatr[k * nDim + k]);
        m = k;
        for (i = k + 1; i < nDim; i++) {
            if (fMaxElem < fabs(pfMatr[i * nDim + k])) {
                fMaxElem = pfMatr[i * nDim + k];
                m = i;
            }
        }

        // permutation of base line (index k) and max element line(index m)
        if (m != k) {
            for (i = k; i < nDim; i++) {
                fAcc = pfMatr[k * nDim + i];
                pfMatr[k * nDim + i] = pfMatr[m * nDim + i];
                pfMatr[m * nDim + i] = fAcc;
            }
            fDet *= (-1.);
        }

        if (pfMatr[k * nDim + k] == 0.) return 0.0;

        // trianglulation of matrix
        for (j = (k + 1); j < nDim; j++) // current row of matrix
        {
            fAcc = -pfMatr[j * nDim + k] / pfMatr[k * nDim + k];
            for (i = k; i < nDim; i++) {
                pfMatr[j * nDim + i] = pfMatr[j * nDim + i] + fAcc * pfMatr[k * nDim + i];
            }
        }
    }

    for (i = 0; i < nDim; i++) {
        fDet *= pfMatr[i * nDim + i]; // diagonal elements multiplication
    }

    return fDet;
}

double Partition(double **a, int s, int end, int n) {
    int i, j, j1, j2;
    double det = 0;
    double **m = NULL;

    det = 0;                      // initialize determinant of sub-matrix

    // for each column in sub-matrix
    for (j1 = s; j1 < end; j1++) {
        // get space for the pointer list
        m = (double **) malloc((n - 1) * sizeof(double *));

        for (i = 0; i < n - 1; i++)
            m[i] = (double *) malloc((n - 1) * sizeof(double));

        for (i = 1; i < n; i++) {
            j2 = 0;

            for (j = 0; j < n; j++) {
                if (j == j1) continue;

                m[i - 1][j2] = a[i][j];


                j2++;
            }
        }
        int dim = n - 1;
        double fMatr[dim * dim];
        for (i = 0; i < dim; i++) {
            for (j = 0; j < dim; j++) {
                fMatr[i * dim + j] = m[i][j];
                // printf("%3.2lf    ",fMatr[i*nDim+j]);
            }
            //printf("\n");
        }

        det += pow(-1.0, 1.0 + j1 + 1.0) * a[0][j1] * MatrixDeterminant(dim, fMatr);

        //free(fMatr);
        for (i = 0; i < n - 1; i++) free(m[i]);

        free(m);

    }

    return (det);
}
