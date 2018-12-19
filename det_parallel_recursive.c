
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

double Determinant(double **a, int start, int end, int n);

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
    double det, StartTime, EndTime, read_StartTime, read_EndTime, print_StartTime, print_EndTime, **matrix, *buffer_send, determinant_of_matrix;
    int n;    /*number of rows and columns in matrix */

    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Status status;
    rc = MPI_Init(&argc, &argv);
    rc |= MPI_Comm_size(MPI_COMM_WORLD, &NumberProcesses);
    rc |= MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(hostname, &len);
    StartTime = MPI_Wtime();
    if (rc != 0)
        printf("error initializing MPI and obtaining task ID information\n");


    NumberWorkers = NumberProcesses - 1;
    //if(NumberWorkers>n)
    //NumberWorkers=n;

/**************************** Cordinator task ************************************/
    if (!rank) {
        printf("%s is ready with  task ID = %d\n", hostname, rank);
        read_StartTime = MPI_Wtime();

        if (argv[1] != NULL) {
            FILE *fp;
            fp = fopen(argv[1], "r");
            double item;
            fscanf(fp, "%d", &n);
            buffer_send = (double *) malloc(sizeof(double) * n * n);
            int i, j;
            for (i = 0; i < n * n; i++) {

                fscanf(fp, "%lf", &buffer_send[i]);
            }
            printf("\n");
        } else {
            printf("Enter the n for n*n matrix : ");
            fflush(stdout);
            scanf("%d", &n);
            buffer_send = (double *) malloc(sizeof(double) * n * n);
            for (i = 0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    buffer_send[i * n + j] = (double) (rand() % 20) + 1;
                    //printf("%3.2lf    ",buffer_send[i*n+j]);
                }
                //printf("\n");
            }

        }

        read_EndTime = MPI_Wtime();

        /**************** Create random number for matrix **************************/
        //NumberWorkers+master=numberworker+1
        printf("Number of tasks = %d\n", NumberWorkers + 1);
        //print_StartTime=MPI_Wtime();         
        //print_EndTime=MPI_Wtime();
        //temp is used for rounding offset

        float temp;
        temp = n / NumberProcesses;
        temp += 0.5;
        offset = temp;

        extra = n % NumberProcesses;
        messagetype = FROM_MASTER;
        //MPI_Bcast (&offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //MPI_Bcast (buffer_send, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        for (destination = 1; destination <= NumberWorkers; destination++) {

            MPI_Send(&n, 1, MPI_INT, destination, messagetype, MPI_COMM_WORLD);
            //MPI_Send(&offset, 1, MPI_INT, destination, messagetype, MPI_COMM_WORLD);
            MPI_Send(buffer_send, n * n, MPI_DOUBLE, destination, messagetype, MPI_COMM_WORLD);

        }
        matrix = (double **) malloc((n) * sizeof(double[n]));
        for (k = 0; k < n; ++k)
            matrix[k] = (double *) malloc((n) * sizeof(double));

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++) {
                matrix[i][j] = buffer_send[i * n + j];
            }


        determinant_of_matrix = Determinant(matrix, 0, offset, n);
        printf("%s calculate it's part with determinant=%3.2lf\n", hostname, determinant_of_matrix);

        for (i = 0; i < n; ++i)
            free(matrix[i]);
        free(matrix);
        free(buffer_send);

        /* wait for results from all worker tasks and computer determinant of matrix */
        messagetype = FROM_WORKER;

        for (i = 1; i <= NumberWorkers; i++) {
            source = i;
            MPI_Recv(&det, 1, MPI_DOUBLE, source, messagetype, MPI_COMM_WORLD, &status);
            determinant_of_matrix += det;
        }
        //end time
        EndTime = MPI_Wtime();
        printf("Elapsed time is %f\n", ((EndTime - StartTime) - (read_EndTime - read_StartTime)));
        printf("Determinant of matrix is :%3.2lf\n", determinant_of_matrix);

    }

/**************************** worker task ************************************/
    if (rank) {
        messagetype = FROM_MASTER;
        MPI_Recv(&n, 1, MPI_INT, MASTER, messagetype, MPI_COMM_WORLD, &status);
        buffer_send = (double *) malloc(sizeof(double) * n * n);
        //MPI_Recv(&offset, 1, MPI_INT, MASTER, messagetype, MPI_COMM_WORLD, &status);
        MPI_Recv(buffer_send, n * n, MPI_DOUBLE, MASTER, messagetype, MPI_COMM_WORLD, &status);
        printf("%s is ready with  task ID = %d\n", hostname, rank);
        //if(NumberWorkers>n)
        //NumberWorkers=n;

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
                matrix[i][j] = buffer_send[i * n + j];
            }


        det = Determinant(matrix, start, end, n);
        printf("%s calculate it's part with determinant=%3.2lf\n", hostname, det);


        int h = 0;
        for (h = 0; h < n; ++h)
            free(matrix[h]);
        free(matrix);
        free(buffer_send);
        messagetype = FROM_WORKER;
        /* 	Send answer to master node.      */
        MPI_Send(&det, 1, MPI_DOUBLE, MASTER, messagetype, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}


double Determinant(double **a, int start, int end, int n) {
    int i, j, j1, j2;
    double det = 0;
    double **m = NULL;


    if (n < 1) {}

    else if (n == 1) {
        det = a[0][0];
    } else if (n == 2) {

        det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    } else {

        det = 0;

        for (j1 = start; j1 < end; j1++) {

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

            det += pow(-1.0, 1.0 + j1 + 1.0) * a[0][j1] * Determinant(m, 0, n - 1, n - 1);

            for (i = 0; i < n - 1; i++) free(m[i]);

            free(m);

        }
    }
    return (det);
}



