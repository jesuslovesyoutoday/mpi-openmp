#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

double LatencyOneNode(int rank)
{
    int sdata;
	int rdata;
    int N = 100000000;

	double start = MPI_Wtime();
    int i;
    for (i = 0; i < N; i++)
    {
        if (rank == 0)
        {
            MPI_Send(&sdata, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
        }
        if (rank == 1)
        {
            MPI_Recv(&rdata, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
                                               MPI_STATUS_IGNORE);
        }
        if (rank == 1)
        {
            MPI_Send(&sdata, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 0)
        {
            MPI_Recv(&rdata, 1, MPI_INT, 1, 0, MPI_COMM_WORLD,
                                               MPI_STATUS_IGNORE);
        }
    }
    double time = MPI_Wtime() - start;

    return time;
}

double LatencyTwoNodes(int rank)
{
    int sdata;
    int rdata;
    int N = 100000;

    double start = MPI_Wtime();
    int i;
    for (i = 0; i < N; i++)
    {
        if (rank == 0)
        {
            MPI_Send(&sdata, 1, MPI_INT, 4, 0, MPI_COMM_WORLD);
        }
        if (rank == 4)
        {
            MPI_Recv(&rdata, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
                                               MPI_STATUS_IGNORE);
        }
        if (rank == 4)
        {
            MPI_Send(&sdata, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 0)
        {
            MPI_Recv(&rdata, 1, MPI_INT, 4, 0, MPI_COMM_WORLD,
                                               MPI_STATUS_IGNORE);
        }
    }
    double time = MPI_Wtime() - start;

    return time;
}

double* ExperimentalOneNode(int rank)
{
    //int count = 100000;
    int S = 100000;
    int count = 1000;
    double* times = (double*)calloc(S, sizeof(double));
    
    int size;
    for (size = 1; size < S; size ++)
    {
        int* sdata = (int*)calloc(size, sizeof(int));
        int* rdata = (int*)calloc(size, sizeof(int));

        double start = MPI_Wtime();
        int i;
        for (i = 0; i < count; i++)
        {
            if (rank == 0)
            {
                MPI_Send(sdata, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
            }
            if (rank == 1)
            {
                MPI_Recv(rdata, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
                                                   MPI_STATUS_IGNORE);
            }
            if (rank == 1)
            {
                MPI_Send(sdata, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
            if (rank == 0)
            {
                MPI_Recv(rdata, 1, MPI_INT, 1, 0, MPI_COMM_WORLD,
                                                   MPI_STATUS_IGNORE);
            }
        }
        double time = MPI_Wtime() - start;
        times[size] = time;

        free(sdata);
        free(rdata); 
    }
    return times;
}

double* ExperimentalTwoNodes(int rank)
{
    int S = 100000;
    //int count = 100000;
    int count = 100;
    double* times = (double*)calloc(S, sizeof(double));
    
    int size;
    for (size = 1; size < S; size = size + 3)
    {
        int* sdata = (int*)calloc(size, sizeof(int));
        int* rdata = (int*)calloc(size, sizeof(int));

        double start = MPI_Wtime();
        int i;
        for (i = 0; i < count; i++)
        {
            if (rank == 0)
            {
                MPI_Send(sdata, 1, MPI_INT, 4, 0, MPI_COMM_WORLD);
            }
            if (rank == 4)
            {
                MPI_Recv(rdata, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
                                                   MPI_STATUS_IGNORE);
            }
            if (rank == 4)
            {
                MPI_Send(sdata, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
            if (rank == 0)
            {
                MPI_Recv(rdata, 1, MPI_INT, 4, 0, MPI_COMM_WORLD,
                                                   MPI_STATUS_IGNORE);
            }
        }
        double time = MPI_Wtime() - start;

        times[size] = time;

        free(sdata);
        free(rdata); 
    }
    return times;
}

double TheorCoeffOneNode(int rank)
{
    int N = 100000000;
    int I = 100;
    int* sdata = (int*)calloc(N, sizeof(int));
    int* rdata = (int*)calloc(N, sizeof(int));

    double start = MPI_Wtime();
    int i;
    for (i = 0; i < I; i++)
    {
        if (rank == 0)
        {
            MPI_Send(sdata, N, MPI_INT, 1, 0, MPI_COMM_WORLD);
        }
        if (rank == 1)
        {
            MPI_Recv(rdata, N, MPI_INT, 0, 0, MPI_COMM_WORLD,
                                                MPI_STATUS_IGNORE);
        }
        if (rank == 1)
        {
            MPI_Send(sdata, N, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 0)
        {
            MPI_Recv(rdata, N, MPI_INT, 1, 0, MPI_COMM_WORLD,
                                                MPI_STATUS_IGNORE);
        }
    }
    double time = MPI_Wtime() - start;

    free(sdata);
    free(rdata);

    return time;
}

double TheorCoeffTwoNodes(int rank)
{
    int N = 100000000;
    int I = 100;
    int* sdata = (int*)calloc(N, sizeof(int));
    int* rdata = (int*)calloc(N, sizeof(int));

    double start = MPI_Wtime();
    int i;
    for (i = 0; i < I; i++)
    {
        if (rank == 0)
        {
            MPI_Send(sdata, N, MPI_INT, 4, 0, MPI_COMM_WORLD);
        }
        if (rank == 4)
        {
            MPI_Recv(rdata, N, MPI_INT, 0, 0, MPI_COMM_WORLD,
                                                MPI_STATUS_IGNORE);
        }
        if (rank == 4)
        {
            MPI_Send(sdata, N, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        if (rank == 0)
        {
            MPI_Recv(rdata, N, MPI_INT, 4, 0, MPI_COMM_WORLD,
                                                MPI_STATUS_IGNORE);
        }
    }
    double time = MPI_Wtime() - start;

    free(sdata);
    free(rdata);

    return time;
}

int main(int argc, char** argv[])
{
    MPI_Init(&argc, argv);

    int numprocess, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    FILE* fptr;
    
    if (rank == 0)
    {
        fptr = fopen("latency_one_0.txt", "w");
    }
    if (rank == 1)
    {   
        fptr = fopen("latency_one_1.txt", "w");
    } 
    if (rank == 0 || rank == 1)
    {
        double time = LatencyOneNode(rank);
        fprintf(fptr, "%f", time);
        fclose(fptr);
    }

    if (rank == 0)
    {
        fptr = fopen("latency_two_0.txt", "w");
    }
    if (rank == 4)
    {   
        fptr = fopen("latency_two_4.txt", "w");
    } 
    if (rank == 0 || rank == 4)
    {
        double time = LatencyTwoNodes(rank);
        fprintf(fptr, "%f", time);
        fclose(fptr);
    }

    if (rank == 0)
    {
        fptr = fopen("teorcoef_one_0.txt", "w");
    }
    if (rank == 1)
    {   
        fptr = fopen("teorcoef_one_1.txt", "w");
    } 
    if (rank == 0 || rank == 1)
    {
        double time = TheorCoeffOneNode(rank);
        fprintf(fptr, "%f", time);
        fclose(fptr);
    }

    if (rank == 0)
    {
        fptr = fopen("teorcoef_two_0.txt", "w");
    }
    if (rank == 4)
    {   
        fptr = fopen("teorcoef_two_4.txt", "w");
    } 
    if (rank == 0 || rank == 4)
    {
        double time = TheorCoeffTwoNodes(rank);
        fprintf(fptr, "%f", time);
        fclose(fptr);
    }
    if (rank == 0)
    {
        fptr = fopen("exp_one_0.txt", "w");
    }
    if (rank == 1)
    {   
        fptr = fopen("exp_one_1.txt", "w");
    } 
    if (rank == 0 || rank == 1)
    {
	double* times = ExperimentalOneNode(rank);
        int i;
        for (i = 0; i < 99999; i++)
        {
            fprintf(fptr, "%f ", times[i]);
        }
        fclose(fptr);
        free(times);
    }

    if (rank == 0)
    {
        fptr = fopen("exp_two_0.txt", "a");
    }
    if (rank == 4)
    {   
        fptr = fopen("exp_two_4.txt", "a");
    } 
    if (rank == 4)
    {
        double* times = ExperimentalTwoNodes(rank);
        int i;
        for (i = 0; i < 99999; i++)
        {
            fprintf(fptr, "%f\n", times[i]);
        }
        fclose(fptr);
        free(times);
    }
    if (rank == 0)
    {
	    double* times = ExperimentalTwoNodes(rank);
        int i;
        for (i = 0; i < 99999; i++)
        {
            fprintf(fptr, "%f\n", times[i]);
        }
        fclose(fptr);
        free(times);
    }
    
    
    MPI_Finalize();	

    return 0;
}
