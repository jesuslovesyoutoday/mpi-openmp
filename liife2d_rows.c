/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, 2016
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <mpi.h>

#define ind(i, j) (((i + l->nx) % l->nx) + ((j + l->ny) % l->ny) * (l->nx))

typedef struct {
	int nx, ny;
	int *u0;
	int *u1;
	int steps;
	int save_steps;
} life_t;

void life_init(const char *path, life_t *l);
void life_free(life_t *l);
void row_step(life_t *l, int num);
void life_step(life_t *l, int rank, int numprocess);
void life_save_vtk(const char *path, life_t *l);
void gathering(life_t *l, int rank, int numprocess);

int main(int argc, char **argv)
{
	if (argc != 2) {
		printf("Usage: %s input file.\n", argv[0]);
		return 0;
	}
	life_t l;
	life_init(argv[1], &l);
	
	int i;
	char buf[100];

    MPI_Init(&argc, &argv);
    int numprocess, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double start = MPI_Wtime();
	for (i = 0; i < l.steps; i++) {
		if ((i % l.save_steps == 0) && rank == numprocess-1) {
			sprintf(buf, "resss/life_%06d.vtk", i);
			printf("Saving step %d to '%s'.\n", i, buf);
			life_save_vtk(buf, &l);
		}
        printf("%d, before step\n", rank);
		life_step(&l, rank, numprocess);
        printf("%d, after step\n", rank);
        gathering(&l, rank, numprocess);
	}
    double time = MPI_Wtime() - start;

    FILE* fptr;
    fptr = fopen("time.txt", "a");
    fprintf(fptr, "%d ", numprocess);  
    fprintf(fptr, "%f\n", time);  

    MPI_Finalize();		

	life_free(&l);
	return 0;
}

/**
 * Загрузить входную конфигурацию.
 * Формат файла, число шагов, как часто сохранять, размер поля, затем идут координаты заполненых клеток:
 * steps
 * save_steps
 * nx ny
 * i1 j2
 * i2 j2
 */
void life_init(const char *path, life_t *l)
{
	FILE *fd = fopen(path, "r");
	assert(fd);
	assert(fscanf(fd, "%d\n", &l->steps));
	assert(fscanf(fd, "%d\n", &l->save_steps));
	printf("Steps %d, save every %d step.\n", l->steps, l->save_steps);
	assert(fscanf(fd, "%d %d\n", &l->nx, &l->ny));
	printf("Field size: %dx%d\n", l->nx, l->ny);

	l->u0 = (int*)calloc(l->nx * l->ny, sizeof(int));
	l->u1 = (int*)calloc(l->nx * l->ny, sizeof(int));
	
	int i, j, r, cnt;
	cnt = 0;
	while ((r = fscanf(fd, "%d %d\n", &i, &j)) != EOF) {
		l->u0[ind(i, j)] = 1;
		cnt++;
	}
	printf("Loaded %d life cells.\n", cnt);
	fclose(fd);
}

void life_free(life_t *l)
{
	free(l->u0);
	free(l->u1);
	l->nx = l->ny = 0;
}

void life_save_vtk(const char *path, life_t *l)
{
	FILE *f;
	int i1, i2, j;
	f = fopen(path, "w");
	assert(f);
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk2d\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %d %d 1\n", l->nx+1, l->ny+1);
	fprintf(f, "SPACING %d %d 0.0\n", 1, 1);
	fprintf(f, "ORIGIN %d %d 0.0\n", 0, 0);
	fprintf(f, "CELL_DATA %d\n", l->nx * l->ny);
	
	fprintf(f, "SCALARS life int 1\n");
	fprintf(f, "LOOKUP_TABLE life_table\n");
	for (i2 = 0; i2 < l->ny; i2++) {
		for (i1 = 0; i1 < l->nx; i1++) {
			fprintf(f, "%d\n", l->u0[ind(i1, i2)]);
		}
	}
	fclose(f);
}

void gathering(life_t *l, int rank, int numprocess)
{
    int i = 0;
    int count = l->ny / numprocess;

    if (rank < numprocess - 1)
    {
        int* sdata = l->u0 + rank * count * l->nx;
        MPI_Send(sdata, l->nx * count, MPI_INT, numprocess-1, 0, MPI_COMM_WORLD);
    }
    else if (rank == numprocess - 1)
    {
        for (i = 0; i < numprocess - 1; i ++)
        {
            int* rdata = l->u0 + i * count * l->nx;
            MPI_Recv(rdata, l->nx * count, MPI_INT, i, 0, MPI_COMM_WORLD,
                                                      MPI_STATUS_IGNORE);
        }
    }
}

void row_step(life_t *l, int num)
{
    int j = num;
    int i;
    for (i = 0; i < l->nx; i++) {
	    int n = 0;
	    n += l->u0[ind(i+1, j)];
	    n += l->u0[ind(i+1, j+1)];
	    n += l->u0[ind(i,   j+1)];
	    n += l->u0[ind(i-1, j)];
	    n += l->u0[ind(i-1, j-1)];
	    n += l->u0[ind(i,   j-1)];
	    n += l->u0[ind(i-1, j+1)];
	    n += l->u0[ind(i+1, j-1)];
	    l->u1[ind(i,j)] = 0;
	    if (n == 3 && l->u0[ind(i,j)] == 0) {
		    l->u1[ind(i,j)] = 1;
	    }
	    if ((n == 3 || n == 2) && l->u0[ind(i,j)] == 1) {
		    l->u1[ind(i,j)] = 1;
	    }
    }
}

void life_step(life_t *l, int rank, int numprocess)
{
    int i, j, N;
    int num = 0;

    int count = l->ny / numprocess;
    int count_last = count + (l->ny % numprocess); 

    if (rank == numprocess - 1)
    {
        N = count_last;
    }
    else
    {
        N = count;
    }

    while(num < N)
    {
        row_step(l, num + rank * count);
        num = num + 1;
    }
    

    
    if (rank != numprocess - 1)
    {
        int* sdata2 = l->u1 + (rank * count + count - 1) * l->nx;
        MPI_Send(sdata2, l->nx, MPI_INT, rank+1, 0, MPI_COMM_WORLD);

        int* rdata2 = l->u1 + (rank * count + count) * l->nx;
        MPI_Recv(rdata2, l->nx, MPI_INT, rank+1, 0, MPI_COMM_WORLD,
                                                  MPI_STATUS_IGNORE);
    }
    if (rank != 0)
    {
        int* sdata1 = l->u1 + rank * count * l->nx;
        MPI_Send(sdata1, l->nx, MPI_INT, rank-1, 0, MPI_COMM_WORLD);

        int* rdata1 = l->u1 + rank * (count - 1) * l->nx;
        MPI_Recv(rdata1, l->nx, MPI_INT, rank-1, 0, MPI_COMM_WORLD,
                                                  MPI_STATUS_IGNORE);
    }
    if (rank == 0)
    {
        int* sdata = l->u1;
        MPI_Send(sdata, l->nx, MPI_INT, numprocess-1, 0, MPI_COMM_WORLD);

        int* rdata = l->u1 + (l->ny - 1) * l->nx;
        MPI_Recv(rdata, l->nx, MPI_INT, numprocess-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank == numprocess-1)
    {
        int* sdata = l->u1 + (l->ny - 1) * l->nx;
        MPI_Send(sdata, l->nx, MPI_INT, 0, 0, MPI_COMM_WORLD);

        int* rdata = l->u1;
        MPI_Recv(rdata, l->nx, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

	int *tmp;
	tmp = l->u0;
	l->u0 = l->u1;
	l->u1 = tmp;

}


