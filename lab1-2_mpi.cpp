#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<mpi.h>

void compute_pi(long long int n, long long int* num_in_cycle, long long int* local_n, int comm_sz, long long int* total_num_in_cycle, MPI_Comm comm, int my_rank) {
	*num_in_cycle = 0;
	*local_n = n / comm_sz;
    double x, y, distance_squared;
    srand(time(NULL));
    for (long long int i = 0; i < *local_n; i++) {
		x = (double)rand() / (double)RAND_MAX;
	    x = x * 2 - 1;
	    y = (double)rand() / (double)RAND_MAX;
	    y = y * 2 - 1;
	    distance_squared = x * x + y * y;
	    if (distance_squared <= 1)
			*num_in_cycle = *num_in_cycle + 1;
    }
    MPI_Reduce(num_in_cycle, total_num_in_cycle, 1, MPI_LONG_LONG, MPI_SUM, 0, comm);
    if (my_rank == 0) {
        double pi = (double)*total_num_in_cycle / (double)n * 4;
        printf("the estimate value of pi is %lf\n", pi);
	}
}

int main(int argc, char** argv) {
    long long int n = 1000;
    long long int num_in_cycle, total_num_in_cycle, local_n;
    int my_rank, comm_sz;
    double begin, end, diff;

    MPI_Comm comm;
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &my_rank);
    if (my_rank == 0) begin = MPI_Wtime();
    MPI_Bcast(&n, 1, MPI_LONG_LONG, 0, comm);

    compute_pi(n, &num_in_cycle, &local_n, comm_sz, &total_num_in_cycle, comm, my_rank);
    
    if (my_rank == 0) {
        double end = MPI_Wtime();
        double diff = end - begin;
        printf("time is %9.7f\n", diff);
    }
    MPI_Finalize();
    return 0;
}
