#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

int isPrime(int num) {
	int flag = 1;
	int s = sqrt(num * 1.0);
	for (int j=2; j<=s; j++) {
		if (num%j == 0) {
			flag = 0;
			break;
		}
	}
	return flag;
}

void main(int argc, char* argv[]) {
	int n = 100000;
	int myid, numprocs, i, cnt, sum, mycnt;
	double begin, end, diff;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (myid == 0) begin = MPI_Wtime();
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	sum = 0;
	for (i = myid*2+1; i<=n; i+=numprocs*2)
		sum += isPrime(i);
	mycnt = sum;
	MPI_Reduce(&mycnt, &cnt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if (myid == 0) {
		printf("range(1,%d) has %d prime number.\n", n, cnt);
		end = MPI_Wtime();
		diff = end - begin;
		printf("%d process time is %9.7f\n", myid, diff);
	}

	MPI_Finalize();
}
