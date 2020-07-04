#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define N 1000000
#define num_cycle 300
#define pos_cnt 10000
#define v0 0
#define v_max 10
#define p 5

typedef struct car
{
	int v;
	int d;
	int pos;
}car;
car car_list[N];

int v_count[v_max+1];
int pos_count[v_max*(num_cycle+N)/pos_cnt + 1]; //count per pos_cnt

int main(int argc, char* argv[]) {
	//初始化
	int i, j;
	for (i = 0; i < N; i++) {
		car_list[i].v = v0;
		car_list[i].d = v_max;
		car_list[i].pos = (v_max+1) * i; //车尾位置
	}
	for (i = 0; i <= v_max; i++)
		v_count[i] = 0;

	int myid, numprocs;
	double begin, end, diff;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	car* mpi_buffer = (car *)malloc(sizeof(car) * N);
	MPI_Buffer_attach(mpi_buffer, sizeof(car) * N);

	if (myid == numprocs-1) begin = MPI_Wtime();
	for (j = 0; j < num_cycle; j++) {
		i = (N / numprocs * myid); //每组第一辆车
		if (myid != 0) {
			MPI_Bsend(&(car_list[i].pos), 1, MPI_INT, myid - 1, j * 10 + myid, MPI_COMM_WORLD);
		}

		for (; i< N/numprocs*(myid + 1)-1; i++) {
			// 更新距离和速度
			car_list[i].d = car_list[i+1].pos - car_list[i].pos - 1;
			if (car_list[i].d > car_list[i].v && car_list[i].v < v_max)
				car_list[i].v++;
			if (car_list[i].d <= car_list[i].v)
				car_list[i].v = car_list[i].d - 1;
			srand(i * N + j);
			if (car_list[i].v > 0) {
				int r = rand() % 10;
				if (r < p) {
					car_list[i].v--;
				}
			}
			// 更新位置
			car_list[i].pos += car_list[i].v;
		}
		// 每组最后一辆车
		if (myid != numprocs - 1) {
			int temp;
			MPI_Status status;
			MPI_Recv(&(temp), 1, MPI_INT, myid+1, j*10+myid+1, MPI_COMM_WORLD, &status);
			//更新距离
			car_list[i].d = temp - car_list[i].pos;
		}
		//更新速度
		if (car_list[i].d > car_list[i].v && car_list[i].v < v_max)
			car_list[i].v++;
		if (car_list[i].d <= car_list[i].v)
			car_list[i].v = car_list[i].d - 1;
		srand((unsigned)time(NULL));
		if (car_list[i].v > 0 && rand() % 10 < p)
			car_list[i].v--;
		//更新位置
		car_list[i].pos += car_list[i].v;

		//MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for (i = 0; i < numprocs - 1; i++) {
		if (myid == i) {
			MPI_Send((car_list + i*N/numprocs), sizeof(car)*N/numprocs, MPI_BYTE, numprocs-1, myid, MPI_COMM_WORLD);
		}
	}
	if (myid == numprocs - 1) {
		end = MPI_Wtime();
		diff = end - begin;
		printf("time is %9.7f\n", diff);

		MPI_Status status;
		for (i = 0; i < numprocs - 1; i++){
			MPI_Recv((car_list + i*N/numprocs), sizeof(car)*N/numprocs, MPI_BYTE, i, i, MPI_COMM_WORLD, &status);
		}
/*
		for (j = 0; j < N; j++) {
			printf("car %5d: v=%2d, pos=%7d, d=%2d\n", j, car_list[j].v, car_list[j].pos, car_list[j].d);
		}

		for (i = 0; i < N; i++) {
			v_count[car_list[i].v]++;
			pos_count[car_list[i].pos / pos_cnt]++;
		}

		for (i = 0; i <= v_max; i++) {
			printf("v=%2d, car_num=%6d\n", i, v_count[i]);
		}

		for (i = 0; i < v_max * (num_cycle + N) / pos_cnt + 1; i++) {
			printf("pos(%7d,%7d), car_num=%6d\n", i*pos_cnt, (i+1)*pos_cnt, pos_count[i]);
		}
*/
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
