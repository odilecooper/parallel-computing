#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#define N 64 // 小球数
#define size 8 // size = sqrt(N)
#define L 0.07 // 正方形区域大小=(size-1)/100
// #define G 6.67E-7
// #define M 10000
#define GM 6.67E-3
#define cycle 200 // 周期数
#define delta_t 0.001 // 每周期时间

typedef struct ball {
	double px, py; // 位置
	double vx, vy; // 速度
	double ax, ay; // 加速度
}ball;

ball ball_list[N];

void compute_force(int index) {
	ball_list[index].ax = 0;
	ball_list[index].ay = 0;
	for (int i = 0; i < N; i++)	{
		if (i != index) {
			double dx = ball_list[i].px - ball_list[index].px;
			double dy = ball_list[i].py - ball_list[index].py;
			double d = (dx * dx + dy * dy); //d^2
			if (d == 0) continue;
			d *= sqrt(d);// d^3
			// printf("%lf %lf %lf\n", dx, dy, d);
			// a.x = Gm * dx / d^3
			ball_list[index].ax += GM * dx / d;
			ball_list[index].ay += GM * dy / d;
			// printf("%lf %lf \n", GM * dx / d, GM * dy / d);
		}
	}
	// printf("%d a: %lf %lf\n",index,ball_list[index].ax,ball_list[index].ay);
}

void compute_velocities(int index) {
	ball_list[index].vx += ball_list[index].ax * delta_t;
	ball_list[index].vy += ball_list[index].ay * delta_t;
	// printf("%d v: %lf %lf\n",index,ball_list[index].vx,ball_list[index].vy);
}

void compute_positions(int index) {
	ball_list[index].px += ball_list[index].vx * delta_t;
	if (ball_list[index].px > L)
		ball_list[index].px = L;
	else if (ball_list[index].px < 0)
		ball_list[index].px = 0;

	ball_list[index].py += ball_list[index].vy * delta_t;
	if (ball_list[index].py > L)
		ball_list[index].py = L;
	else if (ball_list[index].py < 0)
		ball_list[index].py = 0;

	// printf("%d p: %lf %lf\n",index,ball_list[index].px,ball_list[index].py);
}

void print() {
	int table[size][size]={0};
	int i, j;
	for (i = 0; i < N; i++) {
		table[(int)(ball_list[i].px*100)][(int)(ball_list[i].py*100)]++;
	}
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			// printf("(%lf, %lf) ", i*size+j, ball_list[i * size + j].px, ball_list[i * size + j].py);
			printf("%3d", table[i][j]);
		}
		printf("\n");
	}
	printf("end of printing\n\n");
}

int main(int argc, char* argv[]) {
	int i, j;
	for (i = 0; i < N; i++) {
		ball_list[i].px = 0.01 * (i % size);
		ball_list[i].py = 0.01 * (i / size);
		ball_list[i].vx = 0;
		ball_list[i].vy = 0;
		ball_list[i].ax = 0;
		ball_list[i].ay = 0;
		// printf("%d (%lf,%lf)\n", i, ball_list[i].px, ball_list[i].py);
	}

	int myid, numprocs;
	double begin, end, diff;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	double* mpi_buffer = (double*)malloc(sizeof(double) * 7*N*N*cycle);
	MPI_Buffer_attach(mpi_buffer, sizeof(double) * 7*N*N * cycle);

	//模拟开始
	if(myid==0) begin = MPI_Wtime();
	
	for (i = 0; i < cycle; i++) {
		for (j = 0; j < numprocs; j++) {
			if (j != myid)
				MPI_Bsend((ball_list + (N / numprocs) * myid), sizeof(ball) * N / numprocs, MPI_BYTE, j, i * 10 + myid, MPI_COMM_WORLD);
		}

		for (j = 0; j < numprocs; j++) {
			if (j != myid) {
				MPI_Status status;
				MPI_Recv((ball_list + (N / numprocs) * j), sizeof(ball) * N / numprocs, MPI_BYTE, j, i * 10 + j, MPI_COMM_WORLD, &status);
			}
		}

		for (j = (N / numprocs) * myid; j < (N / numprocs) * (myid + 1); j++) {
			compute_force(j);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		for (j = (N / numprocs) * myid; j < (N / numprocs) * (myid + 1); j++) {
			compute_velocities(j);
			compute_positions(j);
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}

		// 模拟结束，输出
		if (myid == 0) {
			end = MPI_Wtime();
			diff = end - begin;
			printf("time:%9.7lf\n", diff);
		}

		if (myid != 0) {
			//printf("sid %d %d %d\n",myid,ball_list+(N/numprocs)*myid,sizeof(ball)*N/numprocs);
			MPI_Send((ball_list + (N / numprocs) * myid), sizeof(ball) * N / numprocs, MPI_BYTE, 0, myid, MPI_COMM_WORLD);
		}

		if (myid == 0) {
			for (int i = 1; i < numprocs; i++) {
				MPI_Status status;
				//printf("rid %d %d %d\n",myid,ball_list+(N/numprocs)*myid,sizeof(ball)*N/numprocs);
				MPI_Recv((ball_list + (N / numprocs) * i), sizeof(ball) * N / numprocs, MPI_BYTE, i, i, MPI_COMM_WORLD, &status);
			}
			print();
		}

	MPI_Finalize();
	return 0;
}