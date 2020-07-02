#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<omp.h>

#define N 1000

int main() {
    int num_in_cycle = 0;
    srand(time(NULL));
    double x, y, distance_point;
    int i;
    double begin, end, diff;
    begin = omp_get_wtime();

#pragma omp parallel for reduction(+:num_in_cycle) private(i,x,y,distance_point)
    for (i = 0; i < N; i++) {
        x = (double)rand() / (double)RAND_MAX;
        y = (double)rand() / (double)RAND_MAX;
        distance_point = x * x + y * y;
        if (distance_point <= 1) {
             num_in_cycle++;
		}
	}
    double estimate_pi = (double)num_in_cycle / N * 4;
    printf("the estimate value of pi is %lf\n", estimate_pi);
    end = omp_get_wtime();
    diff = end - begin;
    printf("time: %9.7f\n", diff);
    return 0;
}