#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <tchar.h>
#include <math.h>
#include <time.h>

#define N 100000

int _tmain(int argc, _TCHAR* argv[]) {
	omp_set_num_threads(2);
	int i, j, num = 0;
	double begin, end, diff;
	begin = omp_get_wtime();

#pragma omp parallel for reduction(+:num) private(j)
	for (i=2; i<=N; i++) {
		int s = sqrt(i * 1.0);
		for (j=2; j<=s; j++) {
			if (i % j == 0) 
				break;
		}
		if (j > s) num++;
	}
	end = omp_get_wtime();
	printf("%d内素数共有%d个\n", N, num);
	diff = end - begin;
	printf("并行时间是%9.7f\n", diff);
	return 0;
}
