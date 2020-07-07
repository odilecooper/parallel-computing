#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <time.h>

int N = 5000000;
int* array;
//int array[16] = {16,1,15,2,14,13,3,4,12,11,10,5,6,7,9,8 };

int cmp(const void* a, const void* b) {
    if (*(int*)a < *(int*)b) return -1;
    if (*(int*)a > * (int*)b) return 1;
    else return 0;
}

void printArray() {
    int i;
    for (i = 0; i < N; i++) {
        printf("%4d ", array[i]);
        if (i%16 == 15) printf("\n");
    }
    printf("\n");
}

void PSRS(int* array, int N) {
    int p, myId;
    int subArraySize, startIndex, endIndex;
    double begin, end, diff;

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

    if (myId == 0) {
//        printf("original array:\n");
//        printArray();
        begin = MPI_Wtime();
    }

    // 获取起始位置和子数组大小
    startIndex = myId * N / p;
    endIndex = (myId + 1) * N / p;
    subArraySize = endIndex - startIndex;

    MPI_Barrier(MPI_COMM_WORLD);

    // 局部排序
    qsort(array + startIndex, subArraySize, sizeof(array[0]), cmp);
//    printArray();

    // 选取样本
    int i;
    int* sample = (int*)malloc(p * sizeof(int));
    for (i = 0; i < p; i++) {
        sample[i] = array[startIndex + (i * (N / (p * p)))];
    }

    if (p == 1) {
        end = MPI_Wtime();
        diff = end - begin;
        printf("time:%9.7lf\n", diff);
//        printArray();
        return;
    }

    int* gatherSamp = (int*)malloc(p * p * sizeof(sample[0]));
    int* privots = (int*)malloc((p - 1) * sizeof(sample[0])); //主元

    // 收集样本
    MPI_Gather(sample, p, MPI_INT, gatherSamp, p, MPI_INT, 0, MPI_COMM_WORLD);
    if (myId == 0) {
        // 样本排序
        qsort(gatherSamp, p * p, sizeof(sample[0]), cmp);
        // 选取主元
        for (i = 0; i < p - 1; i++) {
            privots[i] = gatherSamp[(i + 1) * p];
//            printf("%d ", privots[i]);
        }
//        printf("\n");
    }

    // 主元划分
    MPI_Bcast(privots, p - 1, MPI_INT, 0, MPI_COMM_WORLD);
    int* partitionSizes = (int*)malloc(p * sizeof(int));
    for (i = 0; i < p; i++) {
        partitionSizes[i] = 0;
    }
    int index = 0;
    for (i = 0; i < subArraySize; i++) {
        if (array[startIndex + i] > privots[index]) {
            //当前位置数字大于主元，则进行下一个划分
            index ++;
        }
        if (index == p) {
            //最后一次划分，子数组总长减掉当前位置即可得到最后一个子数组划分的大小
            partitionSizes[p - 1] = subArraySize - i + 1;
            break;
        }
        partitionSizes[index]++; //划分大小自增
    }

    free(gatherSamp);
    free(privots);

    // 全局交换
    int totalSize = 0;
    int* sendDisp = (int*)malloc(p * sizeof(int));
    int* recvDisp = (int*)malloc(p * sizeof(int));
    int* newPartitionSizes = (int*)malloc(p * sizeof(int));

    MPI_Alltoall(partitionSizes, 1, MPI_INT, newPartitionSizes, 1, MPI_INT, MPI_COMM_WORLD);

    // 计算划分的总大小，并给新划分分配空间
    for (i = 0; i < p; i++) {
        totalSize += newPartitionSizes[i];
    }
    int *newPartitions = (int*)malloc(totalSize * sizeof(int));

    sendDisp[0] = 0;
    recvDisp[0] = 0;
    // 计算相对sendbuf的位移，此位移处存放输出到进程的数据
    // 计算相对recvbuf的位移，此位移处存放从进程接收的数据
    for (i = 1; i < p; i++) {
        sendDisp[i] = partitionSizes[i - 1] + sendDisp[i - 1];
        recvDisp[i] = newPartitionSizes[i - 1] + recvDisp[i - 1];
    }

    //发送数据，实现n次点对点通信
    MPI_Alltoallv(&(array[startIndex]), partitionSizes, sendDisp, MPI_INT, newPartitions, newPartitionSizes, recvDisp, MPI_INT, MPI_COMM_WORLD);

    free(sendDisp);
    free(recvDisp);

    int* sortedSubList;
    int* indexes, * partitionEnds, * subListSizes, totalListSize;

    indexes = (int*)malloc(p * sizeof(int));
    partitionEnds = (int*)malloc(p * sizeof(int));
    indexes[0] = 0;
    totalListSize = newPartitionSizes[0];
    for (i = 1; i < p; i++) {
        totalListSize += newPartitionSizes[i];
        indexes[i] = indexes[i - 1] + newPartitionSizes[i - 1];
        partitionEnds[i - 1] = indexes[i];
    }
    partitionEnds[p - 1] = totalListSize;

    sortedSubList = (int*)malloc(totalListSize * sizeof(int));
    subListSizes = (int*)malloc(p * sizeof(int));
    recvDisp = (int*)malloc(p * sizeof(int));

    // 归并排序
    int j;
    for (i = 0; i < totalListSize; i++) {
        int lowest = INT_MAX;
        int ind = -1;
        for (j = 0; j < p; j++) {
            if ((indexes[j] < partitionEnds[j]) && (newPartitions[indexes[j]] < lowest)) {
                lowest = newPartitions[indexes[j]];
                ind = j;
            }
        }
        sortedSubList[i] = lowest;
        indexes[ind] += 1;
    }

    // 发送子列表大小到根进程
    MPI_Gather(&totalListSize, 1, MPI_INT, subListSizes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // 计算根进程上相对recvbuf的偏移量
    if (myId == 0) {
        recvDisp[0] = 0;
        for (i = 1; i < p; i++) {
            recvDisp[i] = subListSizes[i - 1] + recvDisp[i - 1];
        }
    }

    //发送排好序的子列表到根进程
    MPI_Gatherv(sortedSubList, totalListSize, MPI_INT, array, subListSizes, recvDisp, MPI_INT, 0, MPI_COMM_WORLD);

    free(partitionEnds);
    free(sortedSubList);
    free(indexes);
    free(subListSizes);
    free(recvDisp);

    if (myId == 0) {
        end = MPI_Wtime();
        diff = end - begin;
        printf("time:%9.7lf\n", diff);
//        printf("sorted array:\n");
//        printArray();
    }

    if (p > 1) {
        free(newPartitions);
    }
    free(partitionSizes);
    free(newPartitionSizes);
    free(sample);
}

int main(int argc, char* argv[]) {
    array = (int*)malloc(N * sizeof(int));
    int i;
    for (i = 0; i < N; i++) {
        array[i] = (int)rand() % (4*N);
    }

    MPI_Init(&argc, &argv);
    PSRS(array, N);

    free(array);
    MPI_Finalize();
    return 0;
}
