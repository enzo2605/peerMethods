#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <stdarg.h>
#include "peerMethods.h"

void *Malloc(size_t size) {
    void *res = malloc(size);
    if (res == NULL) {
        perror("malloc");
        exit(1);
    }
    return res;
}
       
void *Calloc(size_t nmemb, size_t size) {
    void *res = calloc(nmemb, size);
    if (res == NULL) {
        perror("calloc");
        exit(1);
    }
    return res;
}

void initializeRandomVector(double *vector, int N) {
    int i;
    double randomNumber;

    for (i = 0; i < N; i++) {
        randomNumber = (double)rand() / ((double)RAND_MAX);
        vector[i] = randomNumber;
    }
}

void initializeRandomMatrix(double *matrix, int M, int N) {
    int i, j;
    // generation by columns starting by the index 1
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            matrix[i * N + j] = (double)rand() / ((double)RAND_MAX);
        }
    }
}

int initMatrixByRowWithValuesFromVector(double *matrix, int M, int N, double *vector, int vector_size) {
    if (M * N != vector_size) {
        fprintf(stdout, "\nIt has been impossible to initialize the matrix with the vector passed.");
        fprintf(stdout, "\nThe size of matrix and the vector are incompatible.\nReturned -1\n");
        return -1;
    }

    int k = 0;
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            matrix[j * M + i] = vector[k++];
        }
    }
    return 0;
}

void initVectorWAnotherVector(double *newVector, double *oldVector, int n) {
    int i;
    for (i = 0; i < n; i++) {
        newVector[i] = oldVector[i];
    }
}

void freeEverything(void *arg1, ...) {
    va_list args; // list of arguments
    void *vp; // pointer to the i-th arguments
    // Free the first argument
    free(arg1);
    va_start(args, arg1);
    // Until you reach the end of the arguments, free the argument
    // pointed by vp
    while ((vp = va_arg(args, void *)) != 0) {
        free(vp);
    }
    va_end(args);
}