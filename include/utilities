/**
 * Author: Vincenzo Iannucci
 * Purpose: The library provides an implementation for useful functions
 * like wrapper functions, functions used to dispay results and so on.
 * **/
#ifndef utilities_h
#define utilities_h

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>

#ifdef __cplusplus
extern "C" {
#endif

/*************************************************
 *          Wrapper functions
 ***********************************************/
void *Malloc(size_t size);
void *Calloc(size_t nmemb, size_t size);

/*************************************************
 *          Initialization functions
 ***********************************************/
void initializeRandomVector(double *vector, int N);
void initializeRandomMatrix(double *matrix, int M, int N);
int initMatrixByRowWithValuesFromVector(double *matrix, int M, int N, double *vector, int vector_size);
void initVectorWAnotherVector(double *newVector, double *oldVector, int n);

/*************************************************
 *  Free all the memory dynamically allocated
 ***********************************************/
void freeEverything(void *arg1, ...);

#ifdef __cplusplus
}
#endif
#endif // !utilities_h