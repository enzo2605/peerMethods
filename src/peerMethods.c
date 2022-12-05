#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "cblas.h"
#include "peerMethods.h"

void initReturnStruct(return_values *rv) {
    rv->t = NULL;
    rv->y = NULL;
    rv->yT = NULL;
}

int saveResultsInFile(const char* fileName, return_values result) {
    // Open the file
    FILE *filePtr;
    filePtr = fopen(fileName, "w+");
    // Check possible errors
    if (filePtr == NULL) {
        perror("output file error");
        return 1;
    }
    int numArrays = 3;
    fprintf(filePtr, "%d\n", numArrays);
    // Write y_T
    fprintf(filePtr, "%d\n", result.yT_size);
    for (int i = 0; i < result.yT_size; i++) {
        fprintf(filePtr, "%.15f\n", result.yT[i]);
    }
    // Write y
    fprintf(filePtr, "%d\n", result.y_rows * result.y_cols);
    for (int i = 0; i < result.y_rows; i++) {
        for (int j = 0; j < result.y_cols; j++) {
            fprintf(filePtr, "%.15f\n", result.y[j * result.y_rows + i]);
        }
    }
    // Write t
    // Write y_T
    fprintf(filePtr, "%d\n", result.t_size);
    for (int i = 0; i < result.t_size; i++) {
        fprintf(filePtr, "%.15f\n", result.t[i]);
    }
    // Close the file
    fclose(filePtr);
    return 0;
}

void defineLMatrix(double **L, int *LSize, double Delta_x) {
    // Calculate Ldiff
    int sizeTempDiagOne, sizeTempDiagMinusOne;
    // Multiply identity matrix by a scalar -2.0f
    double *eyeM = eyeD(M);
    cblas_dscal(M * M, -2.0f, eyeM, 1);

    double *onesVector = onesD(M - 1);
    double *tempDiagOne = diagD(onesVector, M - 1, 1, &sizeTempDiagOne);
    double *tempDiagMinusOne = diagD(onesVector, M - 1, -1, &sizeTempDiagMinusOne);
    // Sum eyeM with tempDiagOne
    cblas_daxpy(M * M, 1.0f, eyeM, 1, tempDiagOne, 1);
    // Sum the resulting matrix with tempDiagMinusOne
    cblas_daxpy(M * M, 1.0f, tempDiagOne, 1, tempDiagMinusOne, 1);

    // Copy data to Ldiff matrix
    double *Ldiff = (double *)Calloc(M * M, sizeof(double));
    memcpy(Ldiff, tempDiagMinusOne, M * M * sizeof(double));
    // Multiply Ldiff by scalar alpha
    double alpha = 1.0f / (Delta_x * Delta_x);
    cblas_dscal(M * M, alpha, Ldiff, 1);
    // Using periodic boundary conditions, the matrix LDiff must actually
    // be slightly modified, setting LDiff (M,1) = LDiff (1,M) = 1 / Delta_x^2
    Ldiff[M - 1] = alpha;
    Ldiff[(M - 1) * M] = alpha;

    // Define DLdiff
    double *DLdiff = (double *)Calloc(M * M, sizeof(double));
    memcpy(DLdiff, Ldiff, M * M * sizeof(double));
    cblas_dscal(M * M, D, DLdiff, 1);

    // Define dLdiff
    double *dLdiff = (double *)Calloc(M * M, sizeof(double));
    memcpy(dLdiff, Ldiff, M * M * sizeof(double));
    cblas_dscal(M * M, d, dLdiff, 1);

    // Define matrix L
    *L = threeBlockDiagD(M, Ldiff, DLdiff, dLdiff, LSize);

    // Free all the memory dynamically allocated
    freeEverything(eyeM, onesVector, tempDiagOne, tempDiagMinusOne, Ldiff, DLdiff, dLdiff, (void *)0);
}

double *Sherratt(const double *y0, int y0Size, const double *L, int Lsize, int *sherrattSize) {
    // Declaring three dynamic array representing the three function contained in y0
    double *U = (double *)Calloc(M, sizeof(double));
    double *V = (double *)Calloc(M, sizeof(double));
    double *W = (double *)Calloc(M, sizeof(double));
    for (int i = 0; i < M; i++) {
        U[i] = y0[i];
        V[i] = y0[i + M];
        W[i] = y0[i + 2 * M];
    }

    int newSize = y0Size;
    double *fun = (double *)Calloc(newSize, sizeof(double));
    for (int i = 0; i < M; i++) {
        fun[i] = W[i] * U[i] * (U[i] + H * V[i]) - B1 * U[i] - S * U[i] * V[i];
        fun[i + M] = F * W[i] * V[i] * (U[i] + H * V[i]) - B2 * V[i];
        fun[i + 2 * M] = a - W[i] - W[i] * (U[i] + V[i]) * (U[i] + H * V[i]);
    }

    //printf("H: %f\nB1: %f\nS: %f\nF: %f\nB2: %f\na: %f\n", H, B1, S, F, B2, a);
    //saveVectorsInFile("../parallel/good.txt", 1, fun, y0Size, (void *)0);
    //exit(0);

    // Matrix by vector product
    double *Ly = (double *)Calloc(newSize, sizeof(double));
    cblas_dgemv(CblasRowMajor, CblasNoTrans, Lsize, Lsize, 1, L, Lsize, y0, 1, 0, Ly, 1);

    // Punctual sum of two vectors
    cblas_daxpy(Lsize, 1.0f, fun, 1, Ly, 1);
    *sherrattSize = Lsize;

    // Free all the useless memory allocated
    freeEverything(U, V, W, fun, (void *)0);

    return Ly;
}

double *RungeKutta4th(double h, double t0, const double *y0, int y0Size, const double *L, int Lsize, int *ySize) {
    double c[4] = { 0.0f, 1.0f / 2.0f, 1.0f / 2.0f, 1.0f };
    double A[4][4] = { 
                        { 0.0f, 0.0f, 0.0f, 0.0f }, 
                        { 1.0f / 2.0f, 0.0f, 0.0f, 0.0f }, 
                        { 0.0f, 1.0f / 2.0f, 0.0f, 0.0f }, 
                        { 0.0f, 0.0f, 1.0f, 0.0f }
                    };
    double b[4] = { 1.0f / 6.0f, 1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 6.0f };

    // Compute Y1
    const double *Y1 = y0;
    //printDVector(Y1, y0Size, "Y1");

    // Compute Y2
    int fY1_size;
    double *fY1 = Sherratt(Y1, y0Size, L, Lsize, &fY1_size);
    double *Y2 = (double *)Calloc(fY1_size, sizeof(double));
    memcpy(Y2, fY1, fY1_size * sizeof(double));
    cblas_dscal(fY1_size, h * A[1][0], Y2, 1);
    cblas_daxpy(fY1_size, 1.0f, y0, 1, Y2, 1);
    //printDVector(Y2, y0Size, "Y2");

    // Compute Y3
    int fY2_size;
    double *fY2 = Sherratt(Y2, y0Size, L, Lsize, &fY2_size);
    double *Y3 = (double *)Calloc(fY2_size, sizeof(double));
    memcpy(Y3, fY2, fY2_size * sizeof(double));
    cblas_dscal(fY2_size, h * A[2][1], Y3, 1);
    cblas_daxpy(fY2_size, 1.0f, y0, 1, Y3, 1);
    //printDVector(Y3, y0Size, "Y3");

    // Compute Y4
    int fY3_size;
    double *fY3 = Sherratt(Y3, y0Size, L, Lsize, &fY3_size);
    double *Y4 = (double *)Calloc(fY3_size, sizeof(double));
    memcpy(Y4, fY3, fY3_size * sizeof(double));
    cblas_dscal(fY3_size, h * A[3][2], Y4, 1);
    cblas_daxpy(fY3_size, 1.0f, y0, 1, Y4, 1);
    //printDVector(Y4, y0Size, "Y4");

    // Compute y
    int fY4_size;
    double *fY4;
    fY4 = Sherratt(Y4, y0Size, L, Lsize, &fY4_size);
    //saveVectorsInFile("../parallel/good.txt", 4, fY1, y0Size, fY2, y0Size, fY3, y0Size, fY4, y0Size, (void *)0);
    //exit(0);


    cblas_dscal(fY1_size, h * b[0], fY1, 1);
    cblas_dscal(fY2_size, h * b[1], fY2, 1);
    cblas_dscal(fY3_size, h * b[2], fY3, 1);
    cblas_dscal(fY4_size, h * b[3], fY4, 1);

    double *y = (double *)Calloc(y0Size, sizeof(double));
    cblas_daxpy(y0Size, 1.0f, y0, 1, y, 1);
    cblas_daxpy(y0Size, 1.0f, fY1, 1, y, 1);
    cblas_daxpy(y0Size, 1.0f, fY2, 1, y, 1);
    cblas_daxpy(y0Size, 1.0f, fY3, 1, y, 1);
    cblas_daxpy(y0Size, 1.0f, fY4, 1, y, 1);
    *ySize = y0Size;
    //printDVector(y, *ySize, "y");

    return y;
}

void fPeerClassic_twoStages(int N, double *t_span, int t_span_size, const double *L, int Lsize, const double *y0, int y0_size, return_values *collect_result) {
    /******************************* 
     * Fixing method coefficients 
     * ****************************/
    double b11 = 0.0f, b21 = 0.0f, b12 = 1.0f, b22 = 1.0f, c1 = 0.0f, c2 = 1.0f, r21 = 0.0f;

    double a11 = -((b11 - 2 * b11 * c1 - pow(c1, 2) + b11 * pow(c1, 2)) / (2 * (-1 + c1)));
    double a12 = -((b11 + 2 * c1 - 2 * b11 * c1 - pow(c1, 2) + b11 * pow(c1, 2)) / (2 * (-1 + c1)));
    double a21 = -((-1 + b21 - 2 * b21 * c1 + b21 * pow(c1, 2) + 2 * c1 * r21) / (2 * (-1 + c1)));
    double a22 = -((3 + b21 - 2 * c1 - 2 * b21 * c1 + b21 * pow(c1, 2) - 2 * r21) / (2 * (-1 + c1)));
    //fprintf(stdout, "a11: %f\na12: %f\na21: %f\na22: %f\n", a11, a12, a21, a22);

    double c[STAGES] = { c1, c2 };
    //printDVector(c, STAGES, "c");

    double *A = (double *)Calloc(STAGES * STAGES, sizeof(double));
    double tempA[STAGES * STAGES] = { a11, a12, a21, a22 };
    initMatrixByRowWithValuesFromVector(A, STAGES, STAGES, tempA, STAGES * STAGES);
    //printDMatrix(A, STAGES, STAGES, "A");

    double *B = (double *)Calloc(STAGES * STAGES, sizeof(double));
    double tempB[STAGES * STAGES] = { b11, b12, b21, b22 };
    initMatrixByRowWithValuesFromVector(B, STAGES, STAGES, tempB, STAGES * STAGES);
    //printDMatrix(B, STAGES, STAGES, "B");

    double *R = (double *)Calloc(STAGES * STAGES, sizeof(double));
    double tempR[STAGES * STAGES] = { 0.0f, 0.0f, r21, 0.0f };
    initMatrixByRowWithValuesFromVector(R, STAGES, STAGES, tempR, STAGES * STAGES);
    //printDMatrix(R, STAGES, STAGES, "R");

    /******************************* 
     *  Compute the solution
     * ****************************/
    double h = (t_span[1] - t_span[0]) / N;
    double *t = linspace(t_span[0], t_span[1], N + 1);
    int t_size = N + 1;
    int n = 1;

    int s = STAGES; // Number of stages
    int d1 = y0_size; // Dimension of the problem
    int Y_rows = s * d1, Y_cols = N;
    double *Y = zerosMatrixD(Y_rows, Y_cols);
    //fprintf(stdout, "\nh: %lf\ns: %d\nd1: %d\n", h, s, d1);
    //printDVector(t, N + 1, "t");
    //printDMatrix(Y, Y_rows, Y_cols, "Y");

    /************************************************
     * Runge-Kutta of order four to initialize stages 
     * **********************************************/
    double *FYiRK;
    int FYiRK_size;
    for (int i = 0; i < s; i++) {
        FYiRK = RungeKutta4th(c[i] * h, t[0], y0, y0_size, L, Lsize, &FYiRK_size);
        //printDVector(FYiRK, FYiRK_size, "FYiRK");
        for (int k = 0; k < d1; k++) {
            Y[(n - 1) * Y_rows + (i * d1 + k)] = FYiRK[k];
        }
    }
    //printDMatrix(Y, Y_rows, Y_cols, "Y");
    //printDVector(FYiRK, FYiRK_size, "FYiRK");
    //saveMatrixInFile("../parallel/good.txt", Y, Y_rows, Y_cols);

    int y_rows = d1, y_cols = N + 1;
    double *y = zerosMatrixD(y_rows, y_cols);
    // Fill the firt column of the matrix y
    for (int k = 0; k < d1; k++) {
        y[0 * y_rows + k] = y0[k];
    }
    // Fill the second column of the matrix y
    for (int k = 0; k < d1; k++) {
        y[(n) * y_rows + k] = Y[(n - 1) * Y_rows + ((s - 1) * d1 + k)];
    }
    //printDMatrix(y, y_rows, y_cols, "y");
    //saveMatrixInFile("../parallel/good.txt", y, y_rows, y_cols);

    int Fnm1_size = s * d1;
    double *Fnm1 = zerosD(Fnm1_size);
    double *Yi = zerosD(d1);

    for (int i = 0; i < s; i++) {
        //printDMatrix(Y, Y_rows, Y_cols, "Y");
        for (int k = 0; k < d1; k++) {
            //printf("d1: %d\n", d1);
            Yi[k] = Y[(n - 1) * Y_rows + (i * d1 + k)];
            //int x = (n) * Y_rows + (i * d1 + k);
            //printf("Y[%d]: %lf\n", x, Y[x]);
        }
        int FYi_size;
        double *FYi = Sherratt(Yi, d1, L, Lsize, &FYi_size);
        for (int k = 0; k < d1; k++) {
            Fnm1[i * d1 + k] = FYi[k];
        }
    }
    //printDVector(Yi, d1, "Yi");
    //printDVector(Fnm1, Fnm1_size, "Fnm1");
    //saveVectorsInFile("../parallel/good.txt", 2, Fnm1, Fnm1_size, Yi, d1, (void *)0);

    for (n = 1; n < N; n++) {
        //fprintf(stdout, "\nn: %d \n", n);
        for (int i = 0; i < s; i++) {
            for (int k = 0; k < d1; k++) {
                for (int j = 0; j < s; j++) {
                    Y[(n) * Y_rows + (i * d1 + k)] += B[j * STAGES + i] * Y[(n - 1) * Y_rows + ((j) * d1 + k)] + h * A[j * STAGES + i] * Fnm1[j * d1 + k];
                }
            }
        }
        //printDMatrix(Y, Y_rows, Y_cols, "Y");

        Fnm1 = zerosD(Fnm1_size);
        Yi = zerosD(d1);
        for (int i = 0; i < s; i++) {
            for (int k = 0; k < d1; k++) {
                Yi[k] = Y[(n) * Y_rows + (i * d1 + k)];
            }
            int FYi_size;
            double *FYi = Sherratt(Yi, d1, L, Lsize, &FYi_size);
            for (int k = 0; k < d1; k++) {
                Fnm1[i * d1 + k] = FYi[k];
            }
        }
        //printf("\nn: %d", n);
        //printDVector(Yi, d1, "Yi");
        //printDVector(Fnm1, Fnm1_size, "Fnm1");

        for (int k = 0; k < d1; k++) {
            y[(n + 1) * y_rows + k] = Y[(n) * Y_rows + ((s - 1) * d1 + k)];
        }
    }
    //printDMatrix(Y, Y_rows, Y_cols, "Y");
    //printDMatrix(y, y_rows, y_cols, "y");
    //printDVector(Yi, d1, "Yi");
    //printDVector(Fnm1, Fnm1_size, "Fnm1");

    //fprintf(stdout, "Here\n");
    double *yT = zerosD(d1);
    int yT_size = d1;
    for (int i = 0; i < d1; i++) {
        //fprintf(stdout, "y[(N) * y_cols + i]: %f\n", y[(N) * y_rows + i]);
        yT[i] = y[(N) * y_rows + i];
    }
    //printDVector(yT, yT_size, "yT");

    // After all the calculation, collecting and return the results
    collect_result->y = y;
    collect_result->y_rows = y_rows;
    collect_result->y_cols = y_cols;

    collect_result->t = t;
    collect_result->t_size = t_size;

    collect_result->yT = yT;
    collect_result->yT_size = yT_size;
}