/**
 * @file peerMethods.h
 * @author Vincenzo Iannucci
 * @brief The library provides an implementation for the main function for solving peer method.
 * @version 0.2
 * @date 01-30-2023 (MM-DD-YYYY)
 * @dir peerMethods/include
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#ifndef peerMethods_h
#define peerMethods_h

#ifdef __cplusplus
extern "C" {
#endif

/* Number of stages. */
#define STAGES 2

/* Model parameters declared as constants */
#define a 1.5 
#define B1 0.45 
#define B2 0.3611 
#define F 0.802 
#define H 0.802 
#define S 0.0002 
#define d 500 
#define D 0.802

/* Variable that need to be set in the calling method. */
extern int M;

/** @struct return_values
 *  @brief This struct contains the returning values from peer methods function.
 *  @var foreignstruct::yT 
 *  Member 'yT' represents the final solution of the PDE system.
 *  @var foreignstruct::y
 *  Member 'y'contains the solution at each stage of the PDE system.
 *  @var foreignstruct::t
 *  Member 't' represents the discretization of the temporal grid.
 */
typedef struct {
    double *yT;
    int yT_size;
    double *y;
    int y_rows;
    int y_cols;
    double *t;
    int t_size;
} return_values;

/**
 * @brief Initialize the struct return_values.
 * @param rv pointer to the struct return_values
*/
void initReturnStruct(return_values *rv);

/**
 * @brief Save the struct return_values in a file.
 * @param[in] fileName the name of the file
 * @param[out] rv pointer to the struct return
 * @return 0 if ok, 1 otherwise.
*/
int saveResultsInFile(const char* fileName, return_values result);

/**
 * @brief Spatial semi-discretization of the PDE (Partial Differential Equation).
 * Build the matrix L to perform the spatial semi-discretization of the ODE.
 * @param[out] L returning pointer to the matrix
 * @param[in] LSize return the size of the matrix
 * @param[in] Delta_x the value of the delta
*/
void defineLMatrix(double **L, int *LSize, double Delta_x);

/**
 * @brief Applies the Sherratt method.
 * @param[in] y0 pointer to the y0 vector
 * @param[in] y0Size size of the y0 vector
 * @param[in] L pointer to the matrix L
 * @param[in] LSize size of the matrix
 * @param[out] sherrattSize returing size of the vector calculated by the function
 * @return a pointer the resulting vector after applying the Sherratt method.
*/
double *Sherratt(const double *y0, int y0Size, const double *L, int Lsize, int *sherrattSize);

/**
 * @brief Returns the numerical approximation of the solution of the equation at the
 * next time step, given the current one, by using Runge Kutta 4th order method.
 * @param[in] h constant grid spacing (the space between the values of the discrete grid)
 * @param[in] t0 starting time
 * @param[in] y0 vector representing value of the y function in y(t0)
 * @param[in] y0Size size of the y0 vector
 * @param[in] L pointer to the matrix L
 * @param[in] LSize size of the matrix
 * @param[out] ySize size of the resulting vector
 * @return a vector that approximates the value of the solution at the next time step t0 + h, y(t0 + h)
*/
double *RungeKutta4th(double h, double t0, const double *y0, int y0Size, const double *L, int Lsize, int *ySize);

/**
 * @brief Time discretization using peer methods.
 * The method solve a large system of ODEs (Ordinary Differential Equations) by using peer methods.
 * @param[in] N the size of the time grid
 * @param[in] t_span an array representing the time grid itself
 * @param[in] t_span_size the spatial dimension of the time grid
 * @param[in] L pointer to the matrix L deriving from the spatial semi-discretization
 * @param[in] LSize size of the matrix
 * @param[in] y0 pointer to the y0 vector, containing the initial conditions
 * @param[in] y0Size size of the y0 vector
 * @param[out] collect_result size of the result vector
 * @return a struct containing.
*/
void fPeerClassic_twoStages(int N, double *t_span, int t_span_size, const double *L, int Lsize, const double *y0, int y0_size, return_values *collect_result);

/********************************************************************************
 *                              Wrapper functions
 *******************************************************************************/

/**
 * @brief Function wrapper for malloc() function.
 * @param[in] size Size of the memory allocated
 * @return a pointer to the allocated memory
*/
void *Malloc(size_t size);

/**
 * @brief Function wrapper for calloc() function.
 * @param[in] nmemb number of elements to allocate
 * @param[in] size Size of the memory allocated
 * @return a pointer to the allocated memory
*/
void *Calloc(size_t nmemb, size_t size);

/**********************************************************************************
 *                          Initialization functions
 **********************************************************************************/

/**
 * @brief Initialize a vector with random values. NOTE: the seed must be initialized
 * in the calling method.
 * @param[in] vector pointer to the vector
 * @param[in] N dimension of the vector
*/
void initializeRandomVector(double *vector, int N);

/**
 * @brief Initialize a matrix with random values. NOTE: the seed must be initialized
 * in the calling method.
 * @param[in] matrix pointer to the first element of the matrix
 * @param[in] M number of rows
 * @param[in] N number of columns
*/
void initializeRandomMatrix(double *matrix, int M, int N);

/**
 * @brief Using a vector to initialize the matrix. The matrix and the vector must have the
 * same dimension. For example, if the matrix A has 3 x 4 elements, the vector B must have 12
 * elements.
 * @param[in out] matrix pointer to the matrix
 * @param[in] M rows of the matrix
 * @param[in] N columns of the matrix
 * @param[in] vector pointer to the vector
 * @param[in] vector_size size of the vector (must be equal to M x N)
 * @return 0 if ok, 1 otherwise
*/
int initMatrixByRowWithValuesFromVector(double *matrix, int M, int N, double *vector, int vector_size);

/**
 * @brief Using a vector to initialize another vector. The vectors must have the
 * same dimension.
 * @param[in out] newVector pointer to the new vector
 * @param[in] oldVector pointer to the old vector
 * @param[in] n size of the vectors
*/
void initVectorWAnotherVector(double *newVector, double *oldVector, int n);

/*************************************************
 *  Free all the memory dynamically allocated
 ***********************************************/
/**
 * @brief Free a variable number of pointers.
 * @param[in] arg1 pointer
*/
void freeEverything(void *arg1, ...);

/**********************************************************************************
 *                      MATLAB functions written in C
 **********************************************************************************/

/**
 * @brief Constant times a vector plus a vector.
 * @param N number of elements in input vector(s)
 * @param alpha scalar
 * @param X double precision array
 * @param incX storage spacing between elements of X
 * @param Y double precision array
 * @param incY storage spacing between elements of Y
 */
void peerMethodsDaxpy(const int N, const double alpha, const double *X, const int incX, double *Y, const int incY);

/**
 * @brief Scales a vector by a constant.
 * @param N  number of elements in input vector(s)
 * @param alpha scalar
 * @param X double precision array
 * @param incX storage spacing between elements of X
 */
void peerMethodsDscal(const int N, const double alpha, double *X, const int incX);

/**
 * @brief Performs matrix-vector operations y := y + alpha*A*x.
 * @param M specifies the number of rows of the matrix A
 * @param N specifies the number of columns of the matrix A
 * @param alpha scalar
 * @param A double precision array containing the matrix coefficients
 * @param lda specifies the second dimension of A
 * @param X double precision array containing the elements of the vector x
 * @param incX storage spacing between elements of X
 * @param Y double precision array resulting elements from matrix-vector multiplication
 * @param incY storage spacing between elements of Y
 */
void peerMethodsDgemv(const int M, const int N, const double alpha, const double *A, const int lda, const double *X, 
                      const int incX, double *Y, const int incY);

/**
 * @brief Provide the discretization of an interval starting with first
 * and ending with last.
 * @param first starting of the interval
 * @param last ending of the interval
 * @param step the spacing between each value of the array
 * @param N the size of the final array returned
 * @return pointer to an array of double, representing the discretized interval
 */
double *intervalDiscretization(double first, double last, double step, int *N);

/**
 * @brief Create identity matrix, a matrix in which the values on the 
 * main diagonal have value 1. 
 * @param N size of the output array
 * @return pointer to the new array
*/
double *eyeD(int N);

/**
 * @brief Create array of all ones.
 * @param N
 * @return pointer to the new array
*/
double *onesD(int N);

/**
 * @brief Create array of all zeros.
 * @param N dimension of the array
 * @return pointer to the new array
*/
double *zerosD(int N);

/**
 * @brief Create matrix of all zeros.
 * @param M the number of rows
 * @param N the number of columns
 * @return pointer to the new matrix allocated by rows
*/
double *zerosMatrixD(int M, int N);

/**
 * @brief Create diagonal matrix with all the elements of 
 * vector on the k-th diagonal.
 * @param vector pointer to the vector
 * @param size size of the vector
 * @param k the number of diagonal
 * @param matrix_size the size of the output matrix
 * @return pointer to the new matrix allocated by rows
*/
double *diagD(double *vector, int size, int k, int *matrix_size);

/**
 * @brief Packing three square matrices side by side with the 
 * same dimension into a new big one.
 * @param n number of rows
 * @param A pointer the first matrix
 * @param B pointer the second matrix
 * @param C pointer the third matrix
 * @return pointer to the new matrix
*/
double *packThreeMatrices(int n, double *A, double *B, double *C);

/**
 * @brief Create a square block diagonal matrix made up of three matrices.
 * @param n number of rows
 * @param A pointer the first matrix
 * @param B pointer the second matrix
 * @param C pointer the third matrix
 * @param blockSize the size of the output matrix
 * @return pointer to the output matrix
*/
double *threeBlockDiagD(int n, double *A, double *B, double *C, int *blckSize);

/**
 * @brief Packing three vectors side by side into one.
 * @param n number of rows
 * @param A pointer the first matrix
 * @param B pointer the second matrix
 * @param C pointer the third matrix
 * @param newDimension the size of the output vector
 * @return pointer to the new vector
*/
double *packThreeVectors(int n, double *A, double *B, double *C, int *newDimension);

/**
 * @brief Generate linearly spaced vector. 
 * @example linspace(double x1, double x2, int n) generates n points. The spacing between the points is (x2-x1)/(n-1).
 * @return pointer to the new vector
*/
double *linspace(double x1, double x2, int n);

#ifdef __cplusplus
}
#endif
#endif // !peerMethods_h