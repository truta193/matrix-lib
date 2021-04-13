#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
typedef enum {false, true} bool;

/*
*    A small library for matrix operations
*
*    Author:         https://github.com/truta193
*    Last update:    06/04/2021
*/

/*  Naive matrix multiplication
*
*   If <output> is specified, it will be used as the resulting matrix, returning
*   NULL; else (NULL) a chunk of memory will be malloc'd and the resulting 
*   matrix pointer will be returned by the function.
*
*   Matrix sizes:
*   mat1        = M(m,k)
*   mat2        = M(k,n)
*   mat3/output = M(m,n)
*/

void *matMultiplication(void *mat1, void *mat2, int m, int n, int k, void *output) {
    void *mat3;
    if (output == NULL || output == mat1 || output == mat2) mat3 = malloc(sizeof(float)*m*n); else mat3 = output;
    //mat3 = malloc(sizeof(float)*m*n);
    float *mem1 = (float*)mat1;
    float *mem2 = (float*)mat2;
    float *mem3 = (float*)mat3;

    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) {
        mem3[i*n+j] = 0;
        for (int l = 0; l < k; l++) mem3[i*n+j] += mem1[i*k+l]*mem2[l*n+j];
    };

    if (output == NULL) return mat3;
    if (output == mat1 || output == mat2) {
        memcpy(output, mat3, sizeof(float)*m*n);
        free(mat3);
        return NULL;
    };
    return NULL;
};

/*  Matrix addition.
*   
*   Matrices must be of equal size (mat1,mat2 = M(m,n)). If <output> is 
*   specified, it will be used as the resulting matrix, else (NULL) a chunk of 
*   memory will be malloc'd and the resulting matrix pointer will be returned by
*   the function.
*/
void *matAddition(void *mat1, void *mat2, int m, int n, void *output) {
    void *mat3;
    if (output == NULL) mat3 = malloc(sizeof(float)*m*n); else mat3 = output;

    float *mem1 = (float*)mat1;
    float *mem2 = (float*)mat2;
    float *mem3 = (float*)mat3;   
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) mem3[i*n+j] = mem1[i*n+j] + mem2[i*n+j];
    return mat3;
};

/*  Matrix scaling
*
*   Scales the matrix by <scalar>. If <output> is specified, it will be used as
*   the resulting matrix, else (NULL) a chunk of memory will be malloc'd and the
*   resulting matrix pointer will be returned by the function.
*/
void *matScalar(void *mat, int m, int n, float scalar, void *output) {
    void *mat2;
    if (output == NULL) mat2 = malloc(sizeof(float)*m*n); else mat2 = output;

    float *mem = (float*)mat;
    float *mem2 = (float*)mat2;
    for (int i = 0; i < m*n; i++) mem2[i] = mem[i]*scalar;
    return mat2;
};

/*  Matrix transposing
*
*   Transposes the matrix. If <output> is specified, it will be used as the 
*   resulting matrix, returning NULL; else (NULL) a chunk of memory will be 
*   malloc'd and the resulting matrix pointer will be returned by the function.
*/
void *matTranspose(void *mat, int m, int n, void *output) {
    void *mat2;
    if (output == NULL || output == mat) mat2 = malloc(sizeof(float)*m*n); else mat2 = output;
    float *mem = (float*)mat;
    float *mem2 = (float*)mat2;

    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
        mem2[i*m+j] = mem[j*n+i];

    if (output != NULL) {
        memcpy(output, mat2, sizeof(float)*m*n);
        free(mat2);
        return NULL;
    }
    return mat2;
};
