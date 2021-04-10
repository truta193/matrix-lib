#include <stdio.h>
#include <memory.h>
#include <stdlib.h>
typedef enum {false, true} bool;

bool lsame(char a, char b);
bool lsame(char a, char b) {
    if (a == b) return true;
    return false;
};
#define max(x,y) ((x)>(y)?(x):(y))

//From LAPACK, translated from Fortran to help me understand it better but never actually tested if this works
void sgemm(char transa, char transb, int m, int n, int k,float alpha, float a[4][4], int lda, float b[4][4], int ldb,float beta, float *c[4][4], int ldc) {
    //FOR SIMPLE MULTIPL, BETA = 0, ALPHA = 1, AND MATRICES ARE NOT TRANSPOSED
    //LDA, LDB only used if this is called on a submatrix of a bigger matrix
    
    float temp;
    float zero = 0.0f;
    float one = 1.0f;
    bool nota, notb;
    int i,info,j,l,nrowa,nrowb;

    nota = lsame(transa, 'n');
    notb = lsame(transb, 'n');
    if (nota) nrowa = m; else nrowa = k;
    if (notb) nrowb = k; else nrowb = n;

    info = 0;
    if (!nota && !lsame(transa, 'c') && !lsame(transa, 't')) info = 1; 
        else if (!notb && !lsame(transb, 'c') && !lsame(transb, 't')) info = 2;
            else if (m < 0) info = 3; 
                else if (n < 0) info = 4;
                    else if (k < 0) info = 5;
                        else if (lda < max(1,nrowa)) info = 8;
                            else if (ldb < max(1, nrowb)) info = 10;
                                else if (ldc < max(1,m)) info = 13;
    if (info != 0) {
        // CALL xerbla('SGEMM ', info); (Fortran)
        return;
    };

    if (m == 0 || n == 0 || ((alpha == zero || k == 0 ) && beta == one)) return;
    
    if (alpha == zero) {
        if (beta == zero) 
            for (j = 1; j <=n; j++) for (i = 1; i <= m; i++) *c[i][j] = zero;
        else 
            for (j = 1; j <=n; j++) for (i = 1; i <= m; i++) *c[i][j] *= beta;
    };

    if (notb) {
        if (nota) {
            for (j = 1; j <= n; j++) {
                if (beta == zero) {
                    for (i = 1; i <= m; i++) *c[i][j] = zero; 
                }else if (beta != one) {
                        for (i = 1; i <= m; i++) *c[i][j] *= beta;
                };
                for (l = 1; l <= k; l++) {
                    temp = alpha*b[l][j];
                    for (i = 1; i <= m; i++) *c[i][j] += temp*a[i][l];
                };
            };
        } else {
            for (j = 1; j <= n; j++) for (i = 1; i <=m; i++) {
                temp = zero;
                for (l = 1; i <= k; l++) temp += a[j][i] * b[l][j];
                if (beta == zero) {
                        *c[i][j] = alpha*temp; 
                    } else {
                        *c[i][j] = alpha*temp + beta*(*c[i][j]);
                    };
            };
        };
    } else {
        if (nota) {
            for (j = 1; j <= n; j++) {
                if (beta == zero) {
                    for (i = 1; i <= m; i++) *c[i][j] = zero;
                } else if (beta != one) {
                    for (i = 1; i <= m; i++) *c[i][j] *= beta;
                };
                for (l = 1; l <= k; l++) {
                    temp = alpha*b[j][l];
                    for (i = 1; i <= m; i++) *c[i][j] += temp*a[i][l];
                };
            };
        } else {
            for ( j = 1; j <= n; j++) for (i = 1; i <= m; i++) {
                temp = zero;
                for (l = 1; l <= k; l++) temp += a[l][i]*b[j][l];
                if (beta == zero) {
                    *c[i][j] = alpha*temp; 
                } else {
                    *c[i][j] = alpha*temp + beta*(*c[i][j]);
                };
            };
        };
    };  
};

typedef struct matrF {
    int rows;
    int cols;
    float *data;
} matrF;

float* test(float *mat, int n, int m) {
    void *temp = malloc(sizeof(float)*6);
    float *rat = (float*)temp;
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < m; j++) rat[i*m+j] = 7.0f;
    return rat; 
};

#define tst(matrix, rows, cols) test((matrix)[0], (rows), (cols))

matrF create(int rows, int cols, float *elems) {
    void *buffer = malloc(sizeof(float)*rows*cols);
    memcpy(buffer, elems, sizeof(float)*rows*cols);
    matrF temp = {
        .cols = cols,
        .rows = rows,
        .data = buffer
    };
    return temp;
};

void *readMat(const void *mat, int rows, int cols) {
    for (int i = 0; i < rows*cols; i++) printf("%f\n", ((float*)mat)[i]);
    float a[4][4] = {
        1,1,1,1,
        1,1,1,1,
        1,1,1,1,
        1,1,1,1
    };
    void *temp = malloc(sizeof(float)*rows*cols);
    memcpy(temp, a, sizeof(float)*4*4);
    return temp;
};


void *matMultiplication(void *mat1, void *mat2, int m, int n, int k, void *output) {
    void *mat3;
    if (output == NULL) mat3 = malloc(sizeof(float)*m*n); else mat3 = output;

    float *mem1 = (float*)mat1;
    float *mem2 = (float*)mat2;
    float *mem3 = (float*)mat3;
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) {
        mem3[i*n+j] = 0;
        for (int l = 0; l < k; l++) mem3[i*n+j] += mem1[i*k+l]*mem2[l*n+j];  
    };
    return mat3;
};

void *matAddition(void *mat1, void *mat2, int m, int n, void *output) {
    void *mat3;
    if (output == NULL) mat3 = malloc(sizeof(float)*m*n); else mat3 = output;

    float *mem1 = (float*)mat1;
    float *mem2 = (float*)mat2;
    float *mem3 = (float*)mat3;   
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) mem3[i*n+j] = mem1[i*n+j] + mem2[i*n+j];
    return mat3;
};

void *matScalar(void *mat, int m, int n, float scalar, void *output) {
    void *mat2;
    if (output == NULL) mat2 = malloc(sizeof(float)*m*n); else mat2 = output;

    float *mem = (float*)mat;
    float *mem2 = (float*)mat2;
    for (int i = 0; i < m*n; i++) mem2[i] = mem[i]*scalar;
    return mat2;
};

void *matTranspose(void *mat, int m, int n, void *output) {
    void *mat2;
    if (output == NULL) mat2 = malloc(sizeof(float)*m*n); else mat2 = output;

    float *mem = (float*)mat;
    float *mem2 = (float*)mat2;
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
        mem2[i*m+j] = mem[j*n+i];
    return mat2;
};

//Testing
int main() {
    float mat1[2][4] = {
        2, 7, 5, 1, 
        10, 9, 25, 7
    };
    float mat2[4][3] = {
        5, 9, 1,
        45, 6, 8,
        13, 12, 8,
        7, 11, 15
    };
    float mat11[2][3] = {0};
    float mat3[2][3] = {0};
    float matadd[2][3] = {
        100,100,100,
        100,100,100
    };

    matMultiplication(mat1, mat2, 2, 3, 4, mat3);

    matScalar(mat3, 2, 3, 2.0f, mat3);
    for (int i = 0; i < 2; i++){
         for (int j = 0; j < 3; j++) printf("%.2f ", mat3[i][j]);
         printf("\n");
    };
    
    matAddition(mat3, matadd, 2, 3, mat3);
    for (int i = 0; i < 2; i++){
         for (int j = 0; j < 3; j++) printf("%.2f ", mat3[i][j]);
         printf("\n");
    };

    
    float *temp = matTranspose(mat3, 2, 3, NULL);
    for (int i = 0; i < 3; i++){
         for (int j = 0; j < 2; j++) printf("%.2f ", temp[i*2+j]);
         printf("\n");
    };

    return 0;
}