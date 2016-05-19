#include <stdio.h>
#include "mpi.h"


float internalProduct(float *v1, float *v2, int col) {
    int i;
    float sum = 0.0;

    for (i = 0; i < col; i++)
        sum += v1[i] * v2[i];

    return (sum);

}

void matrixFilling(float *v, int nLin, int nCol, int tlin) {
    int i, j;

    // matriz identidade
    for (i = 0; i < nLin; i++) {
        for (j = 0; j < nCol; j++)
            v[i * nCol + j] = 0.0;
        v[i * nCol + i] = 1.0;
    }

    for (i = nLin; i < tlin; i++) {
        for (j = 0; j < nCol; j++)
            v[i * nCol + j] = 0.0;
    }

}


int main(int argc, char **argv) {
    int rank, size, error;
    float *v1, *v2, *sv1, *my_r, *r;
    int k, nLines, my_lin, maisum;
    int lin, col;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (rank == 0) {
        printf("Dimensoes: lins cols ? ");
        scanf("%d %d", &lin, &col);
        maisum = (lin % size) == 0 ? 0 : 1;
        nLines = ((lin / size) + maisum) * size;
        v1 = new float[nLines * col];
        v2 = new float[col];
        r = new float[nLines];                // vector resultado
        matrixFilling(v1, lin, col, nLines);
        // preenche vector
        for (k = 0; k < col; k++)
            v2[k] = k;
    }

    // broadcast de lin, col a todos os processos
    //
    //


    // reservar espaco para os dados, o ultimo processador podera ficar com zeros no fim
    maisum = (lin % size) == 0 ? 0 : 1;
    my_lin = (lin / size) + maisum;
    sv1 = new float[my_lin * col];

    if (rank > 0)
        v2 = new float[col];


    // distribui a matriz
    //erro =
    printf("processo %d recebeu %d elementos, erro=%d\n", rank, my_lin * col, error);

    // distribui vector
    //erro =
    printf("processo %d recebeu vector de %d elementos, erro=%d\n", rank, col, error);

    // calcula produto matrix-vector
    my_r = new float[my_lin];

    for (k = 0; k < my_lin; k++)
        my_r[k] = internalProduct(&sv1[k * col], v2, col);

    // processo 0 obtem resultado
    //erro =

    if (rank == 0) {
        printf("Vector Resultado\n");
        for (k = 0; k < lin; k++)
            printf("%.1f ", r[k]);
        printf("\n");
    }

    MPI_Finalize();

    delete[] v1;
    delete[] v2;
    delete[] sv1;
    delete[] my_r;
    if (r != NULL) delete[] r;

    return 0;
}
