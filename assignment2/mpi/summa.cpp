#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>

int rank;

int number_of_processes;

int m;
int n;
int k;

void initMatrices(double pha[], double phb[], double phc[], int m_ar, int m_br, MPI_Comm comm_cart){
    int coords[2];
    MPI_Cart_coords(comm_cart, rank, 2, coords);

    int Nj = coords[0];
    int Ni = coords[1];

    long int i,j;

    double f = 1.0;
    for(i=0; i< m_ar; i++){
        for(j=0; j< m_br; j++) {
            //pha[i*m_br + j] = (Ni*m_ar + i) == (Nj*m_ar + j) ? 1.0 : 0.0;
            pha[i*m_br + j] = f;
            f++;
        }
    }

    for(i=0; i< m_br; i++){
        for(j=0; j< m_ar; j++) {
            phb[i*m_ar + j] = (Ni*m_ar + i) == (Nj*m_ar + j) ? 1.0 : 0.0;
            //phb[i*m_ar + j] = (double)(i+1);
        }
    }

    for(i=0; i< m_ar; i++){
        for(j=0; j< m_ar; j++) {
            phc[i*m_ar + j] = 0;
        }
    }
}
void multMatricesLineByLine(int m_ar, int m_br, const int kk, double * pha, double * phb, double * phc){

    for(int i=0; i<m_ar; i++){
        for(int k=0; k<m_br; k++){
            for(int j=0; j<m_ar; j++){
                phc[i*m_ar+j] += (pha[i*m_ar+k] * phb[k*m_br+j]);
            }
        }
    }
}
/*
void multiplyMatrixNaive(const int m, const int n, const int k,const double *A, const double *B, double *C) {

    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < k; ++i) {
            C[j*k + i] = 0.0;
            for (int l = 0; l < n; ++l) {
                C[j*k + i] += A[j*n + l] * B[l*k + i];
            }
        }
    }
}
*/
void sumMatrix(const int m, const int n, double *A, double *B, double *C) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int idx = i*m + j;
            C[idx] = A[idx] + B[idx];
        }
    }
}

void gatherMatrix(int mb, int nb, double *A_loc, int m_a, int n_a, double *A_glob) {
    double * A_tmp = NULL;
    if (rank == 0) {
        A_tmp = (double *) calloc(m_a * n_a, sizeof(double));
    }

    MPI_Gather(A_loc, mb*nb, MPI_DOUBLE, A_tmp, mb*nb, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // only rank 0 has something to do, others are free
    if (rank != 0) return;

    int nblks_m = m_a / mb;
    int nblks_n = n_a / nb;
    int idx = 0;
    for (int blk_i = 0; blk_i < nblks_m; ++blk_i) {
        for (int blk_j = 0; blk_j < nblks_n; ++blk_j) {

            // position in global matrix where current block should start
            int blk_start_row = blk_i * mb;
            int blk_start_col = blk_j * nb;

            for (int i = 0; i < mb; ++i) {
                for (int j = 0; j < nb; ++j) {
                    A_glob[(blk_start_row + i)*n + (blk_start_col + j)] = A_tmp[idx];
                    idx++;
                }
            }

        }
    }

    free(A_tmp);
}

void printMatrix(const int rows, const int cols, const double *matr) {
    //printf("%i\n", rank);
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            printf("%6.2f", matr[j*cols + i]);
        }
        printf("\n");
    }
}

void distributeSquareMatrix(double * matrix, int size, double * local_matrix) {
    int number_blocks = size / number_of_processes;
    int block_size = size / number_blocks;

    double * temp_matrix = NULL;

    if(rank == 0){
        temp_matrix = (double * ) malloc(size * size * sizeof(double));


        int i = 0,j = 0;
        for (int block_i = 0; block_i < number_blocks; block_i++) {
            for (int block_j = 0; block_j < number_blocks; block_j++) {
                int block_row = block_i * block_size;
                int block_col = block_j * block_size;
                for(int local_i = 0; local_i < block_size; local_i++){
                    for(int local_j = 0; local_j < block_size; local_j++){
                        //temp_matrix[pos] = matrix[(block_row + local_i) * size + (block_col + local_j)];
                        //pos++;
                    }
                }
            }
        }
    }
    std::cout << "ORIGINAL" << std::endl;
    printMatrix(size, size, matrix);

    std::cout << "TEMP" << std::endl;
    printMatrix(size, size, temp_matrix);

    //MPI_Scatter(temp_matrix, size * size, MPI_DOUBLE, local_matrix, block_size * block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);



/*
    double * local_matrix = (double * ) malloc(block_size * block_size * sizeof(double));


    //
    for (int block_i = 0; block_i < number_blocks; block_i++) {
        for (int block_j = 0; block_j < number_blocks; block_j++) {

            for(int local_i = 0; local_i < block_size; local_i++){
                for(int local_j = 0; local_j < block_size; local_j++){

                }
            }
        }
    }
*/
}

double validate(const int m, const int n, const double *Csumma, double *Cnaive) {

    double eps = 0.0;
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int idx = i*n + j;
            Cnaive[idx] = fabs(Cnaive[idx] - Csumma[idx]);
            if (eps < Cnaive[idx]) {
                eps = Cnaive[idx];
            }
        }
    }

    return eps;
}

void getTimes(double start_time, double end_time){

    double t = end_time - start_time;
    //std::cout << "rank " << rank << " took " << t << " seconds" << std::endl;

    double max_time = 0.0;

    MPI_Reduce(&t, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank == 0){
      std::cout << "SUMMA took " << max_time << " seconds" << std::endl;
    }

    double * times = NULL;
    if (rank == 0) {
        times = (double *) malloc(sizeof(double) * number_of_processes);
    }

    MPI_Gather(&t, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0){
        std::ostringstream stringStream;
        stringStream << "log_" << number_of_processes << ".txt";
        const std::string tmp = stringStream.str();
        const char* cstr = tmp.c_str();
        std::ofstream logFile;
        logFile.open(cstr, std::ofstream::app | std::ofstream::ate);
        for(int i = 0; i < number_of_processes; i++){
            logFile << m << ";" << i << ";" << times[i] << std::endl;
            std::cout << "Rank  " << i << " took " << times[i] << " seconds."<< std::endl;
        }
        logFile.close();
    }
}

void summa(MPI_Comm comm_cart, const int m_block, const int n_block, const int k_block, double A_local[], double B_local[], double C_local[]) {
    // determine my cart coords
    int coords[2];
    MPI_Cart_coords(comm_cart, rank, 2, coords);

    const int my_row = coords[0];
    const int my_col = coords[1];

    int belongs[2];

    // create row comms for A
    MPI_Comm row_comm;
    belongs[0] = 0;
    belongs[1] = 1;
    MPI_Cart_sub(comm_cart, belongs, &row_comm);

    // create col comms for B
    MPI_Comm col_comm;
    belongs[0] = 1;
    belongs[1] = 0;
    MPI_Cart_sub(comm_cart, belongs, &col_comm);

    /*int row_rank, col_rank;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_rank(col_comm, &col_rank);
    if(rank == 1) std::cout << "Rank: " << rank << "-->(" << my_col << "," << my_row << ") (" << col_rank << "," << row_rank << ")" << std::endl;
    //printf("Rank: %i-->(%i,%i)|(%i,%i)\n", rank, my_col, my_row, col_rank, row_rank);*/

    double * A_saved = (double *) calloc(m_block * n_block, sizeof(double));
    double * B_saved = (double *) calloc(n_block * k_block, sizeof(double));
    double * C_tmp = (double *) calloc(m_block * k_block, sizeof(double));

    memcpy(A_saved, A_local, m_block * n_block * sizeof(double));
    memcpy(B_saved, B_local, n_block * k_block * sizeof(double));

    int number_blocks = n / n_block;
    for(int broadcaster = 0; broadcaster < number_blocks; ++broadcaster){
        if (my_col == broadcaster) {
            memcpy(A_local, A_saved, m_block * n_block * sizeof(double));
        }

        MPI_Bcast(A_local, m_block * n_block, MPI_DOUBLE, broadcaster, row_comm);

        if (my_row == broadcaster) {
            memcpy(B_local, B_saved, n_block * k_block * sizeof(double));
        }

        MPI_Bcast(B_local, n_block * k_block, MPI_DOUBLE, broadcaster, col_comm);

        multMatricesLineByLine(m_block, n_block, k_block, A_local, B_local, C_tmp);

        sumMatrix(m_block, n_block, C_local, C_tmp, C_local);
    }
}

int compute(){
    //Assuming that the processes form a square
    int n_procs_row = sqrt(number_of_processes);
    int n_procs_col = n_procs_row;
    if (n_procs_col * n_procs_row != number_of_processes) {
        std::cerr << "number of proccessors must be a perfect square!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int n_dims = 2;
    int dims[n_dims] = {n_procs_row, n_procs_col};
    int periods[n_dims] = {0, 0};
    int repeat = 0;

    //create comm_groups
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, n_dims, dims, periods, repeat, &comm_cart);

    int m_block = m / n_procs_row;
    int n_block = n / n_procs_col;
    int k_block = k / n_procs_col;

    if (m_block * n_procs_row != m) {
        std::cerr << "m must be dividable by n_procs_row" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (n_block * n_procs_col != n) {
        std::cerr << "n must be dividable by n_procs_col" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    if (k_block * n_procs_col != k) {
        std::cerr << "k must be dividable by n_procs_col" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    double * A_local = (double *) calloc(m_block * n_block, sizeof(double));
    double * B_local = (double *) calloc(n_block * k_block, sizeof(double));
    double * C_local = (double *) calloc(m_block * k_block, sizeof(double));

    double * A = (double *) calloc(m * n, sizeof(double));
    double * B = (double *) calloc(n * k, sizeof(double));

    initMatrices(A_local, B_local, C_local, m_block, n_block, comm_cart);


    /**
    testing the scatter function
    */
    if(rank == 0){
        double * A_1 = (double *) calloc(n * k, sizeof(double));
        double * B_1 = (double *) calloc(n * k, sizeof(double));
        double * C_1 = (double *) calloc(n * k, sizeof(double));

        initMatrices(A_1, B_1, C_1, n, n, comm_cart);

        distributeSquareMatrix(A_1, n, C_1);
    }

    gatherMatrix(m_block, n_block, A_local, m, n, A);
    gatherMatrix(n_block, k_block, B_local, n, k, B);
    /*if (rank == 3) {
    std::cout << "A" << std::endl;
      printMatrix(m_block,m_block, A_local);
      std::cout << "B" << std::endl;
      printMatrix(m_block,m_block, B_local);
      std::cout << "C" << std::endl;
      printMatrix(m_block,m_block, C_local);
    }*/
    double start_time, end_time;

    start_time = MPI_Wtime();

    summa(comm_cart, m_block, n_block, k_block, A_local, B_local, C_local);

    end_time = MPI_Wtime();

    getTimes(start_time, end_time);

    double * C = (double *) calloc(m * n, sizeof(double));
    double * C_naive = (double *) calloc(m * n, sizeof(double));

    gatherMatrix(m_block, k_block, C_local, m, k, C);

    if (rank == 0) {

        multMatricesLineByLine(m, n, k, A, B, C_naive);

        double eps = validate(n, k, C, C_naive);
        if (eps > 1e-4) {
            std::cerr <<  "ERROR: Invalid matrix -> eps = " << eps << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        } else {
            std::cout << "Valid matrix" << std::endl;
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //std::cout << "Rank: " << rank << std::endl;

    /*for(int i = 0; i < 5; i++){
        int j = 1000;
        if(number_of_processes == 9){
            j = 999;
        }
        m = n = k = j * (i + 1);*/
        m = n = k = 9;

        if(rank == 0) {
            std::cout << "======================" << std::endl;
            std::cout << "n: " << n << ", m: " << m << ", k: " << k << std::endl;
        }
        compute();
    //}

    MPI_Finalize();
}

/* EXTRA CODE FOR TESTING

if (my_col == broadcaster) {
    memcpy(A_local, A_saved, m_block * n_block * sizeof(double));
}

if(rank == local){
    std::cout << "=======================" << std::endl;
    std::cout << "=======================" << std::endl;
    std::cout << (my_col == broadcaster ? "SENDING" : "RECEIVING") << std::endl;
    std::cout << "=======================" << std::endl;
    std::cout << "BEFORE A" << std::endl;
    printMatrix(m_block,m_block, A_local);
}
MPI_Bcast(A_local, m_block * n_block, MPI_DOUBLE, broadcaster, row_comm);
if(rank == local){
    std::cout << "AFTER A" << std::endl;
    printMatrix(m_block,m_block, A_local);
}

if (my_row == broadcaster) {
    memcpy(B_local, B_saved, n_block * k_block * sizeof(double));
}

if(rank == local){
    std::cout << "=======================" << std::endl;
    std::cout << "=======================" << std::endl;
    std::cout << (my_row == broadcaster ? "SENDING" : "RECEIVING") << std::endl;
    std::cout << "=======================" << std::endl;
    std::cout << "BEFORE B" << std::endl;
    printMatrix(m_block,m_block, B_local);
}
MPI_Bcast(B_local, n_block * k_block, MPI_DOUBLE, broadcaster, col_comm);
MPI_Barrier(col_comm);
if(rank == local){
    std::cout << "AFTER B" << std::endl;
    printMatrix(m_block,m_block, B_local);
}

multMatricesLineByLine(m_block, n_block, k_block, A_local, B_local, C_tmp);
sumMatrix(m_block, n_block, C_local, C_tmp, C_local);

if(rank == local){
    std::cout << "A" << std::endl;
    printMatrix(m_block,m_block, A_local);
    std::cout << "B" << std::endl;
    printMatrix(m_block,m_block, B_local);
    std::cout << "C" << std::endl;
    printMatrix(m_block,m_block, C_local);
}
*/
