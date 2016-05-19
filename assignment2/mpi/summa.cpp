#include <stdio.h>
#include <math.h>
#include <mpi.h>

int rank;

int number_of_processes;

int m;
int n;
int k;

void initMatrices(double pha[], double phb[], double phc[], int m_ar, int m_br){
    long int i,j;

    for(i=0; i< m_ar; i++){
        for(j=0; j< m_br; j++) {
            pha[i*m_br + j] = (double)1.0;
        }
    }

    for(i=0; i< m_br; i++){
        for(j=0; j< m_ar; j++) {
            phb[i*m_ar + j] = (double)(i+1);
        }
    }

    for(i=0; i< m_ar; i++){
        for(j=0; j< m_ar; j++) {
            phc[i*m_ar + j] = 0;
        }
    }
}

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

void sumMatrix(const int m, const int n, double *A, double *B, double *C) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            int idx = i*m + j;
            C[idx] = A[idx] + B[idx];
        }
    }
}

void gatherMatrix(int mb, int nb, double *A_loc, int m, int n, double *A_glob) {
    double * A_tmp;
    if (rank == 0) {
        A_tmp = new double[m * n];
    }

    MPI_Gather(A_loc, mb*nb, MPI_DOUBLE, A_tmp, mb*nb, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // only rank 0 has something to do, others are free
    if (rank != 0) return;

    int nblks_m = m / mb;
    int nblks_n = n / nb;
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
}

void printMatrix(const int rows, const int cols, const double *matr) {
    printf("%i\n", rank);
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            printf("%6.2f", matr[j*cols + i]);
        }
        printf("\n");
    }

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

void summa(MPI_Comm comm_cart, const int m_block, const int n_block, const int k_block, double A_local[], double B_local[], double C_local[]) {
    // determine my cart coords
    int coords[2];
    MPI_Cart_coords(comm_cart, rank, 2, coords);

    int col_cart = coords[0];
    int row_cart = coords[1];

    int belongs[2];

    // create row comms for A
    MPI_Comm row_comm;
    belongs[0] = 1;
    belongs[1] = 0;
    MPI_Cart_sub(comm_cart, belongs, &row_comm);

    // create col comms for B
    MPI_Comm col_comm;
    belongs[0] = 0;
    belongs[1] = 1;
    MPI_Cart_sub(comm_cart, belongs, &col_comm);

    /*int row_rank, col_rank;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_rank(col_comm, &col_rank);
    std::cout << "Rank: " << rank << "-->(" << col_cart << "," << row_cart << ") (" << col_rank << "," << row_rank << ")" << std::endl;
    //printf("Rank: %i-->(%i,%i)|(%i,%i)\n", rank, col_cart, row_cart, col_rank, row_rank);*/

    double A_saved[m_block * n_block];
    double B_saved[n_block * k_block];
    double C_tmp[m_block * n_block];

    memcpy(A_saved, A_local, m_block * n_block * sizeof(double));
    memcpy(B_saved, B_local, n_block * k_block * sizeof(double));

    for(int k_x = 0; k_x < n / n_block; ++k_x){
        if (col_cart == k_x) {
            memcpy(A_local, A_saved, m_block * n_block * sizeof(double));
        }

        MPI_Bcast(A_local, m_block * n_block, MPI_FLOAT, col_cart, row_comm);

        if (row_cart == k_x) {
            memcpy(B_local, B_saved, n_block * k_block * sizeof(double));
        }

        MPI_Bcast(B_local, n_block * k_block, MPI_FLOAT, row_cart, col_comm);

        multiplyMatrixNaive(m_block, n_block, k_block, A_local, B_local, C_tmp);

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

  double A_local[m_block * n_block];
  double B_local[n_block * k_block];
  double C_local[m_block * k_block];

  initMatrices(A_local, B_local, C_local, m_block, n_block);

  double A[m * n];
  double B[n * k];

  gatherMatrix(m_block, n_block, A_local, m, n, A);
  gatherMatrix(n_block, k_block, B_local, n, k, B);

  double start_time, end_time;

  start_time = MPI_Wtime();

  summa(comm_cart, m_block, n_block, k_block, A_local, B_local, C_local);

  end_time = MPI_Wtime();

  double t = end_time - start_time;
  double max_time = 0.0;

  MPI_Reduce(&t, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  double C[m * n];
  double C_naive[m * n];

  gatherMatrix(m_block, k_block, C_local, m, k, C);

  if (rank == 0) {
      std::cout << "SUMMA took " << max_time << " seconds" << std::endl;
      multiplyMatrixNaive(m, n, k, A, B, C_naive);

      double eps = validate(n, k, C, C_naive);
      if (eps > 1e-4) {
          fprintf(stderr, "ERROR: eps = %f\n", eps);
      } else {
          printf("SUMMA: OK: eps = %f\n", eps);
      }

  }

}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //std::cout << "Rank: " << rank << std::endl;

    //for(int i = 0; i < 6; i++){
        n = 180;
        m = 180;
        k = 180;

        //if(rank == 0) std::cout << "n: " << n << ", m: " << m << ", k: " << k << std::endl;
        compute();
    //}

    MPI_Finalize();
}
