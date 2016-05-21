#include <stdio.h>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <mpi.h>
#include <stdlib.h>
#include <unistd.h>

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

    for(i=0; i< m_ar; i++){
        for(j=0; j< m_br; j++) {
            pha[i*m_br + j] = (Ni*m_ar + i) == (Nj*m_ar + j) ? (Ni*m_ar + i) : 0.0;
        }
    }

    for(i=0; i< m_br; i++){
        for(j=0; j< m_ar; j++) {
            phb[i*m_ar + j] = (Ni*m_ar + i) == (Nj*m_ar + j) ? (Ni*m_ar + i) : 0.0;
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
    for (int j = 0; j < (rows > 10 ? 10 : rows); ++j) {
        for (int i = 0; i < (cols > 10 ? 10 : cols); ++i) {
            printf(" %3.f", matr[j*cols + i]);
        }
        printf("\n");
    }
}
void distributeSquareMatrix(double * matrix, int size, double * local_matrix) {
    //int number_blocks = size / number_of_processes;
    //int block_size = size / number_blocks;

    double * temp_matrix = NULL;

    if(rank == 0){
        temp_matrix = (double * ) malloc(size * size * sizeof(double));

/*
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
        }*/
    }
    std::cout << "ORIGINAL" << std::endl;
    printMatrix(size, size, matrix);

    std::cout << "TEMP" << std::endl;
    printMatrix(size, size, temp_matrix);

    //MPI_Scatter(temp_matrix, size * size, MPI_DOUBLE, local_matrix, block_size * block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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
void getTimes(double start_time, double end_time, double networkElapsedTime, double processingElapsedTime){

    double t = end_time - start_time;
    //std::cout << "rank " << rank << " took " << t << " seconds" << std::endl;

    double max_time = 0.0;

    MPI_Reduce(&t, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank == 0){
    //  std::cout << "SUMMA took " << max_time << " seconds" << std::endl;
    }

    double * times = NULL;
    double * times_net = NULL;
    double * times_proc = NULL;

    if (rank == 0) {
        times = (double *) malloc(sizeof(double) * number_of_processes);
        times_net = (double *) malloc(sizeof(double) * number_of_processes);
        times_proc = (double *) malloc(sizeof(double) * number_of_processes);
    }

    MPI_Gather(&t, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&networkElapsedTime, 1, MPI_DOUBLE, times_net, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&processingElapsedTime, 1, MPI_DOUBLE, times_proc, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0){
        std::ostringstream stringStream;
        stringStream << "log_" << number_of_processes << ".txt";
        const std::string tmp = stringStream.str();
        const char* cstr = tmp.c_str();
        std::ofstream logFile;
        logFile.open(cstr, std::ofstream::app | std::ofstream::ate);
        for(int i = 0; i < number_of_processes; i++){
            logFile << m << ";" << i << ";" << times[i] <<  ";" << times_net[i] << ";" << times_proc[i] << std::endl;
            //std::cout << "Rank  " << i << " took " << times[i] << " seconds."<<  std::endl;
        }
        //logFile << max_time << std::endl;
        logFile.close();
    }
    
    free(times);
    free(times_net);
    free(times_proc);
}

void summa(MPI_Comm comm_cart, const int m_block, const int n_block, const int k_block, double A_local[], double B_local[], double C_local[],
 double * networkElapsedTime, double * processingElapsedTime) {
	// determine my cart coords
    int coords[2];
    MPI_Cart_coords(comm_cart, rank, 2, coords);
    const int my_row = coords[0];
    const int my_col = coords[1];

    // create row comms for A
    MPI_Comm row_comm;
    int belongs_row[2] = {0,1};
    MPI_Cart_sub(comm_cart,belongs_row, &row_comm);

    // create col comms for B
    MPI_Comm col_comm;
    int belongs_col[2] = {1,0};
    MPI_Cart_sub(comm_cart, belongs_col, &col_comm);

	double * A_current, * B_current;
    double * A_saved = (double *) calloc(m_block * n_block, sizeof(double));
    double * B_saved = (double *) calloc(n_block * k_block, sizeof(double));

	*networkElapsedTime = 0;
	*processingElapsedTime = 0;

	// the important stuff.
    int number_blocks = n / n_block;
    for(int broadcaster = 0; broadcaster < number_blocks; ++broadcaster){
        A_current = my_col == broadcaster ? A_local : A_saved;
		B_current = my_row == broadcaster ? B_local : B_saved;

   		double network_start_time = MPI_Wtime();
        MPI_Bcast(A_current, m_block * n_block, MPI_DOUBLE, broadcaster, row_comm);
        MPI_Bcast(B_current, n_block * k_block, MPI_DOUBLE, broadcaster, col_comm);
    	*networkElapsedTime += (MPI_Wtime() - network_start_time);

		double processing_start_time = MPI_Wtime();
        multMatricesLineByLine(m_block, n_block, k_block, A_current, B_current, C_local);
		*processingElapsedTime += (MPI_Wtime() - processing_start_time);
    }
    
	free(A_saved);
	free(B_saved);
}

int ring_next(int current, int ring_size){
	return ((current+1) % ring_size);
}
int ring_prev(int current, int ring_size){
	return (current-1) < 0 ? (ring_size-1) : (current-1);
}


void summa_ring(MPI_Comm comm_cart, const int m_block, const int n_block, const int k_block, double A_local[], double B_local[], double C_local[],
 double * networkElapsedTime, double * processingElapsedTime) {
    // determine my cart coords
    int coords[2];
    MPI_Cart_coords(comm_cart, rank, 2, coords);
    const int my_row = coords[0];
    const int my_col = coords[1];

    // create row comms for A
    MPI_Comm row_comm;
    int belongs_row[2] = {0,1};
    MPI_Cart_sub(comm_cart,belongs_row, &row_comm);

    // create col comms for B
    MPI_Comm col_comm;
    int belongs_col[2] = {1,0};
    MPI_Cart_sub(comm_cart, belongs_col, &col_comm);


	int n_procs_row, n_procs_col; 
	MPI_Comm_size(row_comm, &n_procs_row);
	MPI_Comm_size(col_comm, &n_procs_col);

	double * A_current, * B_current;
    double * A_saved = (double *) calloc(m_block * n_block, sizeof(double));
    double * B_saved = (double *) calloc(n_block * k_block, sizeof(double));

	int row_next = ring_next(my_row, n_procs_row);
	int row_prev = ring_prev(my_row, n_procs_row);
	int col_next = ring_next(my_col, n_procs_col);
	int col_prev = ring_prev(my_col, n_procs_col);
	MPI_Status status;
	
	*networkElapsedTime = 0;
	*processingElapsedTime = 0;
	

	// the important stuff.
	int init_broadcaster = ((my_row+my_col)% n_procs_col); 
	int broadcaster = init_broadcaster;
	do{
		A_current = my_col == broadcaster ? A_local : A_saved;
		B_current = my_row == broadcaster ? B_local : B_saved;
		
		double network_start_time = MPI_Wtime();
        if(((my_row+my_col)% 2) == 0){

			MPI_Send(A_current, m_block * n_block, MPI_DOUBLE, col_next, 0, row_comm);
			MPI_Recv(A_saved, m_block * n_block, MPI_DOUBLE, col_prev, 0, row_comm, &status);	

			MPI_Send(B_current, m_block * n_block, MPI_DOUBLE, row_next, 0, col_comm);
			MPI_Recv(B_saved, m_block * n_block, MPI_DOUBLE, row_prev, 0, col_comm, &status);	
		}else {

			MPI_Recv(A_saved, m_block * n_block, MPI_DOUBLE, col_prev, 0, row_comm, &status);	
			MPI_Send(A_current, m_block * n_block, MPI_DOUBLE, col_next, 0, row_comm);
			
			MPI_Recv(B_saved, m_block * n_block, MPI_DOUBLE, row_prev, 0, col_comm, &status);	
			MPI_Send(B_current, m_block * n_block, MPI_DOUBLE, row_next, 0, col_comm);
		}
		*networkElapsedTime += (MPI_Wtime() - network_start_time);
        
        
        
		double processing_start_time = MPI_Wtime();
        multMatricesLineByLine(m_block, n_block, k_block, A_current, B_current, C_local);
    	*processingElapsedTime += (MPI_Wtime() - processing_start_time);

        
		broadcaster = ring_next(broadcaster, n_procs_row);
	}while(broadcaster != init_broadcaster);

}


void sizeMatrixesValidation(int n_procs_row, int n_procs_col,int m_block, int n_block, int k_block) {
	if (n_procs_col * n_procs_row != number_of_processes) {
        std::cerr << "number of proccessors must be a perfect square!" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

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
}


void compute(){
    //Assuming that the processes form a square
    int n_procs_row = sqrt(number_of_processes);
    int n_procs_col = n_procs_row;
    
    int m_block = m / n_procs_row;
    int n_block = n / n_procs_col;
    int k_block = k / n_procs_col;
    
    //create comm_groups
    int n_dims = 2;
    int dims[n_dims] = {n_procs_row, n_procs_col};
    int periods[n_dims] = {0, 0};
    int repeat = 0;
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, n_dims, dims, periods, repeat, &comm_cart);
	
	sizeMatrixesValidation(n_procs_row, n_procs_col, m_block,n_block,k_block);


    double * A_local = (double *) calloc(m_block * n_block, sizeof(double));
    double * B_local = (double *) calloc(n_block * k_block, sizeof(double));
    double * C_local = (double *) calloc(m_block * k_block, sizeof(double));
    double * A = (double *) calloc(m * n, sizeof(double));
    double * B = (double *) calloc(n * k, sizeof(double));
    double * C = (double *) calloc(m * n, sizeof(double));

    initMatrices(A_local, B_local, C_local, m_block, n_block, comm_cart);

    gatherMatrix(m_block, n_block, A_local, m, n, A);
    gatherMatrix(n_block, k_block, B_local, n, k, B);
    
       
    //testing the scatter function
/*    distributeSquareMatrix(A_1, n, C_1);

    if(rank == 0){
        double * A_1 = (double *) calloc(n * k, sizeof(double));
        double * B_1 = (double *) calloc(n * k, sizeof(double));
        double * C_1 = (double *) calloc(n * k, sizeof(double));

        initMatrices(A_1, B_1, C_1, n, n, comm_cart);

    }
*/


	/*if(rank == 0){
	std::cout << "MATRIX A: " << std::endl;
	printMatrix(m,n,A);
	std::cout << "MATRIX B: " << std::endl;
	printMatrix(m,n,A);
	}*/

	double networkElapsedTime; 
	double processingElapsedTime;
    double start_time = MPI_Wtime();
    summa_ring(comm_cart, m_block, n_block, k_block, A_local, B_local, C_local,&networkElapsedTime,&processingElapsedTime);
    double end_time = MPI_Wtime();
    getTimes(start_time, end_time, networkElapsedTime, processingElapsedTime);

	gatherMatrix(m_block, k_block, C_local, m, k, C);
	if(rank == 0){
		std::cout << "MATRIX C: " << std::endl;
		printMatrix(m,n,C);
	} 

	free(A_local);
	free(B_local);
	free(C_local);
	free(A);
	free(B);
	free(C);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
/*
	m = n = k = atoi(argv[1]);
	if(rank == 0) {
		std::cout << "======================" << std::endl;
		std::cout << "n: " << n << ", m: " << m << ", k: " << k << std::endl;
	}*/

	for(int nn = 1024; nn <= 9216; nn = 1024 + nn){
		for(int s = 0; s < 5; s++){
			if(rank == 0)std::cout << "Size: " << nn << ", iteration: " << s << std::endl; 
			m = n = k = nn;
			compute();
			
			MPI_Barrier(MPI_COMM_WORLD);
			sleep(10);
		}
	}


    MPI_Finalize();
}


