#include <stdio.h>
#include <math.h>
#include <mpi.h>

int main(int argc, char **argv) {
    MPI::Init(argc, argv);
    int number_of_processes = MPI::COMM_WORLD.Get_size();
    int rank = MPI::COMM_WORLD.Get_rank();
    MPI::Group world_group_id = MPI::COMM_WORLD.Get_group();

    //create comm_groups
    

    int block_size = n / sqrt(rank);
    int block[block_size*block_size];
    int buffer[block_size];

    for(int k = 0; k < n-1; k++){

      //Broadcasting line
      if(/* a part of the column is owned by me do*/){
        for(int i = 0; i < processes_by_line - 1; i++){
          buffer = {line_piece_1, line_piece_2, "..."}; //FROM A
          MPI_Isend(&buffer, buffer.size(), MPI_INT, destination, tag, MPI_COMM_WORLD, &request);
          //OU MPI_Ibcast()
        }
      }

      if(/* a part of the line is owned by me do*/){
        for(int j = 0; j < processes_by_column - 1; j++){
          buffer = {column_piece_1, column_piece_2, "..."}; //FROM B
          MPI_Isend(&buffer, buffer.size(), MPI_INT, destination, tag, MPI_COMM_WORLD, &request);
        }
      }

      int block_column_recv[block_size];
      int block_line_recv[block_size];
      //receive from processors of the same row or same column
      MPI_Irecv(&block_column_recv, buffer.size(), MPI_INT, source, MPI_ANY_TAG /*??*/, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);

      MPI_Irecv(&block_line_recv, buffer.size(), MPI_INT, source, MPI_ANY_TAG /*??*/, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);




    }

}
