#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{       

    int myrank;
    
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
    int scalar = -1;
    
    int FLAG = 1;
    MPI_Status status;
    
    MPI_Request R;

    printf("\nOuside before - myrank: %i / scalar: %i\n", myrank, scalar);
 
    scalar = myrank;
    
    if(myrank == 0) {
      printf("\nInside - myrank: %i / scalar: %i\n", myrank, scalar);
      MPI_Isend(&scalar, 1, MPI_INT, 1, FLAG, MPI_COMM_WORLD, &R);
      MPI_Wait(&R, &status);
    }
    
    if(myrank == 1) {
      printf("\nInside before - myrank: %i / scalar: %i\n", myrank, scalar);
      MPI_Recv(&scalar, 1, MPI_INT, 0, FLAG, MPI_COMM_WORLD, &status);
      printf("\nInside after - myrank: %i / scalar: %i\n", myrank, scalar);
    }
    
    printf("\nOuside after - myrank: %i / scalar: %i\n", myrank, scalar);

    MPI_Finalize();
    
    return 0;
  
}
