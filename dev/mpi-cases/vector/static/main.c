#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{       

    int myrank, i;
    
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    int FLAG = 1;
    MPI_Status status; 
    MPI_Request R;
      
    const int VSIZE = 10;
    double vector[VSIZE];
    
    for(i=0; i<VSIZE; i++)
      vector[i] = -1.0;
    
    for(i=0; i<VSIZE; i++)
      printf("\nOuside before - myrank: %i / vector: %f\n", myrank, vector[i]);
    printf("\n");
    
    for(i=0; i<VSIZE; i++)
      vector[i] = -1.0*(-myrank);
    
    if(myrank == 0) {    
      for(i=0; i<VSIZE; i++)
	printf("\nInside - myrank: %i / vector: %f\n", myrank, vector[i]);
      printf("\n");
      MPI_Isend(&vector, VSIZE, MPI_DOUBLE, 1, FLAG, MPI_COMM_WORLD, &R);
      MPI_Wait(&R, &status);
    }
    
    if(myrank == 1) {  
      for(i=0; i<VSIZE; i++)
	printf("\nInside before - myrank: %i / vector: %f\n", myrank, vector[i]);
      printf("\n");
      MPI_Recv(&vector, VSIZE, MPI_DOUBLE, 0, FLAG, MPI_COMM_WORLD, &status);
      for(i=0; i<VSIZE; i++)
	printf("\nInside after - myrank: %i / vector: %f\n", myrank, vector[i]);
      printf("\n");      
    }
    
    for(i=0; i<VSIZE; i++)
      printf("\nOuside after - myrank: %i / vector: %f\n", myrank, vector[i]);
    printf("\n");
    
    MPI_Finalize();
    
    return 0;
  
}
