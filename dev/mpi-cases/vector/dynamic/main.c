#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{       

    int myrank, i;
    
    MPI_Init(&argc, &argv);    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    int FLAG1 = 1, FLAG2 = 2;
    
    MPI_Status status; 
    MPI_Request R;
      
    int vsize;
    
    double *vector;
        
    if(myrank == 0) {   
      
      vsize = 10;
      
      MPI_Send(&vsize, 1, MPI_INT, 1, FLAG1, MPI_COMM_WORLD);
      
      vector = malloc(vsize*(sizeof(double)));
      
      for(i=0; i<vsize; i++)
	vector[i] = 1.0*vsize;
      
      for(i=0; i<vsize; i++)
	printf("\nInside - myrank: %i - vector[%i]: %f\n", myrank, i, vector[i]);
      printf("\n");
      
      MPI_Isend(vector, vsize, MPI_DOUBLE, 1, FLAG2, MPI_COMM_WORLD, &R);
      MPI_Wait(&R, &status);
      
      free(vector);
      vector = NULL;
      
    }
    
    if(myrank == 1) {  
      
      MPI_Recv(&vsize, 1, MPI_INT, 0, FLAG1, MPI_COMM_WORLD, &status);
      
      printf("\nInside after - myrank: %i - vsize: %i\n", myrank, vsize);
      
      vector = malloc(vsize*(sizeof(double)));
      
      MPI_Recv(vector, vsize, MPI_DOUBLE, 0, FLAG2, MPI_COMM_WORLD, &status);
      
      for(i=0; i<vsize; i++)
	printf("\nInside after - myrank: %i - vector[%i]: %f\n", myrank, i, vector[i]);
      printf("\n");  
      
      free(vector);
      vector = NULL;
      
    }
    
    MPI_Finalize();
    
    return 0;
  
}
