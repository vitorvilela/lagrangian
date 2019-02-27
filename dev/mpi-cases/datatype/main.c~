#include "mpi.h"
#include <stdio.h>

#define N 10
#define ITENS 3


typedef struct Particle {
  int id;
  char kind;
  double p[N];
} Particle_t;


int main(int argc, char *argv[])
{      
  
  int myrank, i, j;
    
  int FLAG1 = 1, FLAG2 = 2;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
  Particle_t *ptr_particle;
    
  int nparticles;
  
  MPI_Status status;
  MPI_Request R;
  
  MPI_Datatype MPI_Particle_type;
  MPI_Datatype type[ITENS] = {MPI_INT, MPI_CHAR, MPI_DOUBLE};
  MPI_Aint disp[ITENS], address[ITENS];   
  
  int blocklen[ITENS] = {1, 1, N};
   
 
  if(myrank == 0) {
    
    nparticles = 2;
    
    printf("\nInside - myrank: %i - nparticles: %i\n", myrank, nparticles);
    
    MPI_Send(&nparticles, 1, MPI_INT, 1, FLAG1, MPI_COMM_WORLD);
    
    ptr_particle = malloc(nparticles*sizeof(Particle_t));

    MPI_Address(&ptr_particle[0].id, &address[0]);
    MPI_Address(&ptr_particle[0].kind, &address[1]);
    MPI_Address(&ptr_particle[0].p, &address[2]);
    
    disp[0] = address[0] - address[0];
    disp[1] = address[1] - address[0];
    disp[2] = address[2] - address[0];
    
    MPI_Type_create_struct(ITENS, blocklen, disp, type, &MPI_Particle_type);
    MPI_Type_commit(&MPI_Particle_type);    
    
    for(i=0; i<nparticles; i++) {
      ptr_particle[i].id = i+1;
      ptr_particle[i].kind = 'P';
      for(j=0; j<N; j++)
	ptr_particle[i].p[j] = 1.0*j;     
    }
    
    for(i=0; i<nparticles; i++) {
      printf("\nInside - myrank: %i - ptr_particle[%i].id: %i\n", myrank, i, ptr_particle[i].id);
      printf("\nInside - myrank: %i - ptr_particle[%i].kind: %c\n", myrank, i, ptr_particle[i].kind);
      for(j=0; j<N; j++)
	printf("\n\tInside - myrank: %i - ptr_particle[%i].p[%i]: %f\n", myrank, i, j, ptr_particle[i].p[j]);
    }
    printf("\n");    
    
    MPI_Isend(ptr_particle, nparticles, MPI_Particle_type, 1, FLAG2, MPI_COMM_WORLD, &R);
    MPI_Wait(&R, &status);   
    
    free(ptr_particle);
    ptr_particle = NULL;
    
  }
     
     
  if(myrank == 1) {  
      
    MPI_Recv(&nparticles, 1, MPI_INT, 0, FLAG1, MPI_COMM_WORLD, &status);
      
    printf("\nInside after - myrank: %i - nparticles: %i\n", myrank, nparticles);
      
    ptr_particle = malloc(nparticles*sizeof(Particle_t));
     
    MPI_Address(&ptr_particle[0].id, &address[0]);
    MPI_Address(&ptr_particle[0].kind, &address[1]);
    MPI_Address(&ptr_particle[0].p, &address[2]);
    
    disp[0] = address[0] - address[0];
    disp[1] = address[1] - address[0];
    disp[2] = address[2] - address[0];
    
    MPI_Type_create_struct(ITENS, blocklen, disp, type, &MPI_Particle_type);
    MPI_Type_commit(&MPI_Particle_type); 
    
    MPI_Recv(ptr_particle, nparticles, MPI_Particle_type, 0, FLAG2, MPI_COMM_WORLD, &status);
      
    for(i=0; i<nparticles; i++) {
      printf("\nInside after - myrank: %i - ptr_particle[%i].id: %i\n", myrank, i, ptr_particle[i].id);
      printf("\nInside after - myrank: %i - ptr_particle[%i].kind: %c\n", myrank, i, ptr_particle[i].kind);
      for(j=0; j<N; j++)
	printf("\n\tInside after - myrank: %i - ptr_particle[%i].p[%i]: %f\n", myrank, i, j, ptr_particle[i].p[j]);
    }
    printf("\n");  
      
    free(ptr_particle);
    ptr_particle = NULL;
      
  }
    
  MPI_Finalize();
    
  return 0;
  
}
