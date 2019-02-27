#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "lag_map.h"
#include "handle_stoch_particles.h"
#include "mpi.h"



int main(int argc, char *argv[]) {
     
  int my_rank, i=0;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  
  printf("\n\n\t\t*** START *** PROC %i\n\n", my_rank);

  
  double stime = 0.0; 
  int ct_seed;
  char in[1];  
    
  int seed = -3213853;
  int *ptr_seed = &seed;    

  
  struct timeval tvBegin, tvEnd, tvDiff;   
  gettimeofday(&tvBegin, NULL);
      
  double a1=0, b1=0, a2=0, b2=0, a3=0, b3=0;
  const int size_B = npl*8+1;
  
  int *B;
  B = malloc(17*sizeof(int));
  
  if(my_rank == 0) {    
    B[0] = 2; B[1] = 0; B[2] = 1; B[3] = 1; B[4] = 1; B[5] = 8;
    B[6] = 16; B[7] = 16; B[8] = 1; B[9] = 1; B[10] = 1; B[11] = 21;
    B[12] = 1; B[13] = 16; B[14] = 12; B[15] = 32; B[16] = 2048;
    a1 = 0.0;
    b1 = 0.5;
    a2 = 0.0;
    b2 = 1.0;
    a3 = 0.0;
    b3 = 1.0;
  }
  
  if(my_rank == 1) {
    B[0] = 2; B[1] = 0; B[2] = 9; B[3] = 1; B[4] = 1; B[5] = 8; B[6] = 16;
    B[7] = 16; B[8] = 1; B[9] = 1; B[10] = 17; B[11] = 21; B[12] = 1;
    B[13] = 16; B[14] = 12; B[15] = 32; B[16] = 2048;
    a1 = 0.5;
    b1 = 1.0;
    a2 = 0.0;
    b2 = 1.0;
    a3 = 0.0;
    b3 = 1.0;    
  }
  
  const double A[15] = {ga1, gb1, ga2, gb2, ga3, gb3, 6.25E-2, 3.125E-2, 6.25E-2, 3.125E-2, 6.25E-2, 3.125E-2};
  
  Map_Init(A, B, size_B);
  printf("\nMap_Init() - PROC %i\n", my_rank);
  
  U_Particles_Init(my_rank, a1, a2, a3, b1, b2, b3, A);
  printf("\n U_Particles_Init() - PROC %i\n", my_rank);
  
  Particles_Zeta_Init(1);
  printf("\n Particles_Zeta_Init() - PROC %i\n", my_rank);

  
  //*****************

  MPI_Datatype MPI_Particle_type;

  Create_MPI_Particle_Type(&MPI_Particle_type);

  
  //**************
  
  
  
  
  
  while(stime < FINAL_TIME) {
      
      i++;      
      stime += dt;
      
      ct_seed = i;
      Transport_Stoch_Particles(A, a1, a2, a3, b1, b2, b3, &ct_seed, my_rank);

      if(i % (int)(FINAL_TIME/(dt*nfiles)) == 0) {
	Print_Particles_VTK(&i, &my_rank);

      }
      
    }
    

  
  gettimeofday(&tvEnd, NULL);  
  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
  printf("\n\n\t\tTime, in seconds, for PaSR with NPV (%i) = ", np_total);
  printf("%ld.%06ld\n\n", tvDiff.tv_sec, tvDiff.tv_usec);
 
  printf("\n\n\t\t*** END ***\n\n");
  
  MPI_Finalize();
  
  return 0;
  
}
