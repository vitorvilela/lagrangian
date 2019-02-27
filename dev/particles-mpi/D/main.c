#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "lag_map.h"
#include "handle_stoch_particles.h"



int main(int argc, char *argv[]) {
     
  int my_rank = 0, i = 0;
  
  printf("\n\n\t\t*** START *** PROC %d\n\n", my_rank);
  
  double simulation_time = 0.0; 
  int ct_seed = 0;
  char in[1];  
        
  const int size_B = npl*8+1;
  
  int *B;
  B = malloc(17*sizeof(int));
  
  B[0] = 2; B[1] = 0; B[2] = 1; B[3] = 1; B[4] = 1; B[5] = 16;
  B[6] = 16; B[7] = 16; B[8] = 1; B[9] = 1; B[10] = 1; B[11] = 21;
  B[12] = 1; B[13] = 32; B[14] = 12; B[15] = 32; B[16] = 2048;
 
  const double A[15] = {ga1, gb1, ga2, gb2, ga3, gb3, 6.25E-2, 3.125E-2, 6.25E-2, 3.125E-2, 6.25E-2, 3.125E-2};
    
  Map_Init(A, B, size_B); 
  U_Particles_Init(my_rank, ga1, ga2, ga3, gb1, gb2, gb3, A);
  Particles_Zeta_Init(1);
    
  while(simulation_time < FINAL_TIME) {
    
    i++;      
    simulation_time += dt;
    printf("\nsimulation_time = %lf\n", simulation_time);
    ct_seed = i;
    printf("\nct_seed = %d\n", ct_seed);

    Transport_Stoch_Particles(A, &ct_seed, my_rank);
   
    if(i % (int)(FINAL_TIME/(dt*nfiles)) == 0)
      Print_Particles_VTK(i-1, my_rank);
    
  }
    
 
  printf("\n\n\t\t*** END ***\n\n");
  
//   MPI_Finalize();
  
  return 0;
  
}
