#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "uthash.h"
#include "lag_map.h"

#define eps 0.000001


volume_t *eulerian = NULL;


volume_t* Search_Table(const int level, const double x, const double y, const double z, const double A[], const int levels) {
  
  int l;  
  
  char in[1];
  
  index_t index;  
  
  volume_t *ptr_volume, *sub_volume;    
    
  index.igx = ceil((x-A[0])/A[6]);
  index.igy = ceil((y-A[2])/A[6+levels]);
  index.igz = ceil((z-A[4])/A[6+2*levels]);

  HASH_FIND(hh, eulerian, &index, sizeof(index_t), ptr_volume);
      
  if(ptr_volume == NULL) 
    printf("\n\n\t\t*** ERROR IN MAP - THERE IS NO VOLUME WITH THIS INDEX IN HASH TABLE ***\n\n");
   
  if(level > 1) {
    
    for(l=1; l<level; l++) { 
            
      index.igx = ceil((x-A[0])/A[6+l]);
      index.igy = ceil((y-A[2])/A[6+levels+l]);
      index.igz = ceil((z-A[4])/A[6+2*levels+l]); 
    
      HASH_FIND(hh, ptr_volume->rvolume, &index, sizeof(index_t), sub_volume);
      
      if(sub_volume == NULL) 
	printf("\n\n\t\t*** ERROR IN MAP - THERE IS NO VOLUME WITH THIS INDEX IN HASH TABLE ***\n\n");
     
      ptr_volume = sub_volume;
          
    } 
    
  }
      
  return ptr_volume; 
  
}


void Map_Init(const double A[], const int B[], const int size_B) {  
   
  char in[1];
  
  int i, j, k;
  double x, y, z;
  int rl, igx, igy, igz, nx, ny, nz, iug;
  volume_t *ptr_volume, *base_volume, *tmp_volume;
  const int levels = B[0];
  int p = 0;
     
  
  while(p < size_B-1) {
     
    rl = B[++p];    
    igx = B[++p];
    igy = B[++p];
    igz = B[++p]; 
    nx = B[++p];
    ny = B[++p];
    nz = B[++p]; 
    iug = B[++p];    
    
    for(i=igx; i<igx+nx; i++) { 
      for(j=igy; j<igy+ny; j++) { 
	for(k=igz; k<igz+nz; k++) { 
	     
	  ptr_volume = malloc(sizeof(volume_t)); 
	     
	  ptr_volume->index.igx = i;
	  ptr_volume->index.igy = j;
	  ptr_volume->index.igz = k;  
	  ptr_volume->rl = rl;
	  ptr_volume->iug = iug;
	  ptr_volume->id_particle = NULL;
	  ptr_volume->rvolume = NULL;
	
	  if(rl == 0)     
	    HASH_ADD(hh, eulerian, index, sizeof(index_t), ptr_volume);
	  
	  
	  else {  
	   
	    x = i*A[6+rl];
	    y = j*A[6+levels+rl];
	    z = k*A[6+2*levels+rl];
	   
	    base_volume = Search_Table(rl, x-eps, y-eps, z-eps, A, levels); 
	    
	    HASH_ADD(hh, base_volume->rvolume, index, sizeof(index_t), ptr_volume);
	      
	  }
	    
	}
      }
    }
        
  }
  
}


volume_t* Search_Particle(const double A[], const double x, const double y, const double z, int levels) {
    
  index_t index;
  
  volume_t *ptr_volume, *tmp_volume;
  
  char in[1];
    
  int l = 0;
  
  index.igx = (int) (((x-A[0])/A[6])+1);
  index.igy = (int) (((y-A[2])/A[6+(levels)])+1);
  index.igz = (int) (((z-A[4])/A[6+2*(levels)])+1);
  

  HASH_FIND(hh, eulerian, &index, sizeof(index_t), ptr_volume);

  if(ptr_volume == NULL) 
    printf("\n\n\t\t*** ERROR IN MAP - THERE IS NO VOLUME WITH THIS INDEX IN HASH TABLE ***\n\n");
   
  if(ptr_volume->rvolume == NULL) 
    return ptr_volume; 
  
  else {
  
    do { 
          
      l++;
      
      index.igx = (int) (((x-A[0])/A[6+l])+1);
      index.igy = (int) (((y-A[2])/A[6+(levels)+l])+1);
      index.igz = (int) (((z-A[4])/A[6+2*(levels)+l])+1); 

      HASH_FIND(hh, ptr_volume->rvolume, &index, sizeof(index_t), tmp_volume);
      
      if(tmp_volume == NULL) 
	printf("\n\n\t\t*** ERROR IN MAP - THERE IS NO VOLUME WITH THIS INDEX IN HASH TABLE ***\n\n");
     
      ptr_volume = tmp_volume;
    
    } 
    
    while(ptr_volume->rvolume != NULL);
  
    return ptr_volume;
    
  }

}


