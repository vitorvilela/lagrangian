#ifndef map_H
#define map_h

#define eps 0.000001

typedef struct id_np {
  int id;
  UT_hash_handle hh;
} id_np_t;

typedef struct index {
  int igx, igy, igz;  
} index_t;

typedef struct volume {
  int rl;
  int iug;
  index_t index;
  struct volume *rvolume;
  struct id_np *id_np;  
  UT_hash_handle hh;
} volume_t;


volume_t *eulerian = NULL;


volume_t* SEARCH_TABLE(const int level, const double x, const double y, const double z, const double A[], const int levels) {
  
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


void INITIALIZE_MAP(const double A[], const int B[], const int size_B) {  
   
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
	  ptr_volume->id_np = NULL;
	  ptr_volume->rvolume = NULL;
	
	  if(rl == 0)    
	    HASH_ADD(hh, eulerian, index, sizeof(index_t), ptr_volume);
	     
	  else {  
	   
	    x = i*A[6+rl];
	    y = j*A[6+levels+rl];
	    z = k*A[6+2*levels+rl];        
	   
	    base_volume = SEARCH_TABLE(rl, x-eps, y-eps, z-eps, A, levels); 
	    
	    HASH_ADD(hh, base_volume->rvolume, index, sizeof(index_t), ptr_volume);
	  
	  }
	    
	}
      }
    }
        
  }
  
}


volume_t* SEARCH_PARTICLE(const double A[], const double x, const double y, const double z, const int levels) {
    
  index_t index;
  
  volume_t *ptr_volume, *tmp_volume;
  
  char in[1];      
  
  int l = 0;
      
  index.igx = ceil((x-A[0])/A[6]);
  index.igy = ceil((y-A[2])/A[6+levels]);
  index.igz = ceil((z-A[4])/A[6+2*levels]);   
    
  HASH_FIND(hh, eulerian, &index, sizeof(index_t), ptr_volume);
  
  if(ptr_volume == NULL) 
    printf("\n\n\t\t*** ERROR IN MAP - THERE IS NO VOLUME WITH THIS INDEX IN HASH TABLE ***\n\n");
   
  if(ptr_volume->rvolume == NULL) 
    return ptr_volume;
  
  else {
  
    do { 
          
      l++;
      
      index.igx = ceil((x-A[0])/A[6+l]);
      index.igy = ceil((y-A[2])/A[6+levels+l]);
      index.igz = ceil((z-A[4])/A[6+2*levels+l]); 
    
      HASH_FIND(hh, ptr_volume->rvolume, &index, sizeof(index_t), tmp_volume);
      
      if(tmp_volume == NULL) 
	printf("\n\n\t\t*** ERROR IN MAP - THERE IS NO VOLUME WITH THIS INDEX IN HASH TABLE ***\n\n");
     
      ptr_volume = tmp_volume;
    
    } 
    while(ptr_volume->rvolume != NULL);
  
    return ptr_volume;
    
  }

}


#endif

