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
//   ptr_volume = malloc(sizeof(volume_t));
    
  index.igx = ceil((x-A[0])/A[6]);
    printf("\n index.igx = %i\n", index.igx);
  index.igy = ceil((y-A[2])/A[6+levels]);
  printf("\n index.igy = %i\n", index.igy);
  index.igz = ceil((z-A[4])/A[6+2*levels]);
    printf("\n index.igz = %i\n", index.igz);

    gets(in);
    
  HASH_FIND(hh, eulerian, &index, sizeof(index_t), ptr_volume);
  
  printf("\n HASH_FIND() %i %i %i - %i \n", ptr_volume->index.igx, ptr_volume->index.igy, ptr_volume->index.igz, ptr_volume->rl);

    
  if(ptr_volume == NULL) 
    printf("\n\n\t\t*** ERROR IN MAP - THERE IS NO VOLUME WITH THIS INDEX IN HASH TABLE ***\n\n");
   
  if(level > 1) {
    
    printf("\n level > 1\n");
    
    for(l=1; l<level; l++) { 
      
//       sub_volume = malloc(sizeof(volume_t));
      
      index.igx = ceil((x-A[0])/A[6+l]);
       printf("\n index.igx = %i\n", index.igx);
      index.igy = ceil((y-A[2])/A[6+levels+l]);
       printf("\n index.igy = %i\n", index.igy);
      index.igz = ceil((z-A[4])/A[6+2*levels+l]); 
       printf("\n index.igz = %i\n", index.igz);
    
      HASH_FIND(hh, ptr_volume->rvolume, &index, sizeof(index_t), sub_volume);
      
       printf("\nsub_volume %i %i %i - %i\n", sub_volume->index.igx, sub_volume->index.igy, sub_volume->index.igz, sub_volume->rl);
	  gets(in);
      
      if(sub_volume == NULL) 
	printf("\n\n\t\t*** ERROR IN MAP - THERE IS NO VOLUME WITH THIS INDEX IN HASH TABLE ***\n\n");
     
      ptr_volume = sub_volume;
      
//       sub_volume = NULL;    
//       free(sub_volume);
    
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
	  
	  printf("\nptr_volume %i %i %i - %i\n", ptr_volume->index.igx, ptr_volume->index.igy, ptr_volume->index.igz, ptr_volume->rl);
	  gets(in);  
	  
	  
	
	  if(rl == 0)    {
	    HASH_ADD(hh, eulerian, index, sizeof(index_t), ptr_volume);
	    printf("\nADD MAP ptr_volume %i %i %i - %i\n", ptr_volume->index.igx, ptr_volume->index.igy, ptr_volume->index.igz, ptr_volume->rl);
	    gets(in);    
	  }
	     
	  else {  
	   
	    x = i*A[6+rl];
	    y = j*A[6+levels+rl];
	    z = k*A[6+2*levels+rl];        
	   
	    base_volume = SEARCH_TABLE(rl, x-eps, y-eps, z-eps, A, levels); 
	    
	    HASH_ADD(hh, base_volume->rvolume, index, sizeof(index_t), ptr_volume);
	    printf("\nADD MAP ptr_volume %i %i %i - %i\n", ptr_volume->index.igx, ptr_volume->index.igy, ptr_volume->index.igz, ptr_volume->rl);
	    gets(in); 
	  
	  }
	  
// 	  ptr_volume = NULL;    
// 	  free(ptr_volume); 
	  
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
      
//   ptr_volume = malloc(sizeof(volume_t));
      
  index.igx = ceil((x-A[0])/A[6]);
         printf("\n index.igx = %i\n", index.igx);
  index.igy = ceil((y-A[2])/A[6+levels]);
         printf("\n index.igy = %i\n", index.igy);
  index.igz = ceil((z-A[4])/A[6+2*levels]);   
         printf("\n index.igz = %i\n", index.igz);

    
  HASH_FIND(hh, eulerian, &index, sizeof(index_t), ptr_volume);
  printf("\nHASH_FIND() ptr_volume %i %i %i\n", ptr_volume->index.igx, ptr_volume->index.igy, ptr_volume->index.igz);
  gets(in); 
  
  if(ptr_volume == NULL) 
    printf("\n\n\t\t*** ERROR IN MAP - THERE IS NO VOLUME WITH THIS INDEX IN HASH TABLE ***\n\n");
   
  if(ptr_volume->rvolume == NULL) 
    return ptr_volume;
  
  else {
  
    printf("\nif(ptr_volume->rvolume != NULL)\n");
    gets(in); 
  
    do { 
          
      l++;
     
//       sub_volume = malloc(sizeof(volume_t));
      
      index.igx = ceil((x-A[0])/A[6+l]);
             printf("\n index.igx = %i\n", index.igx);
      index.igy = ceil((y-A[2])/A[6+levels+l]);
             printf("\n index.igy = %i\n", index.igy);
      index.igz = ceil((z-A[4])/A[6+2*levels+l]); 
             printf("\n index.igz = %i\n", index.igz);

    gets(in);
    
      HASH_FIND(hh, ptr_volume->rvolume, &index, sizeof(index_t), tmp_volume);
      
      if(tmp_volume == NULL) 
	printf("\n\n\t\t*** ERROR IN MAP - THERE IS NO VOLUME WITH THIS INDEX IN HASH TABLE ***\n\n");
     
      ptr_volume = tmp_volume;
      
//       sub_volume = NULL;    
//       free(sub_volume);
    
    } 
    while(ptr_volume->rvolume != NULL);
  
    return ptr_volume;
    
  }
  
  printf("\nleaving find particle\n");
  gets(in); 
  
}


#endif

