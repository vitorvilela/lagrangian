#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include "uthash.h"
#include "lag_map.h"
#include "handle_stoch_particles.h"


//Global pointers

notionalParticle_t *lagrangian = NULL;
idStack_t *idstack;


//Stack's functions

idStack_t* createStack(int max_np_proc) {
  
  idStack_t *S;
  S = malloc(sizeof(idStack_t));
  S->id = malloc(max_np_proc*sizeof(int));
  S->size = 0;
  S->capacity = max_np_proc;
  
  return S;
  
}

void pop(idStack_t *S) {
  
  if(S->size == 0)
    printf("Stack is Empty");

  else     
    S->size--;
   
}

int top(idStack_t *S) {
  
  if(S->size == 0)
    printf("Stack is Empty");

  else     
    return S->id[S->size-1];    
  
  
}

void push(idStack_t *S, int id) {
  
  if(S->size == S->capacity) 
    printf("Stack is Full");
  
  else 
    S->id[S->size++] = id;
   
}

void Stack_Init(const int max_np_proc) {
  
  int i;
  
  idstack = createStack(max_np_proc);
  for(i=max_np_proc; i>=1; i--)
    push(idstack, i);
  
}




//Transport's functions

void Manual_Particle_Init(const double A[], const int id, const double Position[]) {
    
  volume_t *ptr_volume;
  id_lagrangian_t *ptr_id_p;
  notionalParticle_t *ptr_particle;
  
  ptr_particle = malloc(sizeof(notionalParticle_t));        
  ptr_particle->id = id;
  ptr_particle->posixo = Position[0];
  ptr_particle->posiyo = Position[1];
  ptr_particle->posizo = Position[2];
  ptr_particle->posix = Position[0];
  ptr_particle->posiy = Position[1];
  ptr_particle->posiz = Position[2];  
  HASH_ADD(hh, lagrangian, id, sizeof(int), ptr_particle);
      
  pop(idstack);
  
  ptr_volume = Search_Particle(A, ptr_particle->posixo, ptr_particle->posiyo, ptr_particle->posizo, npl);
  ptr_id_p = malloc(sizeof(id_lagrangian_t)); 
  ptr_id_p->id = id;      
  HASH_ADD_INT(ptr_volume->id_particle, id, ptr_id_p);
  
}



void U_Particles_Init(const int my_rank, const double a1, const double a2, const double a3, const double b1, const double b2, const double b3, const double A[]) {

  //Function's variables
  int max_np_proc = np_total;  
  long int seed = -3213853;
  long int *ptr_seed = &seed;
  int c, p, stack_id;
  volume_t *ptr_volume;
  id_lagrangian_t *ptr_id_p;
  notionalParticle_t *ptr_particle;
    
  Stack_Init(max_np_proc);
        
  const double lx = b1-a1;
  const double ly = b2-a2;
  const double lz = b3-a3;  
  const double vol_proc = lx*ly*lz;    
  const double vol_total = (A[1]-A[0])*(A[3]-A[2])*(A[5]-A[4]);
  const int np_proc = (int)(np_total*(vol_proc/vol_total));
      
  for(p=1; p<=np_proc; p++) {
  
    ptr_particle = malloc(sizeof(notionalParticle_t));        
    ptr_particle->id = p;
    HASH_ADD(hh, lagrangian, id, sizeof(int), ptr_particle);
  
    //Remove id from idstack **************
    stack_id = top(idstack);
    pop(idstack);
    //*************************************
    
    ptr_particle->proc = my_rank;
    ptr_particle->new_proc = my_rank;
            
    ptr_particle->posixo = a1 + RanTwo(ptr_seed)*lx;
    ptr_particle->posiyo = a2 + RanTwo(ptr_seed)*ly;
    ptr_particle->posizo = a3 + RanTwo(ptr_seed)*lz; 
    ptr_particle->posix = ptr_particle->posixo;
    ptr_particle->posiy = ptr_particle->posiyo;
    ptr_particle->posiz = ptr_particle->posizo;
    ptr_particle->index.igx = ptr_volume->index.igx;
    ptr_particle->index.igy = ptr_volume->index.igy;
    ptr_particle->index.igz = ptr_volume->index.igz;       
    ptr_particle->new_index.igx = ptr_volume->index.igx;
    ptr_particle->new_index.igy = ptr_volume->index.igy;
    ptr_particle->new_index.igz = ptr_volume->index.igz;
    ptr_particle->treact = 0.0;
    ptr_particle->dwall = 0.0;
    ptr_particle->omegap = 0.0;
    ptr_particle->rhop = 0.0;
    ptr_particle->weight = 0.0; 
    for(c=0; c<n_Yi; c++)
      ptr_particle->zetap[c] = 0.0;
     
    ptr_volume = Search_Particle(A, ptr_particle->posixo, ptr_particle->posiyo, ptr_particle->posizo, npl);
   
    ptr_id_p = malloc(sizeof(id_lagrangian_t)); 
    ptr_id_p->id = p;      
    HASH_ADD_INT(ptr_volume->id_particle, id, ptr_id_p); 
     
  }
  
}


void Particles_Zeta_Init(const int k) {
    
  int my_rank = 0;
  
  notionalParticle_t *iter_particle, *tmp_particle;
  
  int c;
  
  if(k == 0) {
    
    HASH_ITER(hh, lagrangian, iter_particle, tmp_particle) { 
      for(c=0; c<n_Yi; c++) {
	iter_particle->zetap[c] = 0.0;
	printf("\nP %d - Zeta %f \n", iter_particle->id, iter_particle->zetap[c]);
      }
    }
  }
  
  else if(k == 1) {
    
    HASH_ITER(hh, lagrangian, iter_particle, tmp_particle) { 
      for(c=0; c<n_Yi; c++) {
	iter_particle->zetap[c] = 0.5*(zeta1_mix+zeta2_mix) + 0.5*(zeta1_mix-zeta2_mix) * tanh((2.0*(iter_particle->posiyo))/mixlgth - 2.0*targ_mix/mixlgth);
	printf("\nP %d (%i) - Zeta %f \n",iter_particle->id, my_rank, iter_particle->zetap[c]);
      }
    }
  }
  
}

void Move_XYZ(notionalParticle_t * const iter_particle, int *ptr_seed) {
      
  double gasX = Gasdev(ptr_seed);
  double gasY = Gasdev(ptr_seed);
  double gasZ = Gasdev(ptr_seed);
  double factor = sqrt((dt)*(2.0*Ni/Sc));
      
  iter_particle->posix = iter_particle->posixo + U*dt + factor*gasX;
  iter_particle->posiy = iter_particle->posiyo + V*dt + factor*gasY;      
  iter_particle->posiz = iter_particle->posizo + W*dt + factor*gasZ;
      
}


void Bndy_Box(const double A[], notionalParticle_t * const iter_particle) {
      
  if(iter_particle->posix < A[0])
    iter_particle->posix = A[0] + (A[0] - iter_particle->posix);
   
  if(iter_particle->posix > A[1]) 
    iter_particle->posix = A[1] - (iter_particle->posix - A[1]);
    
  if(iter_particle->posiy < A[2])
    iter_particle->posiy = A[2] + (A[2] - iter_particle->posiy);

  if(iter_particle->posiy > A[3])
    iter_particle->posiy = A[3] - (iter_particle->posiy - A[3]);
      
  if(iter_particle->posiz < A[4])
    iter_particle->posiz = A[4] + (A[4] - iter_particle->posiz);

  if(iter_particle->posiz > A[5])
    iter_particle->posiz = A[5] - (iter_particle->posiz - A[5]);     
      
}

void Update_New_Index(const double A[], notionalParticle_t * const iter_particle) {

  volume_t *ptr_volume;
  
  ptr_volume = Search_Particle(A, iter_particle->posix, iter_particle->posiy, iter_particle->posiz, npl);
  iter_particle->new_index.igx = ptr_volume->index.igx;
  iter_particle->new_index.igy = ptr_volume->index.igy;
  iter_particle->new_index.igz = ptr_volume->index.igz;
  
  ptr_volume = NULL;

}


void Map_ID(const double A[], notionalParticle_t * const iter_particle) {

  int id;

  volume_t *ptr_volume, *ptr_old_vol;

  id_lagrangian_t *add_id, *rm_id;

  if( (iter_particle->new_index.igx != iter_particle->index.igx) || (iter_particle->new_index.igy != iter_particle->index.igy) || (iter_particle->new_index.igz != iter_particle->index.igz) ) {
  
    ptr_volume = Search_Particle(A, iter_particle->posix, iter_particle->posiy, iter_particle->posiz, npl);
    id = iter_particle->id;
    add_id = malloc(sizeof(id_lagrangian_t)); 
    add_id->id = id;
    HASH_ADD_INT(ptr_volume->id_particle, id, add_id);
      
    ptr_old_vol = Search_Particle(A, iter_particle->posixo, iter_particle->posiyo, iter_particle->posizo, npl);
    HASH_FIND_INT(ptr_old_vol->id_particle, &id, rm_id);
    HASH_DEL(ptr_old_vol->id_particle, rm_id);
    free(rm_id);
    rm_id = NULL; 
  
    ptr_volume = NULL;
    ptr_old_vol = NULL;
    
  }
    
}


void Update_Particle_Position(notionalParticle_t * const iter_particle) {
      
  iter_particle->posixo = iter_particle->posix;
  iter_particle->posiyo = iter_particle->posiy;
  iter_particle->posizo = iter_particle->posiz;
    
  iter_particle->index.igx = iter_particle->new_index.igx;
  iter_particle->index.igy = iter_particle->new_index.igy;
  iter_particle->index.igz = iter_particle->new_index.igz;
  
}

void Transport_Stoch_Particles(const double A[], int *ct_seed, const int my_rank) {
       
  //Function's variables  
  int seed = -3213853;
  int *ptr_seed;
  notionalParticle_t *iter_particle, *tmp_particle;

  //Seed assignment 
  if(*ct_seed == 1)
    ptr_seed = &seed;
  else
    ptr_seed = ct_seed;
    
  //Time's progress  
  HASH_ITER(hh, lagrangian, iter_particle, tmp_particle) {
         
    Move_XYZ(iter_particle, ptr_seed);
    Bndy_Box(A, iter_particle);
    Update_New_Index(A, iter_particle);
    Map_ID(A, iter_particle);
    Update_Particle_Position(iter_particle);  

  }
  
}



//Random functions

float RanOne(int *idum) {
  
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0 || !iy) { 								// Initialize.
  if (-(*idum) < 1) *idum=1; 								// Be sure to prevent idum = 0.
  else *idum = -(*idum);
  for (j=NTAB+7;j>=0;j--) { 								// Load the shuffle table (after 8 warm-ups).
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  if (j < NTAB) iv[j] = *idum;
  }
  iy=iv[0];
  }
  k=(*idum)/IQ; 									// Start here when not initializing.
  *idum=IA*(*idum-k*IQ)-IR*k; 								// Compute idum=(IA*idum) % IM without overif
  if(*idum < 0) *idum += IM; 								// flows by Schrage’s method.
  j=iy/NDIV_RanOne; 										// Will be in the range 0..NTAB-1.
  iy=iv[j]; 										// Outransfut previously stored value and refill the
  iv[j] = *idum; 									// shuffle table.
  if ((temp=AM_RanOne*iy) > RNMX) return RNMX; 						// Because users don’t expect endpoint values.
  else return temp;
    
}


float RanTwo(long int *idum) {

  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0) { 									//Initialize.
    if (-(*idum) < 1) *idum=1;								//Be sure to prevent idum = 0.
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { 								//Load the shuffle table (after 8 warm-ups).
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1; 									//Start here when not initializing.
  *idum=IA1*(*idum-k*IQ1)-k*IR1; 							//Compute idum=(IA1*idum) % IM1 without
  if (*idum < 0) *idum += IM1; 								//overflows by Schrage’s method.
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; 							//Compute idum2=(IA2*idum) % IM2 likewise.
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV; 										//Will be in the range 0..NTAB-1.
  iy=iv[j]-idum2; 									//Here idum is shuffled, idum and idum2 are
  iv[j] = *idum; 									//combined to generate outransfut.
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;							// Because users don’t expect endpoint values.
  else return temp;

}  


float Gasdev(int *idum) {

  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  if (*idum < 0) iset=0; 
  if (iset == 0) {
    do {
      v1=2.0*RanOne(idum)-1.0; 
      v2=2.0*RanOne(idum)-1.0; 
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0); 
    fac=sqrt(-2.0*log(rsq)/rsq);
 gset=v1*fac;
    iset=1;
    return v2*fac;
  } 
  else { 
    iset=0; 
    return gset;    
  } 
  
}



//Miscelaneous functions

void Print_Particles_VTK(const int ct, const int my_rank) {
    
  char string[40];
  
  int part_counter = 0;
  
  notionalParticle_t *iter_particle, *tmp_particle;
  
  part_counter = HASH_COUNT(lagrangian);
  
  sprintf(string, "particles_proc%i_ct%i.vtk", my_rank, ct);
    
  FILE *pvtk;
  pvtk = fopen(string, "w");
  fprintf(pvtk, "# vtk DataFile Version 2.0\nPoints\nASCII\nDATASET POLYDATA\n");
  fprintf(pvtk, "Points %i float\n", part_counter);
       
  HASH_ITER(hh, lagrangian, iter_particle, tmp_particle) {
    fprintf(pvtk, "%1.2e %1.2e %1.2e \n", iter_particle->posixo, iter_particle->posiyo, iter_particle->posizo);
  }
  
  fprintf(pvtk, "POINT_DATA %d\n", part_counter);
  fprintf(pvtk, "SCALARS zeta float\n");
  fprintf(pvtk, "LOOKUP_TABLE default\n");
    
  HASH_ITER(hh, lagrangian, iter_particle, tmp_particle) {
    fprintf(pvtk, "%f \n", iter_particle->zetap[0]);
  }
  
  fclose(pvtk);
     
}


int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1) {
    
  long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
  result->tv_sec = diff / 1000000;
  result->tv_usec = diff % 1000000;

  return (diff<0);
  
}

