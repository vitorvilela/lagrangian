#include <stdio.h>	/* printf */
#include <stdlib.h>	/* malloc */
#include <time.h>	/* time */
#include <sys/time.h>	/* time */
#include "uthash.h"

#define NX 2
#define NPV 20000000

#define max(a, b) ((a)>(b)?(a):(b))
#define min(a, b) ((a)<(b)?(a):(b))

#define bool int
#define true 1
#define false 0





typedef struct {
  int ix;
} index_t;

typedef struct notional_particle {
  UT_hash_handle hh;   
  index_t index;
  struct notional_particle *sub;  
  unsigned int id;  
  unsigned int counter_id; //just for head particles - account number of particles that passed through the volume
  double pdf;
} notional_particle_t;


void INITIALIZE_HEADS();

void ADD_PARTICLE(notional_particle_t *, unsigned int, double); 

void COPY_PARTICLE(notional_particle_t *, notional_particle_t *, unsigned int);

int id_sort(notional_particle_t *, notional_particle_t *);

int timeval_subtract(struct timeval *, struct timeval *, struct timeval *);



notional_particle_t *particles = NULL;

int main(int argc, char *argvp[]) {
    
  int i, p;
  
  unsigned up;
  
  index_t index;
  
  struct timeval tvBegin, tvEnd, tvDiff;
   
  char in[1];
    
  int volume_counter = 0;
      
  notional_particle_t *ptr_particle, *ptr_head, *toHead, *iter_head, *iter_particle, *tmp1, *tmp2;
    
//  gettimeofday(&tvBegin, NULL);
  
  INITIALIZE_HEADS();
  
//  gettimeofday(&tvEnd, NULL);  
//  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
//  printf("\n\n\t\tTime, in seconds, for INITIALIZE_HEADS with NX (%i) = ", NX);
//  printf("%ld.%06ld\n\n", tvDiff.tv_sec, tvDiff.tv_usec);
   
//  gets(in);
    
//  gettimeofday(&tvBegin, NULL);
  
  ptr_head = malloc(sizeof(notional_particle_t));   
  index.ix = 0;  
  HASH_FIND(hh, particles, &index, sizeof(index_t), ptr_head);
  if(ptr_head == NULL) 
    printf("\n\n\t\t*** There is no Head with this Index ***\n\n");
  else
     //printf("\n\n\t\t*** GET_HEAD(%i) ***\n\n", ptr_head->index.ix);
   
  for(p=1; p<=NPV; p++) {
    ADD_PARTICLE(ptr_head, volume_counter+1, 10.0);
    volume_counter += 1;  
  }
      
  ptr_head = NULL;
  free(ptr_head);
 
//  gettimeofday(&tvEnd, NULL);  
//  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
//  printf("\n\n\t\tTime, in seconds, for ADD_PARTICLE with NPV (%i) and with NX (%i) = ", NPV, NX);
//  printf("%ld.%03ld\n\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
//  gets(in);
  
//  gettimeofday(&tvBegin, NULL);
  
//  printf("\n\n\t\t*** SHOW PARTICLES IN THE BOX ***\n\n");
//  HASH_ITER(hh, particles, iter_head, tmp1) {
//      printf("\n\t\t    Head (%i) \n", iter_head->index.ix);
//    HASH_ITER(hh, iter_head->sub, iter_particle, tmp2) {
//      printf("\n\t\t\t    Particle (%i) - PDF(%f) \n\n", iter_particle->id, iter_particle->pdf);
//    }
//  }
      
//  gettimeofday(&tvEnd, NULL);  
//  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
//  printf("\n\n\t\tTime, in seconds, for HASH_ITER with NPV (%i) and with NX (%i) = ", NPV, NX);
//  printf("%ld.%03ld\n\n", tvDiff.tv_sec, tvDiff.tv_usec);    
      
  
//  gets(in);
  
      
//  gettimeofday(&tvBegin, NULL);
  
  gettimeofday(&tvBegin, NULL);
  
  int pp;
  
  for(pp=1; pp<=NPV; pp++) {
  
    ptr_head = malloc(sizeof(notional_particle_t));  
    index.ix = 0;
    HASH_FIND(hh, particles, &index, sizeof(index_t), ptr_head);
    if(ptr_head == NULL) 
      printf("\n\n\t\t*** There is no Head with this Index ***\n\n");
    else
     // printf("\n\n\t\t*** GET_HEAD(%i) ***\n\n", ptr_head->index.ix);
    
 // gets(in);
     
  
  
    ptr_particle = malloc(sizeof(notional_particle_t));  
    p = pp;  
    HASH_FIND_INT(ptr_head->sub, &p, ptr_particle);  
    if(ptr_particle == NULL) 
      printf("\n\n\t\t*** There is no Particle with this id ***\n\n");
    else
     // printf("\n\n\t\t*** GET_PARTICLE(%i, %i) ***\n\n", ptr_head->index.ix, ptr_particle->id);
  
  
//  gettimeofday(&tvEnd, NULL);  
//  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
//  printf("\n\n\t\tTime, in seconds, for GET_PARTICLE with NPV (%i) and with NX (%i) = ", NPV, NX);
//  printf("%ld.%06ld\n\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
  //gets(in);
        
      toHead = malloc(sizeof(notional_particle_t));   
      index.ix = 1;
      HASH_FIND(hh, particles, &index, sizeof(index_t), toHead);
      if(toHead == NULL) 
	printf("\n\n\t\t*** There is no Head with this Index ***\n\n");
      else
	//printf("\n\n\t\t*** GET_HEAD(%i) ***\n\n", toHead->index.ix);
  
   // gets(in);
  
    //HASH_SORT(toHead->sub, id_sort);
    //p = HASH_COUNT(toHead->sub) + 1;
    up = toHead->counter_id + 1;
    
    COPY_PARTICLE(toHead, ptr_particle, p);
  
 
  
  //gets(in);
  
//  printf("\n\n\t\t*** SHOW PARTICLES IN THE BOX ***\n\n");
//  HASH_ITER(hh, particles, iter_head, tmp1) {
//      printf("\n\t\t    Head (%i) \n", iter_head->index.ix);
//    HASH_ITER(hh, iter_head->sub, iter_particle, tmp2) {
//      printf("\n\t\t\t    Particle (%i) - PDF(%f) \n\n", iter_particle->id, iter_particle->pdf);
//    }
//  }
  
//  gets(in);  
  
    HASH_DEL(ptr_head->sub, ptr_particle);
    //printf("\n\n\t\t*** HASH_DEL(%i, %i) ***\n\n", ptr_head->index.ix, ptr_particle->id);
    ptr_particle = NULL;
    free(ptr_particle); 
      
    ptr_head = NULL;
    free(ptr_head);
  
    toHead = NULL;
    free(toHead);
    
    //gets(in); 
    
 //   printf("\n\n\t\t*** SHOW PARTICLES IN THE BOX ***\n\n");
 //   HASH_ITER(hh, particles, iter_head, tmp1) {
 //     printf("\n\t\t    Head (%i) \n", iter_head->index.ix);
 //     HASH_ITER(hh, iter_head->sub, iter_particle, tmp2) {
//	printf("\n\t\t\t    Particle (%i) - PDF(%f) \n\n", iter_particle->id, iter_particle->pdf);
  //    }    
 //   }
  
  }
    
    gettimeofday(&tvEnd, NULL);  
    timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
    printf("\n\n\t\tTime, in seconds, for MOVE_PARTICLE with NPV (%i) and with NX (%i) = ", NPV, NX);
    printf("%ld.%06ld\n\n", tvDiff.tv_sec, tvDiff.tv_usec);
  
    
  gets(in);
  
  printf("\n\n\t\t*** SHOW PARTICLES IN THE BOX ***\n\n");
  HASH_ITER(hh, particles, iter_head, tmp1) {
      printf("\n\t\t    Head (%i) \n", iter_head->index.ix);
    HASH_ITER(hh, iter_head->sub, iter_particle, tmp2) {
      printf("\n\t\t\t    Particle (%i) - PDF(%f) \n\n", iter_particle->id, iter_particle->pdf);
    }
  }

  HASH_ITER(hh, particles, iter_head, tmp1) {
    HASH_ITER(hh, iter_head->sub, iter_particle, tmp2) {
      HASH_DEL(iter_head->sub, iter_particle);
      free(iter_particle);
    }
    HASH_DEL(particles, iter_head);
    free(iter_head);
  }  
 
  return 0;
  
}


int id_sort(notional_particle_t *a, notional_particle_t *b) {
  return (a->id - b->id);
}


int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}


void INITIALIZE_HEADS() {  
  
  //printf("\n\n\t\t*** INITIALIZE_HEADS() ***\n\n");
  
  int i, j;
  
  notional_particle_t *ptr_head; 
  
  for(i=0; i<NX; i++) {   
    
    ptr_head = malloc(sizeof(notional_particle_t));
    ptr_head->sub = NULL;
    ptr_head->index.ix = i;  
    ptr_head->id = 0;
    ptr_head->counter_id = 0;
      
    HASH_ADD(hh, particles, index, sizeof(index_t), ptr_head);   
          
    //printf("\n\t\t    Head Particle (%i)\n", ptr_head->index.ix);
  
    ptr_head = NULL;
    
    free(ptr_head);     
    
  }  
  
}


void ADD_PARTICLE(notional_particle_t *ptr_head, int identifier, double pdf) {
       
  notional_particle_t *particle = malloc(sizeof(notional_particle_t));
  particle->sub = NULL;  
  particle->id = identifier; 
  particle->counter_id = 0; 
  particle->index.ix = ptr_head->index.ix;
  particle->pdf = pdf;
  
  HASH_ADD_INT(ptr_head->sub, id, particle); 
  
  ptr_head->counter_id = ptr_head->counter_id + 1;
  
  //printf("\n\t\t*** HEAD(%i) -> ADD_PARTICLE(%i, %i) ***\n\n", ptr_head->index.ix, particle->index.ix, particle->id);
  
  particle = NULL;
  free(particle);
   
}

void COPY_PARTICLE(notional_particle_t *ptr_head, notional_particle_t *ptr_particle, int identifier) {
  
  notional_particle_t *particle = malloc(sizeof(notional_particle_t));
  particle->sub = NULL;  
  particle->id = identifier;   
  particle->counter_id = 0;
  particle->index.ix = ptr_head->index.ix;
  particle->pdf = ptr_particle->pdf;
    
  HASH_ADD_INT(ptr_head->sub, id, particle); 
  
  ptr_head->counter_id = ptr_head->counter_id + 1;
  
  //printf("\n\n\t\t*** HEAD(%i) -> COPY_PARTICLE(%i, %i) - PDF(%f) ***\n\n", ptr_head->index.ix, particle->index.ix, particle->id, particle->pdf);
  
  particle = NULL;
  free(particle);  
  
}



/*
void print_vtk(int n_particles, int k) {
      struct lag_particles *lag_ptr;

     char string[500];
     int file_number = 0;
    
     sprintf( string, "lag_DPM_paraview_%i.vtk", k );
    

     FILE *pvtk;
     pvtk=fopen(string, "w");
     fprintf(pvtk,"# vtk DataFile Version 2.0\nPoints\nASCII\nDATASET POLYDATA\n");
     fprintf(pvtk,"Points	%i float\n",n_particles);
                   
     
      for(lag_ptr=id_table; lag_ptr != NULL; lag_ptr=(struct lag_particles*)(lag_ptr->hh.next)) {
        fprintf(pvtk,"%6.16lf	%6.16lf	%6.16lf\n",lag_ptr->x1,  lag_ptr->y1, lag_ptr->z1);
      }

     fclose(pvtk);
     
}
*/



