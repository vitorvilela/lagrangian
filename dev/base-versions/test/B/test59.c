#include <stdio.h>
#include <stdlib.h>
#include "uthash.h"

#define NX 1
#define NPV 1

typedef struct {
  int ix;
} index_t;

typedef struct notional_particle {
  UT_hash_handle hh; 
  index_t index;
  struct notional_particle *sub;
  int id;  
} notional_particle_t;


void INITIALIZE_HEADS();

void GET_HEAD(notional_particle_t *, int);



notional_particle_t *particles = NULL;

int main(int argc, char *argvp[]) {
  
  int i, p;
  
  notional_particle_t *ptr_head, *ptr_particle, *tmp1, *tmp2;
  
  
  INITIALIZE_HEADS();
  
//  ptr_head = malloc(sizeof(*ptr_head));
//  ptr_head->sub = NULL;
//  ptr_head->index.ix = 0;  
//  ptr_head->id = 0;
//  HASH_ADD(hh, particles, index, sizeof(index_t), ptr_head);
  
//  ptr_head = NULL;
//  free(ptr_head);
  
  i = 0;
  ptr_head = malloc(sizeof(notional_particle_t));
  GET_HEAD(ptr_head, i);
  
//  ptr_head = malloc(sizeof(*ptr_head));
//  index_t index;
//  index.ix = 0; 
//  HASH_FIND(hh, particles, &index, sizeof(index_t), ptr_head);
  
  ptr_particle = malloc(sizeof(*ptr_particle));
  ptr_particle->sub = NULL;
  ptr_particle->index.ix = 0;  
  ptr_particle->id = 1;
  HASH_ADD_INT(ptr_head->sub, id, ptr_particle);

  ptr_particle = NULL;
  free(ptr_particle);
    
  ptr_particle = malloc(sizeof(*ptr_particle));
  ptr_particle->sub = NULL;
  ptr_particle->index.ix = 0;  
  ptr_particle->id = 2;
  HASH_ADD_INT(ptr_head->sub, id, ptr_particle);
  
  ptr_particle = NULL;
  free(ptr_particle);
  
  ptr_head = NULL;
  free(ptr_head);
  
  
  printf("\n\n\t\t*** SHOW PARTICLES IN THE BOX ***\n\n");
  HASH_ITER(hh, particles, ptr_head, tmp1) {
    printf("\t\t    Head (%i) \n", ptr_head->index.ix);
    HASH_ITER(hh, ptr_head->sub, ptr_particle, tmp2) {
      printf("\n\t\t\t    Particle (%i) \n", ptr_particle->id);
    }
  }
  
  
  HASH_ITER(hh, particles, ptr_head, tmp1) {
    HASH_ITER(hh, ptr_head->sub, ptr_particle, tmp2) {
      HASH_DEL(ptr_head->sub, ptr_particle);
      free(ptr_particle);
    }
    HASH_DEL(particles, ptr_head);
    free(ptr_head);
  }

  return 0;
}



void INITIALIZE_HEADS() {  
  
  printf("\n\n\t\t*** INITIALIZE_HEADS() ***\n");
  
  int i;
  
  notional_particle_t *ptr_head; 
  
  for(i=0; i<NX; i++) {   
    
    ptr_head = malloc(sizeof(notional_particle_t));
    ptr_head->sub = NULL;
    ptr_head->index.ix = i;  
    ptr_head->id = 0;
      
    HASH_ADD(hh, particles, index, sizeof(index_t), ptr_head);   
          
    printf("\n\t\t    Head Particle (%i)\n", ptr_head->index.ix);
  
    ptr_head = NULL;
    
    free(ptr_head);     
    
  }  
  
}


void GET_HEAD(notional_particle_t *ptr_head, int i) {
    
  index_t index;
  index.ix = i;
  
  HASH_FIND(hh, particles, &index, sizeof(index_t), ptr_head);
  
  if(ptr_head == NULL) 
    printf("\n\n\t\t*** There is no Head with This Index ***\n\n"); 
  else
    printf("\n\n\t\t*** GET_HEAD(%i) ***\n", ptr_head->index.ix);
     
}
