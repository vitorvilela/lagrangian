#include <cstdio>	/* printf */
#include <cstdlib>	/* malloc */
#include <ctime>	/* clock */
#include "uthash.h"



#define NX 1
#define NY 1

#define NPV 100000


using namespace std;



typedef struct {
  int ix, iy;
} index_key_t;

typedef struct notional_particle {
  UT_hash_handle hh; 
  index_key_t index_key;
  struct notional_particle *sub;
  int ID_key;    
} notional_particle_t;



/* Functions Heads */

void initialize_heads();

notional_particle_t *get_head(int, int);

void add_particle(notional_particle_t *, int); 

notional_particle_t *get_particle(notional_particle_t *, int);



/* Hash Table Head */

notional_particle_t *particles = NULL;



int main(int argc, char *argvp[]) {
  
  int time1, time2;
  //float diff = 1.0;
  
  clock_t start;
  double diff;
  
  int i, j, p;
   
  char in[1];
    
  int volume_counter = 0;
  
  notional_particle_t *particle, *iter_head, *iter_sub, *tmp1, *tmp2;
   
  initialize_heads();
      
  for(i=1; i<=NX; i++) {
    for(j=1; j<=NY; j++) { 
      iter_head = malloc(sizeof(notional_particle_t)); 
      iter_head = get_head(i, j);
      for(p=1; p<=NPV; p++) {
	add_particle(iter_head, volume_counter+1);  
	volume_counter += 1;  
      }
      volume_counter = 0;  
    }    
      
  }
      
  //time1 = clock();
  start = clock();
  
  
  int test;  
  for(i=1; i<=NPV; i++) {
    test = i;
  }
  
  
 
  
  
   
    
   
  printf("\n\n*** Find a Specific Particle Inside a Given Volume ***\n\n");   
  iter_head = malloc(sizeof(notional_particle_t)); 
  iter_head = get_head(1, 1);
  particle = malloc(sizeof(notional_particle_t)); 
  particle = get_particle(iter_head, NPV/2);
  
  //time2 = clock();
  //diff = (float)(time2-time1)/(float)CLOCKS_PER_SEC;
  
  diff = ( std::clock() - start ) / (double)CLOCKS_PER_SEC;
  
  
  
  //printf("\n\nTime for NPV %i = %f \n\n", NPV, diff);
  //std::cout<<"printf: "<< diff <<'\n';
  
  
  /* Start - Clean up the Hash Table */  
  HASH_ITER(hh, particles, iter_head, tmp1) {
    HASH_ITER(hh, iter_head->sub, iter_sub, tmp2) {
      HASH_DEL(iter_head->sub, iter_sub);
      free(iter_sub);
    }
    HASH_DEL(particles, iter_head);
    free(iter_head);
  }  
  /* End - Clean up the Hash Table */

  
  return 0;
}






/* Initialize All Head Particles as a Function of NX and NY */  

void initialize_heads() {  
  
  printf("\n\n*** Initialize Heads ***\n\n");
  
  int i, j;
  
  notional_particle_t *ptr_head; 
  
  for(i=1; i<=NX; i++)
    for(j=1; j<=NY; j++) {      
      ptr_head = malloc(sizeof(notional_particle_t));
      ptr_head->sub = NULL;
      ptr_head->index_key.ix = i;
      ptr_head->index_key.iy = j;
      ptr_head->ID_key = 0;       
      HASH_ADD(hh, particles, index_key, sizeof(index_key_t), ptr_head);   
      
      printf("\nParticle: %i \n\nix = %i, iy = %i \n", ptr_head->ID_key, ptr_head->index_key.ix, ptr_head->index_key.iy);
  }
  
}



/* Get ptr to Head  */  

notional_particle_t *get_head(int ix, int iy) {
  
  //printf("\n\n*** Get Head ***\n\n");
  
  index_key_t index_key;
  index_key.ix = ix;
  index_key.iy = iy;
  
  notional_particle_t *ptr_head = malloc(sizeof(notional_particle_t));
  HASH_FIND(hh, particles, &index_key, sizeof(index_key_t), ptr_head);
  
  if(ptr_head == NULL) 
    printf("\n\n*** There is no Head with This Index ***\n\n");  
    
  return ptr_head;
  
}

/* Add Particle */

void add_particle(notional_particle_t *ptr_head, int ID_key) {
  
  //printf("\n\n*** Add Particle ***\n\n");
  
  notional_particle_t *particle = malloc(sizeof(notional_particle_t));
  particle->sub = NULL;  
  particle->ID_key = ID_key;
  HASH_ADD_INT(ptr_head->sub, ID_key, particle);  
  particle->index_key.ix = ptr_head->index_key.ix;
  particle->index_key.iy = ptr_head->index_key.iy;
    
  //printf("\n\t\tParticle: %i \n\t\tix = %i, iy = %i \n\n", particle->ID_key, particle->index_key.ix, particle->index_key.iy);
  
}
       
    
/* Get Particle */  

notional_particle_t *get_particle(notional_particle_t *ptr_head, int ID_key) {
  
  printf("\n\n*** Get Particle ***\n\n");
  
  notional_particle_t *ptr_particle = malloc(sizeof(notional_particle_t));  
    
  HASH_FIND_INT(ptr_head->sub, &ID_key, ptr_particle);
  
  if(ptr_particle == NULL) 
    printf("\n\n*** There is no Particle with This ID ***\n\n");
  else
    printf("\n\t\tParticle: %i \n\t\tix = %i, iy = %i \n\n", ptr_particle->ID_key, ptr_particle->index_key.ix, ptr_particle->index_key.iy);   
  
  return ptr_particle;
  
}





