#include <stdio.h>	/* printf */
#include <stdlib.h>	/* malloc */
#include "uthash.h"


#define NX 2
#define NY 2

#define NPV 2



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



/* Hash Table Head */

notional_particle_t *particles = NULL;



int main(int argc, char *argvp[]) {
  
  int i, j, p;
  
  char in[1];
  
  int volume_counter = 0;
  
  notional_particle_t *iter_head, *iter_sub, *tmp1, *tmp2;
   
  initialize_heads();
  
  gets(in);
    
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
      
  gets(in);
  
  
  /* Start - Print Function */  
  HASH_ITER(hh, particles, iter_head, tmp1) {    
    printf("\nParticle: %i \n\nix = %i, iy = %i \n", iter_head->ID_key, iter_head->index_key.ix, iter_head->index_key.iy);
    HASH_ITER(hh, iter_head->sub, iter_sub, tmp2) {
      printf("\n\t\tParticle: %i \n\t\tix = %i, iy = %i \n\n", iter_sub->ID_key, iter_sub->index_key.ix, iter_sub->index_key.iy);
    }
  }  
  /* End - Print Function */
  
  
  
  /* Start - Find a Specific Volume to Make a Loop Over the Particles Inside It */ 
  printf("\n\n*** Find a Specific Volume to Make a Loop Over the Particles Inside It ***\n\n");  
  iter_head = (notional_particle_t *)malloc(sizeof(notional_particle_t));
  iter_sub = (notional_particle_t *)malloc(sizeof(notional_particle_t));  
  tmp1 = (notional_particle_t *)malloc(sizeof(notional_particle_t));
  
  notional_particle_t *ptr_particle = malloc(sizeof(notional_particle_t));
  ptr_particle->index_key.ix = 1;
  ptr_particle->index_key.iy = 2;
  
  HASH_FIND(hh, particles, &ptr_particle->index_key, sizeof(index_key_t), iter_head); 
  
  if(iter_head == NULL)
    printf("\n\nThere is no particle with this index_key\n\n");
  else { 
    printf("\nParticle: %i \n\nix = %i, iy = %i \n", iter_head->ID_key, iter_head->index_key.ix, iter_head->index_key.iy);
    HASH_ITER(hh, iter_head->sub, iter_sub, tmp1) {
      printf("\n\t\tParticle: %i \n\t\tix = %i, iy = %i \n\n", iter_sub->ID_key, iter_sub->index_key.ix, iter_sub->index_key.iy);
    }       
  }
  /* End - Find a Specific Volume to Make a Loop Over the Particles Inside It */
  
  
  
  /* Start - Find a Specific Particle Inside a Given Volume */ 
  printf("\n\n*** Find a Specific Particle Inside a Given Volume ***\n\n"); 
  ptr_particle = (notional_particle_t *)malloc(sizeof(notional_particle_t));
  ptr_particle->index_key.ix = 2;
  ptr_particle->index_key.iy = 2;
  ptr_particle->ID_key = 1;    
  
  iter_head = (notional_particle_t *)malloc(sizeof(notional_particle_t));
  iter_sub = (notional_particle_t *)malloc(sizeof(notional_particle_t));  
  HASH_FIND(hh, particles, &ptr_particle->index_key, sizeof(index_key_t), iter_head);
  if(iter_head == NULL)
    printf("\n\nThere is no particle with this index_key\n\n");
  else {    
    HASH_FIND_INT(iter_head->sub, &ptr_particle->ID_key, iter_sub);
    printf("\nParticle: %i \n\nix = %i, iy = %i \n", iter_head->ID_key, iter_head->index_key.ix, iter_head->index_key.iy);
    if(iter_sub == NULL)
      printf("\n\nThere is no particle with this ID_key\n\n");
    else 
      printf("\n\t\tParticle: %i \n\t\tix = %i, iy = %i \n\n", iter_sub->ID_key, iter_sub->index_key.ix, iter_sub->index_key.iy);       
  }  
  /* End - Find a Specific Particle Inside a Given Volume */ 
  
  
  
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
  
  printf("\n\n*** Get Head ***\n\n");
  
  index_key_t index_key;
  index_key.ix = ix;
  index_key.iy = iy;
  
  notional_particle_t *ptr_head = malloc(sizeof(notional_particle_t));
  HASH_FIND(hh, particles, &index_key, sizeof(index_key_t), ptr_head);
  
  return ptr_head;
}

/* Add Particle */

void add_particle(notional_particle_t *ptr_head, int ID_key) {
  
  printf("\n\n*** Add Particle ***\n\n");
  
  notional_particle_t *particle = malloc(sizeof(notional_particle_t));
  particle->sub = NULL;  
  particle->ID_key = ID_key;
  HASH_ADD_INT(ptr_head->sub, ID_key, particle);  
  particle->index_key.ix = ptr_head->index_key.ix;
  particle->index_key.iy = ptr_head->index_key.iy;
  
  //printf("\n\t\tParticle: %i \n\t\tix = %i, iy = %i \n\n", particle->ID_key, particle->index_key.ix, particle->index_key.iy);
  
}
    
    






