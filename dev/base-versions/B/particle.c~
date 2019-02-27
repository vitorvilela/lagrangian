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

notional_particle_t *particles = NULL;

int main(int argc, char *argvp[]) {
  
  int i, j, p;
  
  int volume_counter = 0;
  
  notional_particle_t *iter_head, *iter_sub, *tmp1, *tmp2;
  
  
  
  /* Start - Initialize Particles */    
  /* Initialize All Head Particles as a Function of NX and NY */  
  printf("\n\n*** Initialize Particles ***\n\n");
  for(i=1; i<=NX; i++)
    for(j=1; j<=NY; j++) {      
      notional_particle_t *ptr_hp = malloc(sizeof(notional_particle_t));
      ptr_hp->sub = NULL;
      ptr_hp->index_key.ix = i;
      ptr_hp->index_key.iy = j;
      ptr_hp->ID_key = 0; 
      HASH_ADD(hh, particles, index_key, sizeof(index_key_t), ptr_hp);
      
      /* Initialize All Particles Inside the Control Volumes as a Function of NPV */
      
      for(p=1; p<=NPV; p++) {
	notional_particle_t *ptr_sp = malloc(sizeof(notional_particle_t));
	ptr_sp->sub = NULL;
	ptr_sp->index_key.ix = i;
	ptr_sp->index_key.iy = j;
	ptr_sp->ID_key = volume_counter+1;
	HASH_ADD_INT(ptr_hp->sub, ID_key, ptr_sp);
	volume_counter += 1;
      }
      volume_counter = 0;    
    } 
  /* End - Initialize Particles */ 

  
  
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
