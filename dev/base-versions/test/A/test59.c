#include <stdio.h>
#include <stdlib.h>
#include "uthash.h"

typedef struct {
  int ix, iy;
} index_key_t;

typedef struct notional_particle {
  UT_hash_handle hh; 
  index_key_t index_key;
  struct notional_particle *sub;
  int particle_ID;  
} notional_particle_t;

notional_particle_t *particles = NULL;

int main(int argc, char *argvp[]) {
  
  notional_particle_t *particle_head, *particle_ht, *tmp1, *tmp2;
  
  notional_particle_t *ptr_hp = malloc(sizeof(*ptr_hp));
  ptr_hp->sub = NULL;
  ptr_hp->index_key.ix = 1;
  ptr_hp->index_key.iy = 1;
  ptr_hp->particle_ID = 0;
  
  HASH_ADD(hh, particles, index_key, sizeof(index_key_t), ptr_hp);
  
  ptr_hp = NULL;
  free(ptr_hp);
  
  ptr_hp = malloc(sizeof(*ptr_hp));
  index_key_t index_key;
  index_key.ix = 1;
  index_key.iy = 1;
  HASH_FIND(hh, particles, &index_key, sizeof(index_key_t), ptr_hp);
  
  notional_particle_t *ptr_sp = malloc(sizeof(*ptr_sp));
  ptr_sp->sub = NULL;
  ptr_sp->index_key.ix = 1;
  ptr_sp->index_key.iy = 1;
  ptr_sp->particle_ID = 1;
  HASH_ADD_INT(ptr_hp->sub, particle_ID, ptr_sp);

  ptr_sp = NULL;
  free(ptr_sp);
    
  ptr_sp = malloc(sizeof(*ptr_sp));
  ptr_sp->sub = NULL;
  ptr_sp->index_key.ix = 1;
  ptr_sp->index_key.iy = 1;
  ptr_sp->particle_ID = 2;
  HASH_ADD_INT(ptr_hp->sub, particle_ID, ptr_sp);
  
  ptr_sp = NULL;
  free(ptr_sp);
  
  ptr_hp = NULL;
  free(ptr_hp);
  
  HASH_ITER(hh, particles, particle_head, tmp1) {
      printf("Particle: %i \n\nix = %i, iy = %i \n", particle_head->particle_ID, particle_head->index_key.ix, particle_head->index_key.iy);
    HASH_ITER(hh, particle_head->sub, particle_ht, tmp2) {
      printf("\t\tParticle: %i \n\n\t\tix = %i, iy = %i \n", particle_ht->particle_ID, particle_ht->index_key.ix, particle_ht->index_key.iy);
    }
  }
  
  
  HASH_ITER(hh, particles, particle_head, tmp1) {
    HASH_ITER(hh, particle_head->sub, particle_ht, tmp2) {
      HASH_DEL(particle_head->sub, particle_ht);
      free(particle_ht);
    }
    HASH_DEL(particles, particle_head);
    free(particle_head);
  }

  return 0;
}
