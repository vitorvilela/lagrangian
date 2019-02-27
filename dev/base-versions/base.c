#include <stdlib.h>	/* malloc */
#include <stddef.h>	/* offsetof */
#include <string.h>	/* memset */
#include <math.h>	/* cos, sin */
#include <stdio.h>	/* printf */
#include "uthash.h"



#define TOLERANCE 0.5
#define NPV 2					/* Number of Particle per Volume */
#define MAX_PER_VOLUME NPV*(1+TOLERANCE)
#define MIN_PER_VOLUME NPV*(1-TOLERANCE)
#define DELTA_PARTICLE 0.1

#define XI -1.0
#define XF 1.0
#define YI -1.0
#define YF 1.0
#define NX 2
#define NY 2

#define DT 0.1
#define TOTAL_TIME 1.0



typedef struct {
  int ix, iy;
} index_key_t;


typedef struct notional_particle {
  UT_hash_handle hh;  
  index_key_t index_key;
  struct notional_particle *in_volume;
  unsigned int particle_ID;  
  double value;
  double x, y;  
} notional_particle_t;



void set_index_key(notional_particle_t *, double, double); 

void set_particle(notional_particle_t *, notional_particle_t *, unsigned int, double, double); 



notional_particle_t *particles = NULL;



int main (int argc, char *argv[]) {
 
  notional_particle_t *ptr_hp, *ptr_sp; //hp: head particle, sp: sub particle
        
  notional_particle_t *particle_head, *particle_ht, *tmp1, *tmp2;
      




  /* Setting 2 Fixed Particles per Volume Inside the Eulerian Domain */
 
  // (2,2)
  // #1 	
  set_particle(ptr_hp, ptr_sp, 1, 0.5, 0.25);
 
  // #2 	
  set_particle(ptr_hp, ptr_sp, 2, 0.25, 0.5);/*
  // (1,2)
  // #3 	
  set_particle(ptr_particle, 3, -0.5, 0.25);
  // #4 	
  set_particle(ptr_particle, 4, -0.25, 0.5);
  // (1,1)
  // #5 	
  set_particle(ptr_particle, 5, -0.5, -0.25);
  // #6 	
  set_particle(ptr_particle, 6, -0.25, -0.5);
  // (2,1)
  // #7 	
  set_particle(ptr_particle, 7, 0.5, -0.25);
  // #8 	
  set_particle(ptr_particle, 8, 0.25, -0.5);

*/

  /* Iterate Over Hash Elements  */
  
  
  HASH_ITER(hh, particles, particle_head, tmp1) {
      printf("Particle: %u \n\nix = %i, iy = %i \n", particle_head->particle_ID, particle_head->index_key.ix, particle_head->index_key.iy);
    HASH_ITER(hh, particle_head->in_volume, particle_ht, tmp2) {
      printf("\t\tParticle: %u \n\nix = %i, iy = %i \n", particle_ht->particle_ID, particle_ht->index_key.ix, particle_ht->index_key.iy);
    }
  }


  /* Clean Up Both Hash Tables */
/*  
  
  HASH_ITER(hh, particles, particle_head, tmp1) {
    HASH_ITER(hh, particle_head->in_volume, particle_ht, tmp2) {
      HASH_DEL(particle_head->in_volume, particle_ht);
      free(particle_ht);
    }
    HASH_DEL(particles, particle_head);
    free(particle_head);
  }
*/


  return 0;
  
}




void set_index_key(notional_particle_t *ptr_particle, double theta, double r) {
  	
  double x, y;
    
  //x = r*cos(theta);
  //y = r*sin(theta);
   
  x = theta;
  y = r;
  
  ptr_particle->x = x;
  ptr_particle->y = y;
  
  if(x <= 0.0) {
    ptr_particle->index_key.ix = 1;
    if(y <= 0.0)
      ptr_particle->index_key.iy = 1;
    else
      ptr_particle->index_key.iy = 2;
  }
  else {
    ptr_particle->index_key.ix = 2;
    if(y <= 0.0)
      ptr_particle->index_key.iy = 1;
    else
      ptr_particle->index_key.iy = 2; 
  }
  
}

void set_particle(notional_particle_t *ptr_particle, notional_particle_t *ptr_sp, unsigned int particle_ID, double theta, double r) {
     
  
  HASH_FIND(hh, particles, &ptr_particle->index_key, sizeof(index_key_t), ptr_particle);

  if(ptr_particle == NULL) {
	ptr_particle = (notional_particle_t *)malloc( sizeof(notional_particle_t) );
  	ptr_particle->particle_ID = particle_ID;
  	ptr_particle->in_volume = NULL;
	set_index_key(ptr_particle, theta, r-DELTA_PARTICLE);
  	ptr_particle->value = 0.0;
  	HASH_ADD(hh, particles, index_key, sizeof(index_key_t), ptr_particle);
  }
  else {
  	ptr_sp = NULL;
        ptr_sp = (notional_particle_t *)malloc( sizeof(notional_particle_t) );
	ptr_sp->particle_ID = particle_ID;
	ptr_sp->in_volume = NULL;
	set_index_key(ptr_sp, theta, r-DELTA_PARTICLE);
	ptr_sp->value = 0.0;
	HASH_ADD_INT(particles, particle_ID, ptr_particle->in_volume);
  }
  
}  

