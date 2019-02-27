#include <stdio.h>	/* printf */
#include <stdlib.h>	/* malloc */
#include <time.h>	/* clock */
#include "uthash.h"

#define X_LENGTH 1.0
#define Y_LENGTH 1.0
#define NX 2
#define NY 1
#define DX X_LENGTH/NX
#define DY Y_LENGTH/NY

#define max(a, b) ((a)>(b)?(a):(b))
#define min(a, b) ((a)<(b)?(a):(b))

#define bool int
#define true 1
#define false 0

#define NPV 1



typedef struct {
  int ix, iy;
} index_t;

typedef struct {
  double xb_i, yb_i, xb_f, yb_f;
} bound_t;

typedef struct {
  double xp, yp;
} position_t;

typedef struct notional_particle {
  UT_hash_handle hh; 
  index_t index;
  bound_t bound;
  position_t position;
  struct notional_particle *sub;
  int last_id;				/* Just for head particle */
  int id;    
  
} notional_particle_t;



/* Functions Heads */

void initialize_heads(double *, double *);

notional_particle_t *get_head(int, int);

void add_particle(notional_particle_t *, int); 

//void move_particle(notional_particle_t *, int, position_t);

notional_particle_t *get_particle(notional_particle_t *, int);

void move_particle(notional_particle_t *, notional_particle_t *);

int id_sort(notional_particle_t *, notional_particle_t *);


/* Hash Table Head */

notional_particle_t *particles = NULL;



int main(int argc, char *argvp[]) {
  
  
  time_t time1, time2;
  
  int test;
  
  double diff;
    
  int i, j, p, t;
   
  char in[1];
    
  int volume_counter = 0;
  
  double X[NX+1], Y[NY+1];
  
  for(i=0; i<=NX; i++)
    X[i] = i*DX;
    
  for(j=0; j<=NY; j++)     
    Y[j] = j*DY; 
    
  notional_particle_t *particle, *ptr_particle, *iter_head, *iter_sub, *tmp1, *tmp2;
  
  index_t index;
   
  initialize_heads(X, Y); 
  gets(in);
    
  for(i=1; i<NX; i++) {  //alterado i<=NX
    for(j=1; j<=NY; j++) { 
      iter_head = malloc(sizeof(notional_particle_t)); 
      iter_head = get_head(i, j);
      for(p=1; p<=NPV; p++) {
	add_particle(iter_head, volume_counter+1);  
	volume_counter += 1;  
	gets(in);
      }
      iter_head->last_id = volume_counter;
      printf("\t\tHeader (%i, %i) -> last_id = %i \n\n", i, j, iter_head->last_id);
      volume_counter = 0;   
      gets(in);
    }          
  }
  
    
  time1 = time(NULL);
         
 // printf("\n\n\t\t*** Find a Specific Particle Inside a Given Volume ***\n\n");   
 // iter_head = malloc(sizeof(notional_particle_t)); 
 // iter_head = get_head(1, 1);
 // particle = malloc(sizeof(notional_particle_t)); 
 // particle = get_particle(iter_head, NPV);
  
  
//  time2 = time(NULL);
    
//  diff = difftime(time2, time1);
  
//  printf("\n\n\t\tTime for NPV %i = %f \n\n", NPV, diff);
  
  gets(in);
  
  printf("\n\n\t\t*** Linear Transport of the Particles in X-Direction ***\n\n");   
  int particle_step = 10;
  double particle_velocity = DX/particle_step;
  printf("\t\tParticle Velocity = %f\n\n", particle_velocity);
  int iterations = particle_step*NX;
  
  notional_particle_t *iterator = malloc(sizeof(notional_particle_t));     
  notional_particle_t *iterator1 = malloc(sizeof(notional_particle_t));     
  
  for(t=1; t<=iterations; t++) {
    printf("\n\n\t\t*** Step: %i ***\n\n", t);
    for(i=1; i<=NX; i++) {
      for(j=1; j<=NY; j++) { 
	iter_head = malloc(sizeof(notional_particle_t)); 
	iter_head = get_head(i, j);
	//for(p=1; p<=iter_head->last_id; p++) {             /* problema com este loop - verificar iteração dentro da ht */
	for(iterator=iter_head->sub; iterator!=NULL; iterator=(struct notional_particle*)(iterator->hh.next)) {
	  
	  ptr_particle = malloc(sizeof(notional_particle_t)); 
	  
	  //int id = iterator->id;
	  
	  //ptr_particle = get_particle(iter_head, id);
	  ptr_particle = iterator;
	  	  
	  /* Move the Particle */
	  printf("\t\tParticle changes its position\n");
	  ptr_particle->position.xp += particle_velocity;
	  printf("\n\t\tParticle: id(%i) - last_id(%i) \n\t\tix = %i, iy = %i \n\t\txp = %f, yp = %f \n\t\txb_i = %f, xb_f = %f\n\t\tyb_i = %f, yb_f = %f\n\n", ptr_particle->id, iter_head->last_id, ptr_particle->index.ix, ptr_particle->index.iy, ptr_particle->position.xp, ptr_particle->position.yp, ptr_particle->bound.xb_i, ptr_particle->bound.xb_f, ptr_particle->bound.yb_i, ptr_particle->bound.yb_f); 
	  gets(in); 
	  
	  /* Verify Index and Update it if Necessary */
	  if(ptr_particle->position.xp > iter_head->bound.xb_f) {
	  printf("\n\n\t\t*** Particle goes out of the box *** \n\n");
	  gets(in);  
	    if(ptr_particle->position.xp > X_LENGTH) {
	      printf("\n\n\t\t*** Particle goes out of the domain *** \n\n");	 
	      gets(in);
	      printf("WAIT");
	      gets(in);
	      printf("\n\n\t\t*** HASH DEL and free *** \n\n");	
	      HASH_DEL(iter_head->sub, ptr_particle);        /* Problema com esta definição de DEL - ptr_particle é iter_head->sub */
	      free(ptr_particle);	      
	      gets(in);
	      printf("\n\n\t\t*** HASH SORT and Count *** \n\n");
	      HASH_SORT(iter_head->sub, id_sort);
	      iter_head->last_id = HASH_COUNT(iter_head->sub);
	      gets(in);
	      printf("\n\n\t\t*** To show the remainder particles in the box *** \n\n");	  
	      for(iterator1=iter_head->sub; iterator1!=NULL; iterator1=(struct notional_particle*)(iterator1->hh.next)) {				
		printf("\n\t\tParticle: id(%i) - last_id(%i) \n\t\tix = %i, iy = %i \n\t\txp = %f, yp = %f \n\t\txb_i = %f, xb_f = %f\n\t\tyb_i = %f, yb_f = %f\n\n", iterator1->id, iter_head->last_id, iterator1->index.ix, iterator1->index.iy, iterator1->position.xp, iterator1->position.yp, iterator1->bound.xb_i, iterator1->bound.xb_f, iterator1->bound.yb_i, iterator1->bound.yb_f); 
	        gets(in);		
	      }
	      gets(in);
	    }	    
	    
	    else {    
	      printf("\n\n\t\t*** Index Update *** \n\n");	
	      gets(in);
	      ptr_particle->index.ix += 1;    
	      printf("\n\t\tParticle: id(%i)\n\t\t(new) ix = %i, iy = %i \n\t\txp = %f, yp = %f \n\t\txb_i = %f, xb_f = %f\n\t\tyb_i = %f, yb_f = %f\n\n", ptr_particle->id, ptr_particle->index.ix, ptr_particle->index.iy, ptr_particle->position.xp, ptr_particle->position.yp, ptr_particle->bound.xb_i, ptr_particle->bound.xb_f, ptr_particle->bound.yb_i, ptr_particle->bound.yb_f); 
	      gets(in);
	      tmp1 = malloc(sizeof(notional_particle_t));		/* Point out to the new address (box) */
	      tmp1 = get_head(ptr_particle->index.ix, ptr_particle->index.iy);
	      printf("\n\t\tNew header last_id(%i) - before particle adition", tmp1->last_id);	      
	      move_particle(tmp1, ptr_particle);
	      printf("\n\t\tNew header last_id(%i) - after particle adition", tmp1->last_id);
	      gets(in);
	      /* Delete the particle (ptr_particle) from old box (iter_head) */ 
	      printf("\n\n\t\t*** HASH DEL and free *** \n\n");	
	      HASH_DEL(iter_head->sub, ptr_particle);
	      free(ptr_particle); 	      
	      printf("\n\n\t\t*** HASH SORT and Count *** \n\n");
	      HASH_SORT(iter_head->sub, id_sort);
	      iter_head->last_id = HASH_COUNT(iter_head->sub);
	      gets(in);
	      
	    }	    
	  
	  }
	  free(ptr_particle);
	}
	
	//free(iter_head);
      }
    }   
  }
  
  
  /* Clean up the Hash Table */  
  HASH_ITER(hh, particles, iter_head, tmp1) {
    HASH_ITER(hh, iter_head->sub, iter_sub, tmp2) {
      HASH_DEL(iter_head->sub, iter_sub);
      free(iter_sub);
    }
    HASH_DEL(particles, iter_head);
    free(iter_head);
  }  
 

  
  return 0;
}



int id_sort(notional_particle_t *a, notional_particle_t *b) {
    return (a->id - b->id);
}


/* Initialize All Head Particles as a Function of NX and NY */  

void initialize_heads(double *ptr_x, double *ptr_y) {  
  
  printf("\n\n*** Initialize Heads ***\n\n");
  
  int i, j;
  
  notional_particle_t *ptr_head; 
  
  for(i=1; i<=NX; i++)
    for(j=1; j<=NY; j++) {      
      ptr_head = malloc(sizeof(notional_particle_t));
      ptr_head->sub = NULL;
      ptr_head->index.ix = i;
      ptr_head->index.iy = j;
      ptr_head->id = 0;
      ptr_head->bound.xb_i = ptr_x[i-1];
      ptr_head->bound.xb_f = ptr_x[i];
      ptr_head->bound.yb_i = ptr_y[j-1];
      ptr_head->bound.yb_f = ptr_y[j];
      ptr_head->position.xp = ptr_x[i-1] + 0.5*DX;	/* The head particle is positioned */
      ptr_head->position.yp = ptr_y[j-1] + 0.5*DY;	/* at the volume control centre */
      
      HASH_ADD(hh, particles, index, sizeof(index_t), ptr_head);   
      
      printf("\nParticle: id(%i) - last_id(%i) \nix = %i, iy = %i \nxp = %f, yp = %f\nxb_i = %f, xb_f = %f\nyb_i = %f, yb_f = %f\n\n", ptr_head->id, ptr_head->last_id, ptr_head->index.ix, ptr_head->index.iy, ptr_head->position.xp, ptr_head->position.yp, ptr_head->bound.xb_i,  ptr_head->bound.xb_f,  ptr_head->bound.yb_i,  ptr_head->bound.yb_f);
  }
  
}



/* Get ptr to Head  */  

notional_particle_t *get_head(int ix, int iy) {
  
  printf("\n\n\t\t*** Get Head (%i, %i) ***\n\n", ix, iy);
  
  index_t index;
  index.ix = ix;
  index.iy = iy;
  
  notional_particle_t *ptr_head = malloc(sizeof(notional_particle_t));
  HASH_FIND(hh, particles, &index, sizeof(index_t), ptr_head);
  
  if(ptr_head == NULL) 
    printf("\n\n*** There is no Head with This Index ***\n\n");  
    
  return ptr_head;
  
}

/* Add Particle */

void add_particle(notional_particle_t *ptr_head, int id) {
  
  printf("\n\n\t\t*** Add Particle ***\n\n");
  
  notional_particle_t *particle = malloc(sizeof(notional_particle_t));
  particle->sub = NULL;  
  particle->id = id;
  particle->last_id = 0;
  HASH_ADD_INT(ptr_head->sub, id, particle);  
  particle->index.ix = ptr_head->index.ix;
  particle->index.iy = ptr_head->index.iy;
  particle->bound.xb_i = ptr_head->bound.xb_i;
  particle->bound.xb_f = ptr_head->bound.xb_f;
  particle->bound.yb_i = ptr_head->bound.yb_i;
  particle->bound.yb_f = ptr_head->bound.yb_f;
  particle->position.xp = ptr_head->position.xp;
  particle->position.yp = ptr_head->position.yp;  
  
  
    
  printf("\n\t\tParticle: id(%i) - last_id of header(%i)\n\t\tix = %i, iy = %i \n\t\txp = %f, yp = %f\n\t\txb_i = %f, xb_f = %f\n\t\tyb_i = %f, yb_f = %f\n\n", particle->id, ptr_head->last_id, particle->index.ix, particle->index.iy, particle->position.xp, particle->position.yp, particle->bound.xb_i, particle->bound.xb_f, particle->bound.yb_i, particle->bound.yb_f);
  free(particle);
}
       
    
void move_particle(notional_particle_t *ptr_head, notional_particle_t *ptr_particle) {
  
  printf("\n\n\t\t*** Move Particle ***\n\n");
  
  ptr_head->last_id = ptr_head->last_id+1;  
  int id = ptr_head->last_id;
  ptr_particle->id = id;
  
  HASH_ADD_INT(ptr_head->sub, id, ptr_particle); 
    
  
  ptr_particle->bound.xb_i = ptr_head->bound.xb_i;
  ptr_particle->bound.xb_f = ptr_head->bound.xb_f;
  ptr_particle->bound.yb_i = ptr_head->bound.yb_i;
  ptr_particle->bound.yb_f = ptr_head->bound.yb_f;
 
    
  printf("\n\t\tParticle Moves to: id(%i) \n\t\tix = %i, iy = %i \n\t\txp = %f, yp = %f\n\t\txb_i = %f, xb_f = %f\n\t\tyb_i = %f, yb_f = %f\n\n", id, ptr_particle->index.ix, ptr_particle->index.iy, ptr_particle->position.xp, ptr_particle->position.yp, ptr_particle->bound.xb_i, ptr_particle->bound.xb_f, ptr_particle->bound.yb_i, ptr_particle->bound.yb_f);
  
}


/* Get Particle */  

notional_particle_t *get_particle(notional_particle_t *ptr_head, int id) {
  
  printf("\n\n\t\t*** Get Particle id(%i) ***\n\n", id);
  
  notional_particle_t *ptr_particle = malloc(sizeof(notional_particle_t));  
    
  HASH_FIND_INT(ptr_head->sub, &id, ptr_particle);
  
  if(ptr_particle == NULL) 
    printf("\n\n\t\t*** There is no Particle with This ID ***\n\n");
  else
    printf("\n\t\tParticle: id(%i) - last_id(%i) \n\t\tix = %i, iy = %i \n\t\txp = %f, yp = %f \n\t\txb_i = %f, xb_f = %f\n\t\tyb_i = %f, yb_f = %f\n\n", ptr_particle->id, ptr_head->last_id, ptr_particle->index.ix, ptr_particle->index.iy, ptr_particle->position.xp, ptr_particle->position.yp, ptr_particle->bound.xb_i, ptr_particle->bound.xb_f, ptr_particle->bound.yb_i, ptr_particle->bound.yb_f);   
  
  return ptr_particle;
  
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



