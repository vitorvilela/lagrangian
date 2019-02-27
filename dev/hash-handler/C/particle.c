#include <stdio.h>	/* printf */
#include <stdlib.h>	/* malloc */
#include <time.h>	/* time */
#include <sys/time.h>	/* time */
#include "uthash.h"

#define NX 2
#define NPV 2

#define max(a, b) ((a)>(b)?(a):(b))
#define min(a, b) ((a)<(b)?(a):(b))

#define bool int
#define true 1
#define false 0





typedef struct {
  int ix;
} index_t;

typedef struct volume {
  index_t index;
  double x0, x1, y0, y1, z0, z1;
  double xc, yc, zc;
  double dx, dy, dz; 
  int last_used_id;
  struct notional_particle *particle;
  UT_hash_handle hh;
} volume_t;

typedef struct notional_particle {
  index_t index;
  int field_id;  
  double pdf;
  UT_hash_handle hh;   
} notional_particle_t;


void INITIALIZE_HEADS();

void ADD_PARTICLE(volume_t *, int, double); 

void COPY_PARTICLE(volume_t *, notional_particle_t *, int);

int id_sort(notional_particle_t *, notional_particle_t *);

int timeval_subtract(struct timeval *, struct timeval *, struct timeval *);



volume_t *eulerian = NULL;

int main(int argc, char *argvp[]) {
    
  int i, p, v;
  
  index_t index;
  
  struct timeval tvBegin, tvEnd, tvDiff;
   
  char in[1];
    
  int particle_counter = 0;
      
  volume_t *ptr_volume, *toVolume, *iter_volume, *tmp1;
  
  notional_particle_t *ptr_particle, *iter_particle, *tmp2;
    
  
  INITIALIZE_HEADS();
  
  gets(in);
  
  for(v=0; v<NX; v++) {
    
    ptr_volume = malloc(sizeof(volume_t));   
    index.ix = v;  
    HASH_FIND(hh, eulerian, &index, sizeof(index_t), ptr_volume);
    if(ptr_volume == NULL) 
      printf("\n\n\t\t*** There is no Volume with this Index ***\n\n");
    else
      printf("\n\n\t\t*** GET VOLUME (%i) ***\n\n", ptr_volume->index.ix);
   
    particle_counter = 0;
    for(p=1; p<=NPV; p++) {
      ADD_PARTICLE(ptr_volume, particle_counter+1, 10.0);
      particle_counter += 1;  
    }
      
    ptr_volume = NULL;
    free(ptr_volume);
  }
  
  gets(in);

  
  //    *****************  MOVE PARTICLE TO OTHER VOLUME  ***********************
  
    int index_from = 0;
    int index_to = 1;
    int move_particle = 2;
    int new_particle_id;
  
    ptr_volume = malloc(sizeof(volume_t));  
    index.ix = index_from;
    HASH_FIND(hh, eulerian, &index, sizeof(index_t), ptr_volume);
    if(ptr_volume == NULL) 
      printf("\n\n\t\t*** There is no Volume with this Index ***\n\n");
    else
      printf("\n\n\t\t*** GET HEAD (%i) ***\n\n", ptr_volume->index.ix);
    
    gets(in);
    
    ptr_particle = malloc(sizeof(notional_particle_t));  
    p = move_particle;  
    HASH_FIND_INT(ptr_volume->particle, &p, ptr_particle);  
    if(ptr_particle == NULL) 
      printf("\n\n\t\t*** There is no Particle with this id ***\n\n");
    else
      printf("\n\n\t\t*** GET PARTICLE (%i, %i) ***\n\n", ptr_volume->index.ix, ptr_particle->field_id);
  
    toVolume = malloc(sizeof(volume_t));   
    index.ix = 1;
    HASH_FIND(hh, eulerian, &index, sizeof(index_t), toVolume);
    if(toVolume == NULL) 
      printf("\n\n\t\t*** There is no Volume with this Index ***\n\n");
    else
      printf("\n\n\t\t*** GET HEAD (%i) ***\n\n", toVolume->index.ix);
     
    new_particle_id = toVolume->last_used_id + 1;
    
    COPY_PARTICLE(toVolume, ptr_particle, new_particle_id);
   
    gets(in);
  
    printf("\n\n\t\t*** SHOW PARTICLES IN THE BOX ***\n\n");
    HASH_ITER(hh, eulerian, iter_volume, tmp1) {
      printf("\n\t\t    Volume (%i) \n", iter_volume->index.ix);
      HASH_ITER(hh, iter_volume->particle, iter_particle, tmp2) {
	printf("\n\t\t\t    Particle (%i) - PDF(%f) \n\n", iter_particle->field_id, iter_particle->pdf);
      }
    }
  
    gets(in);  
  
    HASH_DEL(ptr_volume->particle, ptr_particle);
    printf("\n\n\t\t*** HASH_DEL (%i, %i) ***\n\n", ptr_volume->index.ix, ptr_particle->field_id);
    ptr_particle = NULL;
    free(ptr_particle); 
      
    ptr_volume = NULL;
    free(ptr_volume);
  
    toVolume = NULL;
    free(toVolume);
    
    gets(in); 
    
    printf("\n\n\t\t*** SHOW PARTICLES IN THE BOX ***\n\n");
    HASH_ITER(hh, eulerian, iter_volume, tmp1) {
      printf("\n\t\t    Volume (%i) \n", iter_volume->index.ix);
      HASH_ITER(hh, iter_volume->particle, iter_particle, tmp2) {
	printf("\n\t\t\t    Particle (%i) - PDF(%f) \n\n", iter_particle->field_id, iter_particle->pdf);
      }
    }
  
    //    *****************  MOVE PARTICLE TO OTHER VOLUME  ***********************
  

  HASH_ITER(hh, eulerian, iter_volume, tmp1) {
    HASH_ITER(hh, iter_volume->particle, iter_particle, tmp2) {
      HASH_DEL(iter_volume->particle, iter_particle);
      free(iter_particle);
    }
    HASH_DEL(eulerian, iter_volume);
    free(iter_volume);
  }  
 
  return 0;
  
}


int id_sort(notional_particle_t *a, notional_particle_t *b) {
  return (a->field_id - b->field_id);
}


int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}


void INITIALIZE_HEADS() {  
  
  printf("\n\n\t\t*** INITIALIZE_HEADS() ***\n\n");
  
  int i, j;
  
  volume_t *ptr_volume; 
  
  for(i=0; i<NX; i++) {   
    
    ptr_volume = malloc(sizeof(volume_t));
    ptr_volume->particle = NULL;
    ptr_volume->index.ix = i;  
    ptr_volume->last_used_id = 0;
         
    HASH_ADD(hh, eulerian, index, sizeof(index_t), ptr_volume);   
          
    printf("\n\t\t    Volume (%i)\n", ptr_volume->index.ix);
  
    ptr_volume = NULL;
    
    free(ptr_volume);     
    
  }  
  
}


void ADD_PARTICLE(volume_t *ptr_volume, int id, double pdf) {
       
  notional_particle_t *particle = malloc(sizeof(notional_particle_t));
  particle->field_id = id;
  particle->index.ix = ptr_volume->index.ix;
  particle->pdf = pdf;
  
  HASH_ADD_INT(ptr_volume->particle, field_id, particle); 
  
  ptr_volume->last_used_id = ptr_volume->last_used_id + 1;
  
  printf("\n\t\t*** VOLUME (%i) -> ADD_PARTICLE(%i, %i) ***\n\n", ptr_volume->index.ix, particle->index.ix, particle->field_id);
  
  particle = NULL;
  free(particle);
   
}


void COPY_PARTICLE(volume_t *ptr_volume, notional_particle_t *ptr_particle, int id) {
  
  notional_particle_t *particle = malloc(sizeof(notional_particle_t));
  particle->field_id = id;  
  particle->index.ix = ptr_volume->index.ix;
  particle->pdf = ptr_particle->pdf;
    
  HASH_ADD_INT(ptr_volume->particle, field_id, particle); 
  
  ptr_volume->last_used_id = ptr_volume->last_used_id + 1;
  
  printf("\n\n\t\t*** COPY PARTICLE (%i, %i) - PDF(%f) TO VOLUME (%i) WITH NEW ID (%i) ***\n\n", ptr_particle->index.ix, ptr_particle->field_id, ptr_particle->pdf, particle->index.ix, particle->field_id);
  
  particle = NULL;
  free(particle);  
  
}



/*
void print_vtk(int n_eulerian, int k) {
      struct lag_eulerian *lag_ptr;

     char string[500];
     int file_number = 0;
    
     sprintf( string, "lag_DPM_paraview_%i.vtk", k );
    

     FILE *pvtk;
     pvtk=fopen(string, "w");
     fprintf(pvtk,"# vtk DataFile Version 2.0\nPoints\nASCII\nDATASET POLYDATA\n");
     fprintf(pvtk,"Points	%i float\n",n_eulerian);
                   
     
      for(lag_ptr=id_table; lag_ptr != NULL; lag_ptr=(struct lag_eulerian*)(lag_ptr->hh.next)) {
        fprintf(pvtk,"%6.16lf	%6.16lf	%6.16lf\n",lag_ptr->x1,  lag_ptr->y1, lag_ptr->z1);
      }

     fclose(pvtk);
     
}
*/



