#ifndef lag_mapH
#define lag_mapH

#include "uthash.h"

typedef struct id_particle {
  int id;
  UT_hash_handle hh;
} id_t;

typedef struct index {
  int igx, igy, igz;  
} index_t;

typedef struct volume {
  int rl;
  int iug;
  index_t index;
  struct volume *rvolume;
  struct id_particle *id_particle;  
  UT_hash_handle hh;
} volume_t;

volume_t* Search_Table(const int level, const double x, const double y, const double z, const double A[], const int levels);

void Map_Init(const double A[], const int B[], const int size_B);

volume_t* Search_Particle(const double A[], const double x, const double y, const double z, int levels);


#endif