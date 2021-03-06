#ifndef lag_mapH
#define lag_mapH

#include "uthash.h"

typedef struct id_np {
  int id;
  UT_hash_handle hh;
} id_np_t;

typedef struct index {
  int igx, igy, igz;  
} index_t;

typedef struct volume {
  int rl;
  int iug;
  index_t index;
  struct volume *rvolume;
  struct id_np *id_np;  
  UT_hash_handle hh;
} volume_t;

volume_t* SEARCH_TABLE(const int level, const double x, const double y, const double z, const double A[], const int levels);

void INITIALIZE_MAP(const double A[], const int B[], const int size_B);

volume_t* SEARCH_PARTICLE(const double A[], const double x, const double y, const double z, int levels);


#endif