#ifndef HandleStochParticleH
#define HandleStochParticleH

#include "uthash.h"


#define IA 16807
#define IM 2147483647
#define AM_RanOne (1.0/IM)
#define IQ 127773
#define IR 2836
#define NDIV_RanOne (1+(IM-1)/NTAB)

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 64
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define n_Yi 1

#define ITENS 6

#define FLAG1 1
#define FLAG2 2

#define max(a, b) ((a)>(b)?(a):(b))
#define min(a, b) ((a)<(b)?(a):(b))

#define bool int
#define true 1
#define false 0

#define n_Yi 1
#define rho_cte 1.0

#define zeta1_mix 0.0
#define zeta2_mix 1.0
#define targ_mix 0.5
#define mixlgth 0.1
#define Comega 2.0

#define U 1.0
#define V 0.0
#define W 0.0 

#define Sc 0.7
#define Ni 1E-5

#define FINAL_TIME 1.0
#define dt 0.1

#define ga1 0.0
#define gb1 1.0
#define ga2 0.0
#define gb2 1.0
#define ga3 0.0
#define gb3 1.0

#define npl 2

#define np_total 10

#define nfiles 10
#define n_interv 100


//Structures definition

typedef struct notionalParticle {  
  double posix, posiy, posiz, posixo, posiyo, posizo;
  double dwall, rhop, omegap, weight, treact;
  double zetap[n_Yi];
  int id;
  int proc, new_proc;
  index_t index, new_index;
  UT_hash_handle hh;  
} notionalParticle_t;


typedef struct idStack {  
  unsigned long size;
  unsigned long capacity;
  unsigned long *id;    
} idStack_t;





//Headers of transport functions

void Move_XYZ(notionalParticle_t * const iter_particle, const double A[]);
void Update_New_Index(const double A[], notionalParticle_t * const iter_particle);
void Map_ID(const double A[], notionalParticle_t * const iter_particle);
void Update_Particle_Position(notionalParticle_t * const iter_particle);
void U_Particles_Init(const int my_rank, const double a1, const double a2, const double a3, const double b1, const double b2, const double b3, const double A[]);
void Particles_Zeta_Init(const int);
void Transport_Stoch_Particles(const double A[], int *ct_seed, const int my_rank);
void Bndy_Box(const double A[], notionalParticle_t * const iter_particle);


//Random functions' headers

float RanOne(int *);
float RanTwo(long int *);
float Gasdev(int *);


//Stack functions' headers

idStack_t* createStack(unsigned int max_np_proc);
unsigned int top(idStack_t *S);
void pop(idStack_t *S);
void push(idStack_t *S, unsigned int id);


//Miscelaneous headers

void Print_Particles_VTK(const int, const int); 
int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1);


#endif