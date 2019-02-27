#include <stdio.h>			/* printf */
#include <stdlib.h>			/* malloc */
#include <sys/time.h>			/* timeval */
#include "uthash.h"
#include "map.h"


#define IA 16807
#define IM 2147483647
#define AM_RAN1 (1.0/IM)
#define IQ 127773
#define IR 2836
#define NDIV_RAN1 (1+(IM-1)/NTAB)

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

#define zeros 0.0
#define mid 0.5
#define one 1.0
#define two 2.0

#define max(a, b) ((a)>(b)?(a):(b))
#define min(a, b) ((a)<(b)?(a):(b))

#define bool int
#define true 1
#define false 0

#define start_up false

#define RANGE_ID np_vol_ini*2 (fator of 2 means 50% of probability to find a new available id)

#define n_Yi 1
#define rho_cte 1.0
#define ni_cte 1.0E-3
#define Sc 0.7
#define zeta1_mix 0.0
#define zeta2_mix 1.0
#define targ_mix 0.5
#define mixlgth 0.1
#define Comega 2.0


#define a1 0.0
#define a2 1.0
#define b1 0.0
#define b2 1.0
#define c1 0.0
#define c2 1.0

#define levels 1

#define nx 2
#define ny 2
#define nz 2

#define dx (a2-a1)/nx
#define dy (b2-b1)/ny
#define dz (c2-c1)/nz

#define nvol nx*ny*nz
#define np_vol_ini 100
#define tol_np 1.1
#define np_total np_vol_ini*nvol*tol_np

#define n_interv 100

#define xtotc 1
#define ytotc 1
#define ztotc 1

#define FINAL_TIME 1000.0
#define dte 0.1

#define U 0.0
#define V 0.0
#define W 0.0


typedef struct notionalParticle {
  
  double posix, posiy, posiz, posixo, posiyo, posizo;
  double dwall, rhop, omegap, weight, treact;
  double zetap[n_Yi];
  index_t index;  
  
} notionalParticle_t;




void INITIALIZE_PARTICLES(notionalParticle_t *, const double []);

// void INITIALIZE_ZETA();
// 
// void MOVE_MICROMIXING(int *ptr_idum3);
// 
// void PRINT_VTK(double);
// 
// void PDF_PRINT(double);
// 
// float ran1(*/int *);

float ran2(int *);

// float gasdev(int *idum);

int timeval_subtract(struct timeval *, struct timeval *, struct timeval *);

void data_distribution(char, int [], int, double);








  
int main() {
     
  printf("\n\n\t\t*** START ***\n\n");
  
  double stime = 0.0;
  
  int i = 0;
    
  struct timeval tvBegin, tvEnd, tvDiff;
   
  gettimeofday(&tvBegin, NULL);
    
  
  const int size_B = levels*8+1;
  const int B[9] = {levels, 0, 1, 1, 1, nx, ny, nz, 1};
  const double A[9] = {a1, a2, b1, b2, c1, c2, dx, dy, dz};
  
  INITIALIZE_MAP(A, B, size_B);

  notionalParticle_t *notionalParticleArray;
  
  notionalParticleArray = malloc(np_total*sizeof(notionalParticle_t));
  if(notionalParticleArray == NULL) {
    printf("\nout of memory\n");
    exit;
  }   
    
  INITIALIZE_PARTICLES(notionalParticleArray, A);
    
//   INITIALIZE_ZETA();
    
    //PRINT_VTK(0.0);
    
    //PDF_PRINT(0.0);
            
    
    /*
    while(stime <= FINAL_TIME) {
      i++;     
      stime += dte;
      MOVE_MICROMIXING(ptr_idum3);      
      //if(i % 1000 == 0) {
	//PDF_PRINT(stime);
	//PRINT_VTK(stime);      
      //}
    }*/
    
//   }
  
  gettimeofday(&tvEnd, NULL);  
  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
  printf("\n\n\t\tTime, in seconds, for PDF with NPV (%i) = ", np_vol_ini);
  printf("%ld.%06ld\n\n", tvDiff.tv_sec, tvDiff.tv_usec);
 
  //printf("\n\n\t\t*** END ***\n\n");
  
  return 0;
  
}





int timeval_subtract(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}



void INITIALIZE_PARTICLES(notionalParticle_t *notionalParticleArray, const double A[]) {
  
  char in;
  
  int idum3 = -3213853;
  int *ptr_idum3 = &idum3;

  int m, c;
  
  double x0, y0, z0;
      
//   const int intervals = 10;
  
//   const double delta_print = (a2-a1)/intervals;
  
//   int counter_x[intervals];
//   for(c=0; c<intervals; c++) counter_x[c] = 0;
//   int counter_y[intervals];
//   for(c=0; c<intervals; c++) counter_y[c] = 0;
//   int counter_z[intervals];
//   for(c=0; c<intervals; c++) counter_z[c] = 0;
 
  volume_t *iter_volume, *tmp_volume;
  
  id_np_t *ptr_particle;
  
  HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {     
    
    x0 = (iter_volume->index.igx-1)*A[6];
    y0 = (iter_volume->index.igy-1)*A[6+levels];
    z0 = (iter_volume->index.igz-1)*A[6+2*levels];
        
    for(m=1; m<=np_vol_ini; m++) {
      
      printf("\n\t\tm = %i of np_vol_ini = %i", m, np_vol_ini);
      
      notionalParticleArray[m].index.igx = iter_volume->index.igx;
      notionalParticleArray[m].index.igy = iter_volume->index.igy;
      notionalParticleArray[m].index.igz = iter_volume->index.igz;
   
      notionalParticleArray[m].posixo = x0 + ran2(ptr_idum3)*(A[6]);
      notionalParticleArray[m].posiyo = y0 + ran2(ptr_idum3)*(A[6+levels]);
      notionalParticleArray[m].posizo = z0 + ran2(ptr_idum3)*(A[6+2*levels]);
      
      notionalParticleArray[m].posix = notionalParticleArray[m].posixo;
      notionalParticleArray[m].posiy = notionalParticleArray[m].posiyo;
      notionalParticleArray[m].posiz = notionalParticleArray[m].posizo;
                  
      notionalParticleArray[m].treact = zeros;
      notionalParticleArray[m].dwall = zeros;
      notionalParticleArray[m].omegap = zeros;
      notionalParticleArray[m].rhop = rho_cte;
      notionalParticleArray[m].weight = (rho_cte*dx*dy*dz)/np_vol_ini;
      
      ptr_particle = malloc(sizeof(id_np_t)); 
      ptr_particle->id = m;      
      HASH_ADD_INT(iter_volume->id_np, id, ptr_particle); 
      
     
      
//       // COUNTER
//       for(c=0; c<intervals; c++) {
// 	if(ptr_particle->posixo <= delta_print*(c+1)) {	      
// 	  counter_x[c]++;
// 	  break;
// 	}
//       }  
// 	  
//       for(c=0; c<intervals; c++) {
// 	if(ptr_particle->posiyo <= delta_print*(c+1)) {	      
// 	  counter_y[c]++;
// 	  break;
// 	}
//       } 
// 	  
//       for(c=0; c<intervals; c++) {
// 	if(ptr_particle->posizo <= delta_print*(c+1)) {	      
// 	  counter_z[c]++;
// 	  break;
// 	}
//       } 
	  
     
      printf("\n\n\t\t*** ADD PARTICLE (%i) XYZ (%f, %f, %f) TO VOLUME (%i, %i, %i) ***\n\n", ptr_particle->id, notionalParticleArray[m].posixo, notionalParticleArray[m].posiyo, notionalParticleArray[m].posixo, notionalParticleArray[m].index.igx, notionalParticleArray[m].index.igy, notionalParticleArray[m].index.igz);
    
      ptr_particle = NULL;
      free(ptr_particle); 
            
    }    
    
//     data_distribution('x', counter_x, intervals, delta_print);
//     data_distribution('y', counter_y, intervals, delta_print);
//     data_distribution('z', counter_z, intervals, delta_print);
    
  }
}

     
  
/*

void INITIALIZE_ZETA() {
    
  int c;
  
  volume_t *iter_volume, *tmp_volume;
  
  partpdf_t *iter_particle, *tmp_particle;
    
  HASH_ITER(hh, eulerian, iter_volume, tmp_volume) { 
    HASH_ITER(hh, iter_volume->particle, iter_particle, tmp_particle) { 
      for(c=0; c<n_Yi; c++) 
	iter_particle->zetap[c] = mid*(zeta1_mix+zeta2_mix) + mid*(zeta1_mix-zeta2_mix) * tanh((two*iter_particle->posizo)/mixlgth - two*targ_mix/mixlgth);
    }
  }
  
}





//  Subroutine involving the transport of particles in compositional domain.
//  This subroutine deals with the transport of particles in compositional domain due to the
//  conditional molecular mixing. In general, the conditional molecular mixing term can be
//  decomposed in two contributions: (i) the first one representing the effects of molecular diffusion
//  on the PDF transport in physical space, and (ii) the second one (micro-mixing) representing the
//  effects of the scalar dissipation rate on the PDF transport in the scalar field sample space.
//  This specific subroutine deals with the micro-mixing term referred above. Two different models
//  are then make available through this subroutine, the Interaction by Exchange with the Mean (IEM)
//  model and the COalescence and DIspersion one (Curl). The IEM model involves basically first the
//  computation of the turbulent frequency and then the computation of the micro-mixing term using
//  this parameter. At the end, limits are imposed on scalar parameter values determined using this
//  model.
//
//  From Fluids (@author J. M. Vedovoto)



void MOVE_MICROMIXING(int *ptr_idum3) {
  
  double delta = zeros;  
  double zetap_S_mean = zeros;  
  double weight_sum = zeros;
  double omegam = zeros;
  
  double grad_gx, grad_gy, grad_gz;
  
  volume_t *iter_volume, *tmp_volume;
  
  partpdf_t *iter_particle, *tmp_particle;
    
  // Loop for Mean calculation - There is just one volume in PaSR - Be carefully when working with particle transport in eulerian domain
  HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {
      HASH_ITER(hh, iter_volume->particle, iter_particle, tmp_particle) {
	
	zetap_S_mean += (iter_particle->zetap[0])*(iter_particle->weight);
	weight_sum += iter_particle->weight;
	
	grad_gx = (1/Sc)*(iter_volume->ni)/(iter_volume->dx);
	grad_gy = (1/Sc)*(iter_volume->ni)/(iter_volume->dy);
	grad_gz = (1/Sc)*(iter_volume->ni)/(iter_volume->dz);
	
	iter_particle->posixo = iter_particle->posix;
	iter_particle->posiyo = iter_particle->posiy;
	iter_particle->posizo = iter_particle->posiz;
	iter_particle->posix = iter_particle->posixo + U*dte + dte*grad_gx + sqrt(dte*(two/Sc)*(iter_volume->ni))*gasdev(ptr_idum3);
	iter_particle->posiy = iter_particle->posiyo + V*dte + dte*grad_gy + sqrt(dte*(two/Sc)*(iter_volume->ni))*gasdev(ptr_idum3);      
	iter_particle->posiz = iter_particle->posizo + W*dte + dte*grad_gz + sqrt(dte*(two/Sc)*(iter_volume->ni))*gasdev(ptr_idum3);
	
	//Wall Boundary Conditions  
	
	if(iter_particle->posix < a1)
	  iter_particle->posix = a1 + (a1 - iter_particle->posix);
	if(iter_particle->posix > a2)
	  iter_particle->posix = a2 - (iter_particle->posix - a2);
      
	if(iter_particle->posiy < b1)
	  iter_particle->posiy = b1 + (b1 - iter_particle->posiy);
	if(iter_particle->posiy > b2)
	  iter_particle->posiy = b2 - (iter_particle->posiy - b2);
      
	if(iter_particle->posiz < c1)
	  iter_particle->posiz = c1 + (c1 - iter_particle->posiz);
	if(iter_particle->posiz > c2)
	  iter_particle->posiz = c2 - (iter_particle->posiz - c2);
      }
  }
    
  zetap_S_mean = zetap_S_mean/weight_sum;
  
  HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {
      HASH_ITER(hh, iter_volume->particle, iter_particle, tmp_particle) {
	delta = pow((iter_volume->dx)*(iter_volume->dy)*(iter_volume->dz), 1/3);
	// Turbulent frequency - Colluci 1998 and Fox 2005
	omegam = Comega*((iter_volume->ni)/Sc)/pow(delta, 2.0);
	// Micromixing model IEM
	iter_particle->zetap[0] = iter_particle->zetap[0] - omegam*dte*(iter_particle->zetap[0] - zetap_S_mean);	
      }
  }
  
}


void PDF_PRINT(double stime) {
  
  int i;
  
  char file_name[50];  
  sprintf(file_name, "PDF(%1.2f).dat", stime);  
  FILE *ptr_file;  
  ptr_file = fopen(file_name, "w");  
  fprintf(ptr_file, "COORDENADA (m) \t PDF \n"); 
  
  double cdf[n_interv];    
  double pdf[n_interv];
  
  for(i=0; i<n_interv; i++) {
    cdf[i] = 0.0;
    pdf[i] = 0.0;
  }
    
  
  //Nodes Between Intervals
  double interv[n_interv+1];
  double d_in = zeta2_mix/(n_interv);
  
  interv[0] = zeta1_mix;
  
  for(i=1; i<(n_interv+1); i++)     
    interv[i] = interv[i-1] + d_in;
  
  volume_t *iter_volume, *tmp_volume;
  
  partpdf_t *iter_particle, *tmp_particle;
     
  HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {
    HASH_ITER(hh, iter_volume->particle, iter_particle, tmp_particle) {
	
      for(i=0; i<n_interv; i++) {
	if(iter_particle->zetap[0] < interv[i+1]) {
	  cdf[i]++;  
	  break;
	}
      }
    }
  }
  
  // PDF INSTRUMENTATION
  for(i=0; i<n_interv; i++) {
    pdf[i] = cdf[i]/np_vol_ini;
    fprintf(ptr_file, "%4.6f \t %4.6f \n", mid*(interv[i]+interv[i+1]), pdf[i]); 
  }
  
  fclose(ptr_file);
    
}
*/

    

void data_distribution(char type, int counter[], int intervals, double delta_print) {
  
  int c;
  
  char file_name[30];
  
  sprintf(file_name, "ran2-%c-npv-%d.dat", type, np_vol_ini);
  
  FILE *file;
  
  file = fopen(file_name, "w");
  
  fprintf(file, "Mid Point \t Particles \n");
  
  for(c=0; c<intervals; c++)
    fprintf(file, "%f \t %d \n", (c+0.5)*delta_print, counter[c]);
  
  fclose(file);
  
}



//  “Minimal” random number generator of Park and Miller with Bays-Durham shuffle and added
// safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint
// values). Call with idum a negative integer to initialize; thereafter, do not alter idum between
// successive deviates in a sequence. RNMX should approximate the largest floating value that is
// less than 1. 
/*
float ran1(int *idum) {
  
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0 || !iy) { 								// Initialize.
  if (-(*idum) < 1) *idum=1; 								// Be sure to prevent idum = 0.
  else *idum = -(*idum);
  for (j=NTAB+7;j>=0;j--) { 								// Load the shuffle table (after 8 warm-ups).
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  if (j < NTAB) iv[j] = *idum;
  }
  iy=iv[0];
  }
  k=(*idum)/IQ; 									// Start here when not initializing.
  *idum=IA*(*idum-k*IQ)-IR*k; 								// Compute idum=(IA*idum) % IM without overif
  if(*idum < 0) *idum += IM; 								// flows by Schrage’s method.
  j=iy/NDIV_RAN1; 										// Will be in the range 0..NTAB-1.
  iy=iv[j]; 										// Output previously stored value and refill the
  iv[j] = *idum; 									// shuffle table.
  if ((temp=AM_RAN1*iy) > RNMX) return RNMX; 						// Because users don’t expect endpoint values.
  else return temp;
    
}



// Long period (> 2 × 1018) random number generator of L’Ecuyer with Bays-Durham shuffle
// and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
// the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
// idum between successive deviates in a sequence. RNMX should approximate the largest floating
// value that is less than 1.
*/

float ran2(int *idum) {

  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0) { 									//Initialize.
    if (-(*idum) < 1) *idum=1;								//Be sure to prevent idum = 0.
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { 								//Load the shuffle table (after 8 warm-ups).
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1; 									//Start here when not initializing.
  *idum=IA1*(*idum-k*IQ1)-k*IR1; 							//Compute idum=(IA1*idum) % IM1 without
  if (*idum < 0) *idum += IM1; 								//overflows by Schrage’s method.
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; 							//Compute idum2=(IA2*idum) % IM2 likewise.
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV; 										//Will be in the range 0..NTAB-1.
  iy=iv[j]-idum2; 									//Here idum is shuffled, idum and idum2 are
  iv[j] = *idum; 									//combined to generate output.
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;							// Because users don’t expect endpoint values.
  else return temp;

}  

// Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum)
// as the source of uniform deviates.
/*
float gasdev(int *idum) {

  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  if (*idum < 0) iset=0; 
  if (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0; 
      v2=2.0*ran1(idum)-1.0; 
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0); 
    fac=sqrt(-2.0*log(rsq)/rsq);
 gset=v1*fac;
    iset=1;
    return v2*fac;
  } 
  else { 
    iset=0; 
    return gset;    
  } 
  
}

  

void PRINT_VTK(double k) {
  
  volume_t *ptr_volume; 
  partpdf_t *iter_particle, *tmp_particle;
  index_t index;

  
  char string[40];
  int file_number = 0;
   
  sprintf(string, "particles_%2.4f.vtk", k);
    
  FILE *pvtk;
  pvtk = fopen(string, "w");
  fprintf(pvtk, "# vtk DataFile Version 2.0\nPoints\nASCII\nDATASET POLYDATA\n");
  fprintf(pvtk, "Points %i float\n", np_vol_ini);

                 
  ptr_volume = malloc(sizeof(volume_t));
  ptr_volume = malloc(sizeof(volume_t));
  
  index.indix = 1;
  index.indjy = 1;
  index.indkz = 1;
  
  HASH_FIND(hh, eulerian, &index, sizeof(index_t), ptr_volume);
  
  HASH_ITER(hh, ptr_volume->particle, iter_particle, tmp_particle)    
   fprintf(pvtk, "%1.2e %1.2e %1.2e \n", iter_particle->posix, iter_particle->posiy, iter_particle->posiz);
  
  fprintf(pvtk, "POINT_DATA %d\n", np_vol_ini);
  fprintf(pvtk, "SCALARS zeta float\n");
  fprintf(pvtk, "LOOKUP_TABLE default\n");
   
  HASH_ITER(hh, ptr_volume->particle, iter_particle, tmp_particle)    
   fprintf(pvtk, "%f \n", iter_particle->zetap[0]);

   
  fclose(pvtk);
     
}*/



// void ADD_PARTICLE(volume_t *ptr_volume, int id, double pdf) {
//        
//   partpdf_t *particle = malloc(sizeof(partpdf_t));
//   particle->id = id;
//   particle->index.ix = ptr_volume->index.ix;
//   particle->pdf = pdf;
//   
//   HASH_ADD_INT(ptr_volume->particle, id, particle); 
//   
//   ptr_volume->last_used_id = ptr_volume->last_used_id + 1;
//   
//   printf("\n\t\t*** VOLUME (%i) -> ADD_PARTICLE(%i, %i) ***\n\n", ptr_volume->index.ix, particle->index.ix, particle->id);
//   
//   particle = NULL;
//   free(particle);
//    
// }


// void COPY_PARTICLE(volume_t *ptr_volume, partpdf_t *ptr_particle, int id) {
//   
//   partpdf_t *particle = malloc(sizeof(partpdf_t));
//   particle->id = id;  
//   particle->index.ix = ptr_volume->index.ix;
//   particle->pdf = ptr_particle->pdf;
//     
//   HASH_ADD_INT(ptr_volume->particle, id, particle); 
//   
//   ptr_volume->last_used_id = ptr_volume->last_used_id + 1;
//   
//   printf("\n\n\t\t*** COPY PARTICLE (%i, %i) - PDF(%f) TO VOLUME (%i) WITH NEW ID (%i) ***\n\n", ptr_particle->index.ix, ptr_particle->id, ptr_particle->pdf, particle->index.ix, particle->id);
//   
//   particle = NULL;
//   free(particle);  
//   
// }


// tmp_particle = malloc(sizeof(partpdf_t));

//     do {
//       new_particle_id = rand()%RANGE_ID;
//       HASH_FIND_INT(toVolume->particle, &new_particle_id, tmp_particle);
//       printf("\n\t\t\t    try new_particle_id (%i) \n\n", new_particle_id);
//       gets(in);
//     }
//     while(tmp_particle != NULL);
// 
//     printf("\n\t\t\t    find new_particle_id (%i) \n\n", new_particle_id);
//     gets(in);
