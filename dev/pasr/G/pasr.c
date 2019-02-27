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
#define targ_mix1 0.1
#define targ_mix2 0.9
#define mixlgth 0.01
#define Comega 2.0

#define FINAL_TIME 3000
#define dte 0.1

#define a1 0.0
#define a2 1.0
#define b1 0.0
#define b2 1.0
#define c1 0.0
#define c2 1.0

#define levels 1

#define nx 5
#define ny 5
#define nz 5

#define dx (a2-a1)/nx
#define dy (b2-b1)/ny
#define dz (c2-c1)/nz

#define nvol nx*ny*nz

#define np_vol_ini 100

#define tol_np 1.1
#define np_total np_vol_ini*nvol*tol_np  

#define U 1.0
#define V 0.0
#define W 0.0

#define nfiles 30
#define n_interv 100

typedef struct notionalParticle {
  
  double posix, posiy, posiz, posixo, posiyo, posizo;
  double dwall, rhop, omegap, weight, treact;
  double zetap[n_Yi];
  index_t index;
  index_t inew;
  
} notionalParticle_t;




void INITIALIZE_PARTICLES(notionalParticle_t *, const double []);

void INITIALIZE_ZETA(notionalParticle_t *);

void MOVE_MICROMIXING(notionalParticle_t [], const double [], const int, int *);
  
void PRINT_VTK(notionalParticle_t *, double);

void PDF_PRINT(notionalParticle_t *, double);

float ran1(int *);

float ran2(int *);

float gasdev(int *idum);

int timeval_subtract(struct timeval *, struct timeval *, struct timeval *);

void data_distribution(char, int [], int, double);

void PRINT_GRID(); 





  
int main() {
     
  printf("\n\n\t\t*** START ***\n\n");
    
  
//   char file_name3[20];     
//   FILE *ptr_file3; 
//    sprintf(file_name3, "PWay(%i).dat", (int)(np_vol_ini*nvol/2));  
//   ptr_file3 = fopen(file_name3, "w");  
//   fprintf(ptr_file3, "Time \t X \t Y \t Z \n");
  
  
  double stime = 0.0;
  
  char in[1];
  
  int i = 0;
  
  int idum3 = -3213853;
  int *ptr_idum3 = &idum3;
    
  struct timeval tvBegin, tvEnd, tvDiff;
   
  gettimeofday(&tvBegin, NULL);
  
  volume_t *iter_volume, *tmp_volume, *svolume;
  id_np_t *iter_particle, *tmp_particle;
  
  
//   volume_t svolume;
    
  PRINT_GRID();
  
  gets(in);
  
  const int size_B = levels*8+1;
  
  const int B[9] = {levels, 0, 1, 1, 1, nx, ny, nz, 1};
  const double A[9] = {a1, a2, b1, b2, c1, c2, dx, dy, dz};
  
//   const int B[25] = {levels, 0, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 3, 2, 1, 2, 5, 2, 1, 1, 5, 2, 1, 2, 9};
//   const double A[15] = {a1, a2, b1, b2, c1, c2, 0.5, 0.25, 0.125, 1.0, 1.0, 1.0, 0.5, 0.25, 0.125};
//   
  INITIALIZE_MAP(A, B, size_B);
  
//   printf("\n INITIALIZE_MAP(A, B, size_B);\n");
//   gets(in);
  
//   HASH_ITER(hh, eulerian, iter_volume, tmp_volume)
//     printf("\n iter_volume %i %i %i \n", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz);

//   gets(in); 
  
  notionalParticle_t *notionalParticleArray;
  
  notionalParticleArray = malloc(((int)np_total)*sizeof(notionalParticle_t));
  if(notionalParticleArray == NULL) {
    printf("\nout of memory\n");
    exit;
  }   
    
  INITIALIZE_PARTICLES(notionalParticleArray, A);
    
//   printf("\n INITIALIZE_PARTICLES(notionalParticleArray, A);\n");
//   gets(in);
      
//   svolume = SEARCH_PARTICLE(A, 0.01, 0.5, 0.501, levels);
//   printf("\n  SEARCH_PARTICLE(0.01,0.5,0.501) = %i %i %i - %i\n", svolume->index.igx, svolume->index.igy, svolume->index.igz, svolume->rl);

  
//   gets(in);

  
  INITIALIZE_ZETA(notionalParticleArray);
    
//   printf("\n INITIALIZE_ZETA(notionalParticleArray);\n");
//   gets(in);
  
//   PRINT_VTK(notionalParticleArray, 0.0);  
  
//   gets(in);
  
//   PDF_PRINT(notionalParticleArray, 0.0);    
//     
//   
//   gets(in);
  
  
    while(stime < FINAL_TIME) {

      i++;      
      stime += dte;
      
//       printf("\n\n\t\t*** TIME %f***\n\n", stime);

      MOVE_MICROMIXING(notionalParticleArray, A, levels, ptr_idum3);

      if(i % (int)(FINAL_TIME/(dte*nfiles)) == 0) {

	PDF_PRINT(notionalParticleArray, stime);

// 	PRINT_VTK(notionalParticleArray, stime); 
//       }
// 	gets(in);


//        HASH_ITER(hh, eulerian, iter_volume, tmp_volume) { 
// 	  HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) { 
// 	    if(iter_particle->id == (int)(np_vol_ini*nvol/2)) {
// 	      fprintf(ptr_file3, "%f \t %f \t %f \t %f \n", stime, notionalParticleArray[iter_particle->id].posix, notionalParticleArray[iter_particle->id].posiy, notionalParticleArray[iter_particle->id].posiz);  
// 	    }      
// 	  }
       }
    }
    
//     fclose(ptr_file3);  

    

  
  gettimeofday(&tvEnd, NULL);  
  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
  printf("\n\n\t\tTime, in seconds, for PaSR with NPV (%i) = ", np_vol_ini);
  printf("%ld.%06ld\n\n", tvDiff.tv_sec, tvDiff.tv_usec);
 
  printf("\n\n\t\t*** END ***\n\n");
  
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
  
  char in[1];
  
  int idum3 = -3213853;
  int *ptr_idum3 = &idum3;

  int m, c;
  
  int part_counter = 0;
  
  double x0, y0, z0;
      
  const int intervals = 10;
  
  const double delta_print = (a2-a1)/intervals;
  
//   int counter_x[intervals];
//   for(c=0; c<intervals; c++) counter_x[c] = 0;
//   int counter_y[intervals];
//   for(c=0; c<intervals; c++) counter_y[c] = 0;
//   int counter_z[intervals];
//   for(c=0; c<intervals; c++) counter_z[c] = 0;
 
  volume_t *iter_volume, *tmp_volume;
  
  id_np_t *ptr_particle, *iter_particle, *tmp_particle;
  
  HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {     
    
    x0 = (iter_volume->index.igx-1)*A[6];
    y0 = (iter_volume->index.igy-1)*A[6+levels];
    z0 = (iter_volume->index.igz-1)*A[6+2*levels];
        
//     printf("\n\t\tx0 = %f", x0);
//     printf("\n\t\ty0 = %f", y0);
//     printf("\n\t\tz0 = %f", z0);    

    for(m=1; m<=np_vol_ini; m++) {
      
      part_counter++;
      
//       printf("\n\t\tm = %i of np_vol_ini = %i", m, np_vol_ini);
      
      notionalParticleArray[part_counter].index.igx = iter_volume->index.igx;
      notionalParticleArray[part_counter].index.igy = iter_volume->index.igy;
      notionalParticleArray[part_counter].index.igz = iter_volume->index.igz;
   
      notionalParticleArray[part_counter].inew.igx = notionalParticleArray[part_counter].index.igx;
      notionalParticleArray[part_counter].inew.igy = notionalParticleArray[part_counter].index.igy;
      notionalParticleArray[part_counter].inew.igz = notionalParticleArray[part_counter].index.igz;
            
      notionalParticleArray[part_counter].posixo = x0 + ran2(ptr_idum3)*(A[6]);
      notionalParticleArray[part_counter].posiyo = y0 + ran2(ptr_idum3)*(A[6+levels]);
      notionalParticleArray[part_counter].posizo = z0 + ran2(ptr_idum3)*(A[6+2*levels]);
      
//       printf("\n\t\tposixo = %f", notionalParticleArray[part_counter].posixo);
//       printf("\n\t\tposiyo = %f", notionalParticleArray[part_counter].posiyo);
//       printf("\n\t\tposizo = %f", notionalParticleArray[part_counter].posizo);
//       
      notionalParticleArray[part_counter].posix = notionalParticleArray[part_counter].posixo;
      notionalParticleArray[part_counter].posiy = notionalParticleArray[part_counter].posiyo;
      notionalParticleArray[part_counter].posiz = notionalParticleArray[part_counter].posizo;
                  
      notionalParticleArray[part_counter].treact = zeros;
      notionalParticleArray[part_counter].dwall = zeros;
      notionalParticleArray[part_counter].omegap = zeros;
      notionalParticleArray[part_counter].rhop = rho_cte;
      notionalParticleArray[part_counter].weight = (rho_cte*dx*dy*dz)/np_vol_ini; //Atualizar temporalmente em função de np_vol se alterar
      
      ptr_particle = malloc(sizeof(id_np_t)); 
      ptr_particle->id = part_counter;      
      HASH_ADD_INT(iter_volume->id_np, id, ptr_particle); 
         
//       printf("\n\n\t\t*** ADD PARTICLE (%i) XYZ (%f, %f, %f) TO VOLUME (%i, %i, %i) ***\n\n", ptr_particle->id, notionalParticleArray[part_counter].posixo, notionalParticleArray[part_counter].posiyo, notionalParticleArray[part_counter].posixo, notionalParticleArray[part_counter].index.igx, notionalParticleArray[part_counter].index.igy, notionalParticleArray[part_counter].index.igz);
    
      ptr_particle = NULL;
      free(ptr_particle); 
            
    }   
  }
    
//   HASH_ITER(hh, eulerian, iter_volume, tmp_volume) { 
//     HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) { 
//     
//       for(c=0; c<intervals; c++) {
// 	if(notionalParticleArray[iter_particle->id].posixo <= delta_print*(c+1)) {      
// 	  counter_x[c]++;
// 	  break;
// 	}
//       }  
// 	  
//       for(c=0; c<intervals; c++) {
// 	if(notionalParticleArray[iter_particle->id].posiyo <= delta_print*(c+1)) {	      
// 	  counter_y[c]++;
// 	  break;
// 	}
//       } 
// 	  
//       for(c=0; c<intervals; c++) {
// 	if(notionalParticleArray[iter_particle->id].posizo <= delta_print*(c+1)) {	      
// 	  counter_z[c]++;
// 	  break;
// 	}
//       } 
//       
//     }
//   }
  
//   data_distribution('x', counter_x, intervals, delta_print);
//   data_distribution('y', counter_y, intervals, delta_print);
//   data_distribution('z', counter_z, intervals, delta_print);
  
}

     


void INITIALIZE_ZETA(notionalParticle_t *notionalParticleArray) {
    
  int c;
  
  volume_t *iter_volume, *tmp_volume;
  
  id_np_t *iter_particle, *tmp_particle;
    
  HASH_ITER(hh, eulerian, iter_volume, tmp_volume) { 
    HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) { 
      for(c=0; c<n_Yi; c++) 
	notionalParticleArray[iter_particle->id].zetap[c] = (mid*(zeta2_mix+zeta1_mix) + mid*(zeta2_mix-zeta1_mix) * tanh((two*(notionalParticleArray[iter_particle->id].posizo))/mixlgth - two*targ_mix1/mixlgth)) + (mid*(zeta1_mix+zeta2_mix) + mid*(zeta1_mix-zeta2_mix) * tanh((two*(notionalParticleArray[iter_particle->id].posizo))/mixlgth - two*targ_mix2/mixlgth)) - 1.0;
// 	notionalParticleArray[iter_particle->id].zetap[c] = mid*(zeta1_mix+zeta2_mix) + mid*(zeta1_mix-zeta2_mix) * tanh((two*(notionalParticleArray[iter_particle->id].posizo))/mixlgth - two*targ_mix/mixlgth);
    }
  }
  
}

void MOVE_MICROMIXING(notionalParticle_t notionalParticleArray[], const double A[], const int ls, int *ptr_idum3) {

  int part_counter;
  
  double delta, zetap_S_mean, weight_sum, omegam, grad_gx, grad_gy, grad_gz;
  
  char in[1];
  
  volume_t *iter_volume, *tmp_volume, *move_volume;
  
  volume_t *svolume;
  
  id_np_t *iter_particle, *tmp_particle, *tmp2_particle, *move_particle;
 
//   printf("\n MICRO \n");
  
  HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {
    
//       printf("\n iter_volume %i %i %i - %i\n", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, iter_volume->rl);
//       gets(in);
      
      part_counter = HASH_COUNT(iter_volume->id_np);
//       printf("\n part_counter %i\n", part_counter);
//       gets(in);

      delta = zeros;  
      zetap_S_mean = zeros; 
      weight_sum = zeros;
      omegam = zeros;

      HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) {
	
// 	printf("\n iter_particle %i\n", iter_particle->id);
// 	gets(in);
	
	notionalParticleArray[iter_particle->id].weight = (rho_cte*dx*dy*dz)/part_counter; 
	zetap_S_mean += (notionalParticleArray[iter_particle->id].zetap[0])*(notionalParticleArray[iter_particle->id].weight);
	weight_sum += notionalParticleArray[iter_particle->id].weight;

	grad_gx = (1/Sc)*(ni_cte)/A[6];
	grad_gy = (1/Sc)*(ni_cte)/A[6+levels];
	grad_gz = (1/Sc)*(ni_cte)/A[6+2*levels];

	notionalParticleArray[iter_particle->id].posixo = notionalParticleArray[iter_particle->id].posix;
	notionalParticleArray[iter_particle->id].posiyo = notionalParticleArray[iter_particle->id].posiy;
	notionalParticleArray[iter_particle->id].posizo = notionalParticleArray[iter_particle->id].posiz;
	notionalParticleArray[iter_particle->id].posix = notionalParticleArray[iter_particle->id].posixo + U*dte + 0.0*dte*grad_gx + sqrt(dte*(two/Sc)*(ni_cte))*gasdev(ptr_idum3);
	notionalParticleArray[iter_particle->id].posiy = notionalParticleArray[iter_particle->id].posiyo + V*dte + 0.0*dte*grad_gy + sqrt(dte*(two/Sc)*(ni_cte))*gasdev(ptr_idum3);      
	notionalParticleArray[iter_particle->id].posiz = notionalParticleArray[iter_particle->id].posizo + W*dte + 0.0*dte*grad_gz + sqrt(dte*(two/Sc)*(ni_cte))*gasdev(ptr_idum3);
	
// 	if(notionalParticleArray[iter_particle->id].posix < a1)
// 	  notionalParticleArray[iter_particle->id].posix = a1 + (a1 - notionalParticleArray[iter_particle->id].posix);
// 	if(notionalParticleArray[iter_particle->id].posix > a2)
// 	  notionalParticleArray[iter_particle->id].posix = a2 - (notionalParticleArray[iter_particle->id].posix - a2);

	if(notionalParticleArray[iter_particle->id].posiy < b1)
	  notionalParticleArray[iter_particle->id].posiy = b1 + (b1 - notionalParticleArray[iter_particle->id].posiy);
	if(notionalParticleArray[iter_particle->id].posiy > b2)
	  notionalParticleArray[iter_particle->id].posiy = b2 - (notionalParticleArray[iter_particle->id].posiy - b2);
      
	if(notionalParticleArray[iter_particle->id].posiz < c1)
	  notionalParticleArray[iter_particle->id].posiz = c1 + (c1 - notionalParticleArray[iter_particle->id].posiz);
	if(notionalParticleArray[iter_particle->id].posiz > c2)
	  notionalParticleArray[iter_particle->id].posiz = c2 - (notionalParticleArray[iter_particle->id].posiz - c2);

	//Periodic Boundary Conditions in X direction - U positive
	if(notionalParticleArray[iter_particle->id].posix > a2) 
	  notionalParticleArray[iter_particle->id].posix = notionalParticleArray[iter_particle->id].posix - (a2-a1);
	  
	svolume = SEARCH_PARTICLE(A, notionalParticleArray[iter_particle->id].posix, notionalParticleArray[iter_particle->id].posiy, notionalParticleArray[iter_particle->id].posiz, ls);
	notionalParticleArray[iter_particle->id].inew.igx = svolume->index.igx;
	notionalParticleArray[iter_particle->id].inew.igy = svolume->index.igy;
	notionalParticleArray[iter_particle->id].inew.igz = svolume->index.igz;
      }      
      zetap_S_mean = zetap_S_mean/weight_sum;
      HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) {
	delta = pow(A[6]*A[6+levels]*A[6+2*levels], 1/3);
	omegam = Comega*(ni_cte/Sc)/pow(delta, 2.0);
	notionalParticleArray[iter_particle->id].zetap[0] = notionalParticleArray[iter_particle->id].zetap[0] - omegam*dte*(notionalParticleArray[iter_particle->id].zetap[0] - zetap_S_mean);
      }         
  }  
  
//   printf("\n MOVE \n");
//   gets(in);

  
  HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {
    
//     printf("\n iter_volume %i %i %i - %i\n", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, iter_volume->rl);
//     gets(in);
    
       
//  for(iter_particle=iter_volume->id_np; iter_particle!=NULL; iter_particle=iter_particle->hh.next) {

    iter_particle=iter_volume->id_np;
   
    while(iter_particle != NULL) {
      
//       printf("\n iter_particle %i\n", iter_particle->id);
//       gets(in);
      
      if( (notionalParticleArray[iter_particle->id].inew.igx != notionalParticleArray[iter_particle->id].index.igx) || (notionalParticleArray[iter_particle->id].inew.igy != notionalParticleArray[iter_particle->id].index.igy) || (notionalParticleArray[iter_particle->id].inew.igz != notionalParticleArray[iter_particle->id].index.igz) ) {
	
	notionalParticleArray[iter_particle->id].index.igx = notionalParticleArray[iter_particle->id].inew.igx;
	notionalParticleArray[iter_particle->id].index.igy = notionalParticleArray[iter_particle->id].inew.igy;
	notionalParticleArray[iter_particle->id].index.igz = notionalParticleArray[iter_particle->id].inew.igz;
	
// 	printf("\n particle must be moved \n");
// 	gets(in);

	tmp_particle = malloc(sizeof(id_np_t));
	tmp_particle->id = iter_particle->id;
	tmp2_particle = iter_particle;

	
// 	printf("\n tmp_particle %i\n", tmp_particle->id);
// 	gets(in);
	
	
	iter_particle=iter_particle->hh.next;
// 	printf("\n iter_particle->hh.next %i \n", iter_particle->id);

// 	move_particle = malloc(sizeof(id_np_t));
// 	move_particle->id = iter_particle->id;
// 	move_volume = malloc(sizeof(volume_t));
	svolume = SEARCH_PARTICLE(A, notionalParticleArray[tmp_particle->id].posix, notionalParticleArray[tmp_particle->id].posiy, notionalParticleArray[tmp_particle->id].posiz, ls);
	move_volume = svolume; //RETIRAR DEPOIS, SVOLUME é PONTEIRO AGORA, PODE SER USADO NO LUGAR DO MOVE_VOLUME ou VICE VERSA
	
// 	printf("\n TO move_volume %i %i %i - %i\n", move_volume->index.igx, move_volume->index.igy, move_volume->index.igz, move_volume->rl);
// 	gets(in);
	
	HASH_ADD_INT(move_volume->id_np, id, tmp_particle);
	
	
// 	move_particle = NULL;
// 	free(move_particle);
// 	move_volume = NULL;
// 	free(move_volume);
// 	notionalParticleArray[iter_particle->id].index = notionalParticleArray[iter_particle->id].inew;
// 	tmp_particle = malloc(sizeof(id_np_t));
// 	tmp_particle = iter_particle;
	

	
	HASH_DEL(iter_volume->id_np, tmp2_particle);
	free(tmp2_particle);
      }   
      
      else
	iter_particle=iter_particle->hh.next;

      
    }  
  }  
}


void PDF_PRINT(notionalParticle_t *notionalParticleArray, double stime) {
  
  int i, j, k;
  
  int part_counter;
  
  double weight_sum;
  
  char in[1];
    
  volume_t *iter_volume, *tmp_volume;
  
  id_np_t *iter_particle, *tmp_particle;
  
  char file_name[100]; 
  FILE *ptr_file; 
  
  char file_name2[100]; 
  FILE *ptr_file2; 
  
//   char file_name3[100]; 
//   FILE *ptr_file3; 
  
  double cdf[n_interv];    
  double pdf[n_interv];
  
  double interv[n_interv+1];
  double d_in = (zeta2_mix-zeta1_mix)/n_interv;
    
  interv[0] = zeta1_mix;
  for(i=1; i<(n_interv+1); i++)     
    interv[i] = interv[i-1] + d_in;
  
  sprintf(file_name2, "PV(%1.2f).dat", stime);  
  ptr_file2 = fopen(file_name2, "w");  
  fprintf(ptr_file2, "VOLUME \t NP \n");
  
//Peso não pode ser medido aqui porque carrega o peso de cada partícula do tempo anterior, porém com a quantidade de partículas do tempo atual
//   sprintf(file_name3, "WEIGHT(%1.2f).dat", stime);  
//   ptr_file3 = fopen(file_name3, "w");  
//   fprintf(ptr_file3, "VOLUME \t WEIGHT \n");

  HASH_ITER(hh, eulerian, iter_volume, tmp_volume) {     

    part_counter = HASH_COUNT(iter_volume->id_np);  
    
    sprintf(file_name, "PDF_%i%i%i(%1.2f).dat", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, stime);  
    ptr_file = fopen(file_name, "w");  
    fprintf(ptr_file, "NORM \t PDF \n"); 
        
    weight_sum = zeros;
    

    for(i=0; i<n_interv; i++) {
      cdf[i] = 0.0;
      pdf[i] = 0.0;
    }

    HASH_ITER(hh, iter_volume->id_np, iter_particle, tmp_particle) {
        
      weight_sum += notionalParticleArray[iter_particle->id].weight;
      
      for(i=0; i<n_interv; i++) {
	if(notionalParticleArray[iter_particle->id].zetap[0] < interv[i+1]) {
	  cdf[i]++;  
	  break;
	}
      }
      
    }
    
    for(i=0; i<n_interv; i++) {
      pdf[i] = cdf[i]/part_counter;
      fprintf(ptr_file, "%4.6f \t %4.6f \n", mid*(interv[i]+interv[i+1]), pdf[i]); 
    }
    
    
//     printf("\nPRINT %i-%i-%i \t %i %f \n", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, part_counter, weight_sum);   

    fprintf(ptr_file2, "%i-%i-%i \t %i \n", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, part_counter);  
//     fprintf(ptr_file3, "%i-%i-%i \t %f \n", iter_volume->index.igx, iter_volume->index.igy, iter_volume->index.igz, weight_sum);   
    
    
    fclose(ptr_file);  
        
  }
  
  fclose(ptr_file2);
//   fclose(ptr_file3);
    
}


  

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

  

void PRINT_VTK(notionalParticle_t *notionalParticleArray, double k) {
  
  int i = 0;
    
  char string[40];
  int file_number = 0;
   
  sprintf(string, "particles_%2.4f.vtk", k);
    
  FILE *pvtk;
  pvtk = fopen(string, "w");
  fprintf(pvtk, "# vtk DataFile Version 2.0\nPoints\nASCII\nDATASET POLYDATA\n");
  fprintf(pvtk, "Points %i float\n", np_vol_ini*nvol);
          
  for(i=1; i<=np_vol_ini*nvol; i++) 
    fprintf(pvtk, "%1.2e %1.2e %1.2e \n", notionalParticleArray[i].posix, notionalParticleArray[i].posiy, notionalParticleArray[i].posiz);
      
  fprintf(pvtk, "POINT_DATA %d\n", np_vol_ini*nvol);
  fprintf(pvtk, "SCALARS zeta float\n");
  fprintf(pvtk, "LOOKUP_TABLE default\n");
    
  for(i=1; i<=np_vol_ini*nvol; i++) 
    fprintf(pvtk, "%f \n", notionalParticleArray[i].zetap[0]);
  
  fclose(pvtk);
     
}


void PRINT_GRID() {
   
  int i, j, k;
  
  double x, y, z;
  
  char string[10];
  int file_number = 0;
   
  sprintf(string, "grid.vtk");
    
  FILE *pvtk;
  pvtk = fopen(string, "w");
  fprintf(pvtk, "# vtk DataFile Version 3.0\nStructured 3D Dataset\nASCII\nDATASET STRUCTURED_GRID\n");
  fprintf(pvtk, "DIMENSIONS %d %d %d\n", nx+1, ny+1, nz+1);
  fprintf(pvtk, "POINTS %d float\n", (nx+1)*(ny+1)*(nz+1));
  
  x = 0;
  y = 0;
  z = 0;
  
  for(i=1; i<=nx+1; i++)
    for(j=1; j<=ny+1; j++)
      for(k=1; k<=nz+1; k++) {
	x = (i-1)*dx;
	y = (j-1)*dy;
	z = (k-1)*dz;
	fprintf(pvtk, "%f %f %f\n", x, y, z);
      }
  
  
   
  fclose(pvtk);
     
}


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


