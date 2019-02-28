#include <stdio.h>					/* printf */
#include <stdlib.h>					/* malloc */
#include <time.h>					/* time */
#include <sys/time.h>					/* timeval */
#include "uthash.h"
#include "map.h"

#define max(a, b) ((a)>(b)?(a):(b))
#define min(a, b) ((a)<(b)?(a):(b))

#define a1 0.0
#define a2 1.0
#define b1 0.0
#define b2 1.0
#define c1 0.0
#define c2 1.0

int timeval_subtract(struct timeval *, struct timeval *, struct timeval *);



char in[1];

int main() {
     
  printf("\n\n\t\t*** START ***\n\n");
  
  struct timeval tvBegin, tvEnd, tvDiff;
   
  gettimeofday(&tvBegin, NULL);
    
  printf("\n\n\t\t*** sizeof(vol) = %d ***\n\n", sizeof(volume_t));
  
  int levels = 7;
  
  const int size_B = levels*8+1;;
  
  const int size_A = levels*3+6;
  
  const double A[27] = {a1, a2, b1, b2, c1, c2, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125};
  
  const int B[57] = {levels, 0, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 4, 4, 4, 9, 2, 1, 1, 1, 8, 8, 8, 73, 3, 1, 1, 1, 16, 16, 16, 585, 4, 1, 1, 1, 32, 32, 32, 4681, 5, 1, 1, 1, 64, 64, 64, 37449, 6, 1, 1, 1, 128, 128, 128, 299593};
  
  
  INITIALIZE_MAP(A, B, size_B);
 
  
  volume_t volume1, volume2;
   
  double x1 = 0.49;
  double y1 = 0.49;
  double z1 = 0.51;
  
  printf("\n\n\t\t*** SEARCHING FOR PARTICLE'S POSITION x(%f), y(%f), z(%f) ***\n\n", x1, y1, z1);
  
  volume1 = SEARCH_PARTICLE(A, x1, y1, z1, levels);
  
  printf("\n\n\t\tP1 IS IN VOLUME WITH INDICES i(%d) j(%d) k(%d) OF rl(%d) \n\n", volume1.index.igx, volume1.index.igy, volume1.index.igz, volume1.rl);
      
 
  gettimeofday(&tvEnd, NULL);  
  timeval_subtract(&tvDiff, &tvEnd, &tvBegin);
  printf("\n\n\t\tTime, in seconds, for main() = ");
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



     






