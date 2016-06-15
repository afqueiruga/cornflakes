#include "toctic.h"
#include <sys/time.h>
#include <stdlib.h>
void toctic(const char * str) {
  static double lasttic = -1.0;
  struct timeval t;
  gettimeofday(&t, NULL);
  double toc =  1.0*t.tv_sec + 1e-6*t.tv_usec;
  if(lasttic > 0.0) {
    printf("Time for %s : %lf s \n",str,toc-lasttic);
  }
  lasttic = toc;
}
void toctic_ts(const char * str, double * lasttic) {
  struct timeval t;
  gettimeofday(&t, NULL);
  double toc =  1.0*t.tv_sec + 1e-6*t.tv_usec;
  if(*lasttic > 0.0) {
    printf("Time for %s : %lf s \n",str,toc-*lasttic);
  }
  *lasttic = toc;
}
