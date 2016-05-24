#ifdef USE_LIS
#ifndef __LIS_TARGET_H
#define __LIS_TARGET_H

#include "target.h"
#include <lis.h>

void Target_LIS_From_Obj(target_t * self, int rank, void * obj);
typedef struct Target_LIS_data_t {
  int own;
  LIS_VECTOR R;
  LIS_MATRIX K;
} Target_LIS_data_t;
#define Target_LIS_Data(x) ((Target_LIS_data_t*)((x)->data))

#endif
#endif
