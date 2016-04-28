#ifndef __TARGET_H
#define __TARGET_H

#include "kernel.h"

/* A polymorphic class for targets with different backends */
typedef struct _TARGET_VTABLE_t _TARGET_VTABLE_t;
typedef struct target_t {
  int rank;
  void * data;
  const _TARGET_VTABLE_t * vtable;
} target_t;

/*
 * The constructor for the built in type, meant to be used with 
 * Scipy in python-cornflakes
 */
void Target_New_Default(target_t * self, int rank);

/* The member methods */
void Target_Push(target_t * self, int n, int * dofs, real_t * vals);
void Target_Destroy(target_t * self);
void Target_Wipe(target_t * self);

#endif
