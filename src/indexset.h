#ifndef __INDEXSET_H
#define __INDEXSET_H
typedef int index_t;
typedef struct indexset_str {
  int n;
  int nalloc;
  index_t * table;
} indexset_t;

void IndexSet_New(indexset_t * self, int nalloc);
int IndexSet_Insert(indexset_t * self, index_t i);
void IndexSet_Destroy(indexset_t * self);

#endif
