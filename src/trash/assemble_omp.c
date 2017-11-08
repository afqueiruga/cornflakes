#include "assemble_omp.h"

#include "cfmat_default.h"
#include "cfdata_default.h"

#include <omp.h>
#include <stdlib.h>
#include "toctic.h"


#ifdef WIP

typedef struct obj_table_entry_str {
  void * einzel;
  void * doppel;
} obj_table_entry_t;
typedef struct obj_table_str {
  int N;
  obj_table_entry_t tbl[30]; //Meh, better not be this many!
} obj_table_t;

obj_table_entry_t * In_Table(obj_table_t * self, void * addr) {
  int i;
  for(i=0;i<self->N;i++) {
    if(self->table[i].einzel == addr)
      return self->table+i;
  }
  return NULL;
}
void Add_To_Table(obj_table_t * self, void * einzel, void * doppel) {
  self->tbl[self->N].einzel = einzel;
  self->tbl[self->N].doppel = doppel;
  self->N++;
}


// This is the exact opposite of this type discipline!
extern const _CFMAT_VTABLE_t CFMat_CSR_vtable;
extern const _CFMAT_VTABLE_t CFMat_BC_vtable;
extern const _CFDATA_VTABLE_t cfdata_default_vtable;
extern const _CFDATA_VTABLE_t cfdata_bc_vtable;

cfdata_t * Data_Descend(cfdata_t * node, obj_table_t * tbl);
cfmat_t * Mat_Descend(cfdata_t * node, obj_table_t * tbl);

cfdata_t * Data_Descend(cfdata_t * node, obj_table_t * tbl) {
  // Check if I'm in the table
  if(intable) {
    return table entry;
  }
  // Otherwise, make a copy of me based on my type.
  // Descend as needed.
  switch(node->vtable) {
  case cfdata_bc_vtable:
    
    break;
  case cfdata_default_vtable:
  default:
    break;
  }
}
void create_target_doppel(target_t * doppel, target_t * einzel, int ntarg) {
  /* Make a copy of the target data structure with doppelgangers.
     The problem is the damn BC chains...*/
  int i;
  for(i=0;i<ntarg;i++) {
    
  }
}
void merge_target_doppel(target_t * einzel, target_t * doppel) {
  /* Loop through and add them all back in */
}
#endif

void assemble_omp(kernel_t * ke, hypergraph_t * hg,
	      dofmap_t ** dofmaps, cfdata_t ** data,
	      target_t * att)
{
  
  target_t * glob_att;
  #pragma omp parallel shared(glob_att)
  {
    int tid = omp_get_thread_num();
    int tmax = omp_get_num_threads();
    double lasttic = -1.0;
    if(tid==0) {
      glob_att = malloc(sizeof(target_t)*tmax*ke->noutp);
    }
    #pragma omp barrier
    printf("Hi, I'm %d \n",tid);
    /* Make tempory storage owned by just this thread. Big storage cost! */
    target_t * my_att = glob_att + ke->noutp * tid;
    for(int i=0; i<ke->noutp; i++ ) {
      if(att[i].rank==2) {
	cfmat_t * myt = malloc(sizeof(cfmat_t));
	CFMat_Default_New(myt, i,ke,hg, att[i].K->N);
	Target_New_From_Ptr(my_att+i, 2, myt);
      } else {
	cfdata_t * myt = malloc(sizeof(cfdata_t));
	CFData_Default_New(myt, att[i].R->N);
	Target_New_From_Ptr(my_att+i, 1, myt);
      }
    }
    toctic_ts("Starting\n",&lasttic);
    int i,j, hex,hx;
    hyperedges_t * he;
    /* Loop over the graph sets */
    for(he = hg->he; he < hg->he+hg->n_types ; he++) {
      /* Allocate the local vectors for this size edge */
      int len_ker_in = kernel_inps_len(ke, he->l_edge);
      int len_ker_out = kernel_outps_len(ke, he->l_edge);
      real_t ker_in[ len_ker_in];
      real_t ker_out[len_ker_out];

      /* Loop over the edges */
      hypervertex_t * edge;
      
      toctic_ts("he",&lasttic);
      #pragma omp for private(hex)
      for(hex=0; hex<he->n_edge; hex++) {
	
	edge = Hyperedges_Get_Edge(he, hex);
	
	/* Collect the data */
	collect(ker_in, ke, edge,he->l_edge, dofmaps,data); // TODO: Optimize by moving some overheard outside of loop
	
	/* Calculate the kernel */
	for(i=0;i<len_ker_out;i++) ker_out[i] = 0.0;
	ke->eval(he->l_edge, ker_in, ker_out);
	
	/* Push the data */
	place_targets(my_att, ke, ker_out,len_ker_out,
	    		  dofmaps, edge, he->l_edge);
      } // end loop hex
    } // end loop he
    toctic_ts("looping",&lasttic);
    #pragma omp barrier
    /* Merge into arguements in serial att */
    if(tid==0) {
      for(int t=0; t<tmax; t++) {
	target_t * his_att = glob_att + ke->noutp * t;

	for(i=0; i<ke->noutp; i++) {
	  if(his_att[i].rank==2) {
	    CFMat_Default_data_t * trips = CFMat_Default_Data(his_att[i].K);
	    int n_trip = trips->IIiter - trips->II;
	    for(j=0;j<n_trip;j++) {
	      CFMat_Set_Value(att[i].K, trips->II[j],trips->JJ[j], trips->V[j]);
	    }
	    CFMat_Destroy(his_att[i].K);
	    free(his_att[i].K);
	    
	  } else {
	    real_t * v = CFData_Default_Data(his_att[i].R);
	    for(j=0; j<att[i].R->N; j++) {
	      CFData_Place(att[i].R, 1,&j,v+j);
	    }
	    CFData_Destroy(his_att[i].R);
	    free(his_att[i].R);
	  }
	  //Target_Destroy(his_att+i);
	}
      }
      free(glob_att);
    }
    toctic_ts("reassembly",&lasttic);
  }// end parallel block
  
}

