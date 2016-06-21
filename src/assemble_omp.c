#include "assemble_omp.h"

#include "cfmat_default.h"
#include "cfdata_default.h"

#include <omp.h>
#include <stdlib.h>
#include "toctic.h"
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
      //printf("B %d %d\n", len_ker_in, len_ker_out);
      /* Loop over the edges */
      hypervertex_t * edge;
      
      toctic_ts("he",&lasttic);
      #pragma omp for private(hex)
      for(hex=0; hex<he->n_edge; hex++) {
	//printf("%d %d\n",tid,hex);
	edge = Hyperedges_Get_Edge(he, hex);
	//printf("e %d\n",hex);
	/* Collect the data */
	collect(ker_in, ke, edge,he->l_edge, dofmaps,data); // TODO: Optimize by moving some overheard outside of loop
	//printf("did\n");
	/* Calculate the kernel */
	//printf("in:"); for(i=0;i<len_ker_in;i++) printf("%lf ",ker_in[i]); printf("\n");
	for(i=0;i<len_ker_out;i++) ker_out[i] = 0.0;
	ke->eval(he->l_edge, ker_in, ker_out);
	//printf("out:"); for(i=0;i<len_ker_out;i++) printf("%lf ",ker_out[i]); printf("\n");
	//printf("eval\n");
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

