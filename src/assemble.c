#include "assemble.h"

#include <stdlib.h>
#include <stdio.h>

void collect(real_t * ker_in, kernel_t * ke, hypervertex_t* edge, int l_edge,
	     dofmap_t ** dms, cfdata_t ** data)
{
  hypervertex_t V;
  int i, j, k, fnum, mnum, maxlen;
  dofmap_t * dmap;
  cfdata_t * datum;
  real_t * ker_in_iter = ker_in;
  
  hypervertex_t select[l_edge];
  int nselect, dim;
  
  for(i=0;i<ke->ninp;i++) {
    
    fnum = ke->inp[i].field_number;
    mnum = ke->inp[i].map_num;
    //printf("inp %d on %d\n",i,fnum);
    dmap = dms[mnum];
    //printf("%lx : %d \n",dmap, dmap->U.strided.stride);
    datum = data[fnum];
    maxlen = Dofmap_Max_Len(dmap);
    int dofs[maxlen];
    int ndof;

    k_map_t  kmap = ke->maps[ mnum ];
    kmap(edge,l_edge, select,&nselect, &dim);
    //printf("map: "); for(j=0;j<nselect;j++) printf("%d ",select[j]); printf("\n");
    for(j=0; j<nselect; j++) {
      V = select[j];
      Dofmap_Get(dmap, V, dofs,&ndof);
      CFData_Get_Values(datum, ndof,dofs, ker_in_iter);
      //for(k=0;k<ndof;k++) printf("%d ",dofs[k]);//ker_in_iter[k] = datum[ dofs[k] ];
      ker_in_iter += ndof;
    }
    
  } // End loop over inps

}
void place_targets(target_t * att,
		   kernel_t * ke,
		   real_t * ker_out, int len_ker_out,
		   dofmap_t ** dms,
		   hypervertex_t * edge, int l_edge)
{
  int i,j,k, t,m;
  dofmap_t * dmap;
  int mnum;
  hypervertex_t V;
  real_t * ker_out_iter = ker_out;
  int ndof, maxlen;

  hypervertex_t select[l_edge];
  int nselect, dim;
  
  // Loop over the targets
  for(t=0; t<ke->noutp; t++) {
    // Make the array of all of the DOFs for simplicity
    int nalldofs = kernel_outp_ndof(ke, ke->outp + t, l_edge); // BUG, YES the segfault
    int alldofs[nalldofs]; //TODO: This should be in a routine
    int iter=0;
    for(m=0; m<ke->outp[t].nmap; m++) {
      mnum = ke->outp[t].map_nums[m];
      //printf("t %d mnum %d\n",t,mnum);
      k_map_t kmap = ke->maps[ mnum ];
      kmap(edge,l_edge, select,&nselect, &dim);
      
      dmap = dms[mnum];
      maxlen = Dofmap_Max_Len(dmap);
      int dofs[maxlen];

      for(j=0; j<nselect; j++) { 
	V = select[j];
	Dofmap_Get(dmap, V, dofs,&ndof);
	for(k=0;k<ndof;k++) {
	    alldofs[iter+k] = dofs[k];
	}
	iter+=ndof;
      }
    } // end map loop
    //for(m=0;m<nalldofs;m++) printf("%d ",alldofs[m]); printf("\n");
    // Now assemble into att[t]:

    {
    ker_out_iter = Target_Place(att+t, nalldofs,alldofs, ker_out_iter);
    }
  } // end target loop
}


#include "cfmat_default.h"
#include "cfdata_default.h"
#include <omp.h>
#include "toctic.h"
void assemble(kernel_t * ke, hypergraph_t * hg,
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

