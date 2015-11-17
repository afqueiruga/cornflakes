#include "assemble.h"

#include <stdio.h>

void collect(real_t * ker_in, kernel_t * ke, hypervertex_t* edge, int l_edge,
	     dofmap_t ** dms, real_t ** data)
{
  hypervertex_t V;
  int i, j, k, fnum, maxlen;
  dofmap_t * dmap;
  real_t * datum;
  real_t * ker_in_iter = ker_in;
  // TODO: Fringe case of global and variable length data.
  for(i=0;i<ke->ninp;i++) {
    fnum = ke->inp[i].field_number;
    dmap = dms[fnum];
    datum = data[fnum];
    maxlen = Dofmap_Max_Len(dmap);
    int dofs[maxlen];
    int ndof;
    for(j=0;j<l_edge;j++) { // BUG
      V = edge[j];
      Dofmap_Get(dmap, V, dofs,&ndof);
      for(k=0;k<ndof;k++) ker_in_iter[k] = datum[ dofs[k] ];
      ker_in_iter += ndof;
    }
  }

}

real_t * place_target(assemble_target_t * att,
		      int * dofs, int n,
		      real_t * ker_out) {
  int i,j;
  switch(att->rank) {
  case 2:
    // Fill this block
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
	att->IIiter[n*i + j ] = dofs[i];
	att->JJiter[n*i + j ] = dofs[j];
	att->Viter [n*i + j ] = ker_out[ n*i + j];
      }
    }
    // Advance our iterators.
    att->IIiter += n*n;
    att->JJiter += n*n;
    att->Viter += n*n;
    return ker_out + n*n;
    break;
  case 1:
    for(i=0;i<n;i++) {
      att->V[dofs[i]] += ker_out[i];
    }
    return ker_out + n;
    break;
  default:
    for(i=0;i<n;i++) {
      att->V[0] += ker_out[i];
    }
    return ker_out + 1;
    break;
  }

}
void place_targets(assemble_target_t ** atts,
		   kernel_t * ke,
		   real_t * ker_out, int len_ker_out,
		   dofmap_t ** dms,
		   hypervertex_t * edge, int l_edge)
{
  int i,j,k, t,d;
  dofmap_t * dmap;
  int fnum;
  hypervertex_t V;
  real_t * ker_out_iter;
  int ndof, maxlen;
  // Loop over the targets
  for(t=0; t<ke->noutp; t++) {
    // Make the array of all of the DOFs for simplicity
    int nalldofs = kernel_outp_len(ke->outps + t,l_edge);
    int alldofs[nalldofs]; //TODO: This should be in a routine
    int iter=0;
    for(d=0; d<ke->outps[t].ndof; d++) {
      fnum = ke->outps[t].dofs[d].field_number;
      dmap = dms[fnum];
      maxlen = Dofmap_Max_Len(dmap);
      int dofs[maxlen];
      // Loop over the vertices
      for(j=0;j<l_edge;j++) { // BUG
	V = edge[j];
	Dofmap_Get(dmap, V, dofs,&ndof);
	for(k=0;k<ndof;k++) {
	  alldofs[iter+k] = dofs[k];
	}
	iter+=ndof;
      }
    }
    // Now assemble in:
    ker_out_iter = place_target(atts[t], alldofs,nalldofs, ker_out_iter);
  } // end target loop
}

void assemble_targets(kernel_t * ke, hypergraph_t * hg,
		      dofmap_t ** dofmaps, real_t ** data,
		      assemble_target_t * att)
{
  int i,j, hex,hx;
  hyperedges_t * he;
  
  /* Loop over the graph sets */
  for(he = hg->he; he < hg->he+hg->n_types ; he++) {
    /* Allocate the loca vectors for this size edge */
    int len_ker_in = kernel_inp_len(ke, he->l_edge);
    int len_ker_out = kernel_outps_len(ke, he->l_edge);
    real_t ker_in[ len_ker_in];
    real_t ker_out[len_ker_out];

    /* Loop over the edges */
    hypervertex_t * edge;
    for(hex=0; hex<he->n_edge; hex++) {
      edge = Hyperedges_Get_Edge(he, hex);

      /* Collect the data */
      collect(ker_in, ke, edge,he->l_edge, dofmaps,data); // TODO: Optimize by moving some overheard outside of loop

      /* Calculate the kernel */
      ke->eval(he->l_edge, ker_in, ker_out);

      /* Push the data */
      
    }
  }

  
}




void assemble_targets_dep(int ntarget, assemble_target_t * att,
			  kernel_t * ke, hyperedges_t * hg,
			  int * outmap, real_t ** data)
{
  #if 0
  int i,j,hx;
  /* II,JJ,KK better have size len_loc_out^2 * he->n_edge */
  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, he->l_edge);
  int len_loc_out = kernel_outp_len(ke, he->l_edge);
  real_t loc_in[len_loc_in];
  int out_alloc = 0; // THE KERNEL SHOULD ENCODE THIS INFORMATION....
  for(i=0;i<ntarget;i++) {
    switch(att[i].rank) {
    case 0: out_alloc += 1; break;
    case 1: out_alloc += len_loc_out; break;
    case 2: out_alloc += len_loc_out*len_loc_out; break; // FOR NOW: ONLY SQUARE MATRICES
    }
  }
  real_t loc_out[out_alloc];

  /*
   * Loop over the edges in the graph
   */
  hypervertex_t * edge;
  for(hx=0; hx<he->n_edge; hx++) {
    edge = Hyperedges_Get_Edge(he,hx);
    //for(i=0;i<he->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,he->l_edge, data);

    //for(i=0;i<len_loc_in;i++) printf("%lf ",loc_in[i]); printf("\n");
    /* Call the kernel */
    ke->assem( loc_in, loc_out );
    //for(i=0;i<len_loc_out;i++) printf("%lf ",loc_out[i]); printf("\n");
    real_t * iter_out = loc_out;
    for(i=0;i<ntarget;i++) {
      iter_out = push_target(att+i,len_loc_out, hx,outmap+len_loc_out*hx, iter_out);
    }
  }
  #endif
}




/************************************
 ***         DEPRECATED:          ***
 ************************************/


void assemble_vector(real_t * R, kernel_t * ke,hyperedges_t * he, int * outmap, real_t ** data)
{
#if 0
  int hx,i,j,A;
  /*
   * R will be accumalated onto!
   */

  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, he->l_edge);
  int len_loc_out = kernel_outp_len(ke, he->l_edge);
  real_t loc_in[len_loc_in]; 
  real_t loc_out[len_loc_out];

  //printf("local sizes are %d %d\n",len_loc_in,len_loc_out);

  //printf("%lf %lf\n", data[0][10],data[1][10]);
  /*
   * Loop over the edges in the graph
   */
  int * edge;
  for(hx=0; hx<he->n_edge; hx++) {
    edge = Hyperedges_Get_Edge(he,hx);
    //for(i=0;i<he->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,he->l_edge, data);

    //for(i=0;i<len_loc_in;i++) printf("%lf ",loc_in[i]); printf("\n");
    /* Call the kernel */
    ke->assem( loc_in, loc_out );
    
    /* Push the data */
    for(i=0;i<len_loc_out;i++) {
      R[ outmap[len_loc_out*hx + i] ] += loc_out[i];
    }
  }
  #endif
}


void assemble_matrix(int * II, int * JJ, real_t * KK,
		     kernel_t * ke, hyperedges_t * he,
		     int * outmap, real_t ** data)
{
  #if 0
  int i,j,hx;
  /* II,JJ,KK better have size len_loc_out^2 * he->n_edge */
  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, he->l_edge);
  int len_loc_out = kernel_outp_len(ke, he->l_edge);
  real_t loc_in[len_loc_in]; 
  // Don't need this guy, we just plop it into KK
  //real_t loc_out[len_loc_out*len_loc_out]; // FOR NOW: ONLY SQUARE MATRICES

  /*
   * Loop over the edges in the graph
   */
  hypervertex_t * edge;
  
  for(hx=0; hx<he->n_edge; hx++) {
    edge = Hyperedges_Get_Edge(he,hx);
    //for(i=0;i<he->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,he->l_edge, data);

    //for(i=0;i<len_loc_in;i++) printf("%lf ",loc_in[i]); printf("\n");
    /* Call the kernel */
    ke->assem( loc_in, KK+len_loc_out*len_loc_out*hx );
    
    /* Push the II,JJ indices */
    for(i=0;i<len_loc_out;i++) {
      for(j=0;j<len_loc_out;j++) {
	II[ len_loc_out*len_loc_out*hx + len_loc_out*i + j] = outmap[len_loc_out*hx + i ];
	JJ[ len_loc_out*len_loc_out*hx + len_loc_out*i + j] = outmap[len_loc_out*hx + j ];
      }
    }
  }

 #endif 
}


void assemble_vector_matrix(real_t * R,
			    int * II, int * JJ, real_t * KK,
			    kernel_t * ke, hyperedges_t * he,
			    int * outmap, real_t ** data)
{
  #if 0
    int i,j,hx;
  /* II,JJ,KK better have size len_loc_out^2 * he->n_edge */
  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, he->l_edge);
  int len_loc_out = kernel_outp_len(ke, he->l_edge);
  real_t loc_in[len_loc_in]; 
  real_t loc_out[len_loc_out + len_loc_out*len_loc_out]; // FOR NOW: ONLY SQUARE MATRICES

  /*
   * Loop over the edges in the graph
   */
  hypervertex_t * edge;
  
  for(hx=0; hx<he->n_edge; hx++) {
    edge = Hyperedges_Get_Edge(he,hx);
    //for(i=0;i<he->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,he->l_edge, data);

    //for(i=0;i<len_loc_in;i++) printf("%lf ",loc_in[i]); printf("\n");
    /* Call the kernel */
    ke->assem( loc_in, loc_out );
    
    /* Push the II,JJ indices */
    for(i=0;i<len_loc_out;i++) {
      for(j=0;j<len_loc_out;j++) {
	II[ len_loc_out*len_loc_out*hx + len_loc_out*i + j] = outmap[len_loc_out*hx + i ];
	JJ[ len_loc_out*len_loc_out*hx + len_loc_out*i + j] = outmap[len_loc_out*hx + j ];
	KK[len_loc_out*len_loc_out*hx + len_loc_out*i + j] = loc_out[len_loc_out + len_loc_out*i + j];
      }
      R[outmap[len_loc_out*hx + i]] += loc_out[i];
    }
    
  }
#endif
}

