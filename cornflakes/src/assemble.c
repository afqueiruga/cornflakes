#include "assemble.h"

#include <stdlib.h>
#include <stdio.h>

void collect(real_t * ker_in, kernel_t * ke, hypervertex_t* edge, int l_edge,
	     dofmap_t ** dms, real_t ** data)
{
  hypervertex_t V;
  int i, j, k, fnum, mnum, maxlen;
  dofmap_t * dmap;
  real_t * datum;
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
    
    for(j=0; j<nselect; j++) {
      V = edge[j];
	Dofmap_Get(dmap, V, dofs,&ndof);
	for(k=0;k<ndof;k++) ker_in_iter[k] = datum[ dofs[k] ];
	ker_in_iter += ndof;
    }
    
  } // End loop over inps

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
void place_targets(assemble_target_t * att,
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
	V = edge[j];
	Dofmap_Get(dmap, V, dofs,&ndof);
	for(k=0;k<ndof;k++) {
	    alldofs[iter+k] = dofs[k];
	}
	iter+=ndof;
      }
    } // end map loop
    //for(m=0;m<nalldofs;m++) printf("%d ",alldofs[m]); printf("\n");
    // Now assemble into att[t]:
    ker_out_iter = place_target(att+t, alldofs,nalldofs, ker_out_iter);
  } // end target loop
}


void setup_targets(kernel_t * ke, assemble_target_t * att, hypergraph_t * hg, int ndof) {
  int j,i, matsize, oplen;
  for(j=0;j<ke->noutp;j++) {
    switch(ke->outp[j].rank) {
    case 0:
      att[j].rank = 0;
      att[j].V = malloc( sizeof(double) );
      break;
    case 1:
      att[j].rank = 1;
      att[j].V = malloc( ndof*sizeof(double) );
      break;
    case 2:
      matsize = 0;
      for( i=0; i<hg->n_types; i++ ) {
	oplen = kernel_outp_len(ke,ke->outp+j,hg->he[i].l_edge);
	matsize += hg->he[i].n_edge * oplen;
      }
      att[j].rank = 2;
      att[j].V = malloc( matsize * sizeof(double) );
      att[j].II = malloc( matsize * sizeof(int) );
      att[j].JJ = malloc( matsize * sizeof(int) );
      
      break;
    default:
      printf("Cannot handle higher rank\n!");
    }
  }
}
void destroy_targets(kernel_t * ke, assemble_target_t * att) {
  int j;
  for(j=0;j<ke->noutp;j++) {
    free(att[j].V);
    if(att[j].rank>=2) {
      free(att[j].II);
      free(att[j].JJ);
    }
  }
}

void assemble_targets(kernel_t * ke, hypergraph_t * hg,
		      dofmap_t ** dofmaps, real_t ** data,
		      assemble_target_t * att)
{
  int i,j, hex,hx;
  hyperedges_t * he;
  /* Reset the iterators */
  for(i=0;i<ke->noutp;i++) {
    if(att[i].rank>=2) {
      att[i].Viter = att[i].V;
      att[i].IIiter = att[i].II;
      att[i].JJiter = att[i].JJ;
    }
  }
  /* Loop over the graph sets */
  for(he = hg->he; he < hg->he+hg->n_types ; he++) {
    /* Allocate the loca vectors for this size edge */
    int len_ker_in = kernel_inps_len(ke, he->l_edge);
    int len_ker_out = kernel_outps_len(ke, he->l_edge);
    real_t ker_in[ len_ker_in];
    real_t ker_out[len_ker_out];
    //printf("B %d %d\n", len_ker_in, len_ker_out);
    /* Loop over the edges */
    hypervertex_t * edge;
    for(hex=0; hex<he->n_edge; hex++) {
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
      place_targets(att, ke, ker_out,len_ker_out,
		    dofmaps, edge, he->l_edge);
    }
  }

  
}
