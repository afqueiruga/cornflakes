#include "assemble.h"

#include <stdio.h>

void collect(real_t * loc_in, kernel_t * ke, int hx, int* edge, int l_edge, real_t ** data) {
  int j,i,A;
  int offs = 0;
  for(j=0; j<ke->ninp; j++) {
    switch(ke->inp[j].loc) {
    case LOC_NODE:
      for(A=0;A<l_edge;A++) {
	for(i=0; i<ke->inp[j].len; i++) {
	  loc_in[offs + A*ke->inp[j].len + i] = data[j][ edge[A]*ke->inp[j].len + i] ;
	}
      }
      offs += ke->inp[j].len*l_edge;
      break;
    case LOC_EDGE:
      for(i=0; i<ke->inp[j].len; i++) {
	loc_in[offs + i] = data[j][ hx*ke->inp[j].len + i ];
      }
      offs += ke->inp[j].len;
      break;
    default:
      for(i=0; i<ke->inp[j].len; i++) {
	loc_in[offs + i] = data[j][ i ];
      }
      offs += ke->inp[j].len;
      break;
    }
  }
}

void assemble_vector(real_t * R, kernel_t * ke,hypergraph_t * hg, int * outmap, real_t ** data)
{
  int hx,i,j,A;
  /*
   * R will be accumalated onto!
   */

  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, hg->l_edge);
  int len_loc_out = kernel_outp_len(ke, hg->l_edge);
  real_t loc_in[len_loc_in]; 
  real_t loc_out[len_loc_out];

  //printf("local sizes are %d %d\n",len_loc_in,len_loc_out);

  //printf("%lf %lf\n", data[0][10],data[1][10]);
  /*
   * Loop over the edges in the graph
   */
  int * edge;
  for(hx=0; hx<hg->n_edge; hx++) {
    edge = Hypergraph_Get_Edge(hg,hx);
    //for(i=0;i<hg->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,hg->l_edge, data);

    //for(i=0;i<len_loc_in;i++) printf("%lf ",loc_in[i]); printf("\n");
    /* Call the kernel */
    ke->assem( loc_in, loc_out );
    
    /* Push the data */
    for(i=0;i<len_loc_out;i++) {
      R[ outmap[len_loc_out*hx + i] ] += loc_out[i];
    }
  }
  
}


void assemble_matrix(int * II, int * JJ, real_t * KK,
		     kernel_t * ke, hypergraph_t * hg,
		     int * outmap, real_t ** data)
{
  int i,j,hx;
  /* II,JJ,KK better have size len_loc_out^2 * hg->n_edge */
  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, hg->l_edge);
  int len_loc_out = kernel_outp_len(ke, hg->l_edge);
  real_t loc_in[len_loc_in]; 
  // Don't need this guy, we just plop it into KK
  //real_t loc_out[len_loc_out*len_loc_out]; // FOR NOW: ONLY SQUARE MATRICES

  /*
   * Loop over the edges in the graph
   */
  int * edge;
  
  for(hx=0; hx<hg->n_edge; hx++) {
    edge = Hypergraph_Get_Edge(hg,hx);
    //for(i=0;i<hg->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,hg->l_edge, data);

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

  
}


void assemble_vector_matrix(real_t * R,
			    int * II, int * JJ, real_t * KK,
			    kernel_t * ke, hypergraph_t * hg,
			    int * outmap, real_t ** data)
{
    int i,j,hx;
  /* II,JJ,KK better have size len_loc_out^2 * hg->n_edge */
  /*
   * Allocate local vectors
   */
  int len_loc_in = kernel_inp_len(ke, hg->l_edge);
  int len_loc_out = kernel_outp_len(ke, hg->l_edge);
  real_t loc_in[len_loc_in]; 
  real_t loc_out[len_loc_out + len_loc_out*len_loc_out]; // FOR NOW: ONLY SQUARE MATRICES

  /*
   * Loop over the edges in the graph
   */
  int * edge;
  
  for(hx=0; hx<hg->n_edge; hx++) {
    edge = Hypergraph_Get_Edge(hg,hx);
    //for(i=0;i<hg->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    collect(loc_in, ke, hx,edge,hg->l_edge, data);

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
}
