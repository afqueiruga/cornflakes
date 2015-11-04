#include "assemble.h"

#include <stdio.h>

void assemble_vector(real_t * R, kernel_t * ke,hypergraph_t * hg, int * outmap, real_t ** data)
{
  int hx,i,j,A;
  /*
   * R will be accumalated onto!
   */

  /*
   * Allocate local vectors
   */
  int len_loc_in = 0;
  for(i=0; i<ke->ninp; i++) {
    if(ke->inp[i].loc == LOC_NODE) {
      len_loc_in += ke->inp[i].len*hg->l_edge;
    } else {
      len_loc_in += ke->inp[i].len;
    }
  }
  int len_loc_out = 0;
  for(i=0; i<ke->noutp; i++) {
    if(ke->outp[i].loc == LOC_NODE) {
      len_loc_out += ke->outp[i].len*hg->l_edge;
    } else {
      len_loc_out += ke->outp[i].len;
    }
  }
  real_t loc_in[len_loc_in]; 
  real_t loc_out[len_loc_out];

  printf("local sizes are %d %d\n",len_loc_in,len_loc_out);

  printf("%lf %lf\n", data[0][10],data[1][10]);
  /*
   * Loop over the edges in the graph
   */
  int * edge;
  for(hx=0; hx<hg->n_edge; hx++) {
    edge = Hypergraph_Get_Edge(hg,hx);
    for(i=0;i<hg->l_edge;i++) printf("%d ",edge[i]); printf("\n");

    /* Pull data */
    int offs = 0;
    for(j=0; j<ke->ninp; j++) {
      switch(ke->inp[j].loc) {
      case LOC_NODE:
	for(A=0;A<hg->l_edge;A++) {
	  for(i=0; i<ke->inp[j].len; i++) {
	    loc_in[offs + A*ke->inp[j].len + i] = data[j][ edge[A]*ke->inp[j].len + i] ;
	  }
	}
	offs += ke->inp[j].len*hg->l_edge;
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

    for(i=0;i<len_loc_in;i++) printf("%lf ",loc_in[i]); printf("\n");
    /* Call the kernel */

    /* Push the data */
    
  }
  
}
