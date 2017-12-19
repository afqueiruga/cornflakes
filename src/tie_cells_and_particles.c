#include "tie_cells_and_particles.h"

#include "target.h"
#include "stdlib.h"
#include "math.h"
#include "spatialhash.h"
// Old code is packed in here for just the one version...
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
    dmap = dms[mnum];
    datum = data[fnum];
    maxlen = Dofmap_Max_Len(dmap);
    int dofs[maxlen];
    int ndof;

    k_map_t  kmap = ke->maps[ mnum ];
    kmap(edge,l_edge, select,&nselect, &dim);
    for(j=0; j<nselect; j++) {
      V = select[j];
      Dofmap_Get(dmap, V, dofs,&ndof);
      CFData_Get_Values(datum, ndof,dofs, ker_in_iter);
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
  for(t=0; t<ke->noutp; t++) {
    int nalldofs = kernel_outp_ndof(ke, ke->outp + t, l_edge);
    int alldofs[nalldofs];
    int iter=0;
    for(m=0; m<ke->outp[t].nmap; m++) {
      mnum = ke->outp[t].map_nums[m];
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
    {
    ker_out_iter = Target_Place(att+t, nalldofs,alldofs, ker_out_iter);
    }
  } // end target loop
}
void assemble(kernel_t * ke, hypergraph_t * hg,
              dofmap_t ** dofmaps, cfdata_t ** data,
              target_t * att)
{
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
    for(hex=0; hex<he->n_edge; hex++) {
      edge = Hyperedges_Get_Edge(he, hex);
      /* Collect the data */
      collect(ker_in, ke, edge,he->l_edge, dofmaps,data); // TODO: Optimize by moving some oveheard outside of loop
      /* Calculate the kernel */
      for(i=0;i<len_ker_out;i++) ker_out[i] = 0.0;
      ke->eval(he->l_edge, ker_in, ker_out);
      /* Push the data */
      place_targets(att, ke, ker_out,len_ker_out,
                    dofmaps, edge, he->l_edge);
    }
  }
}

/*
 * This routine assigns the particles in x to a cell in mesh
 * If the particles isn't inside of a mesh, it would be in the output.
 * It returns a graph that looks like
 * [ N1 ... don't care... Nn ] P0 P1 P2 ... Pn
 * where it is agnostic to the edge of the graph mesh, and just copies it over.
 * 
 * All it needs is a mesh, (doesn't really need to be a mesh) and 3 kerenls that 
 * gives the circumradius (or an upper bound), gives the centroid (or something 
 * close enough), and tests wether or not a list of points is inside of it.
 * They all need to have the same dofmap call signatures.
 */
#include "cfdata_default.h"
void Tie_Cells_and_Particles(hypergraph_t * hgnew,
			     hypergraph_t * mesh,
			     kernel_t * ke_circum,
			     kernel_t * ke_centroid,
			     kernel_t * ke_inside,
			     dofmap_t ** dofmaps,
			     cfdata_t ** data,
			     int Npart, int dim, real_t * x,
			     hypervertex_t * PV)
{
  int i, j;
  int n_cells = mesh->he[0].n_edge; // TODO: BUG. Hypergraph needs a total_edge_count method

  /* Setup storage for the results */
  target_t att[1];
  cfdata_t hs;
  CFData_Default_New(&hs, n_cells);
  Target_New_From_Ptr(att+0, 1,&hs);
  //setup_targets(ke_circum,att, mesh,n_cells);
  
  /* Calculate the maximum circum radius */
  assemble(ke_circum,mesh, dofmaps,data, att);
  real_t hmax = 0.0;
  for(int i=0; i<n_cells; i++) {
    if(CFData_Default_Data(&hs)[i] > hmax) {
      hmax = CFData_Default_Data(&hs)[i];
    }
  }
  hmax = 2.0*sqrt(hmax);
  //printf("hmax: %lf\n",hmax);
  /* Build a hash using this circumradius and fill it with the particles */
  spatialhash_t sh;
  Build_New_Hash(&sh, Npart,dim,x, hmax); // There's a crash here
  int * found = calloc( sizeof(int),Npart);
  
  /* Now loop against the cells */
  hypervertex_t * edge;
  hyperedges_t * he;
  int hex;
  int listbuf = 20;
  int * list = calloc(sizeof(int),listbuf); 
  for(he=mesh->he; he<mesh->he+mesh->n_types ; he++) {
    int l_edge = he->l_edge;
    int len_ker_centroid = kernel_inps_len(ke_centroid, l_edge); 
    real_t ker_in[ len_ker_centroid ];
    real_t centroid[dim]; // Kernel doesn't have the choice... the output has to look like this, per the contract.
    
    for(hex=0; hex<he->n_edge; hex++) {
      edge = Hyperedges_Get_Edge(he, hex);

      /* Calculate the centroid of the cell */
      collect(ker_in,ke_centroid,edge,l_edge, dofmaps,data);
      for(i=0;i<dim;i++) centroid[i]=0.0;
      ke_centroid->eval(l_edge, ker_in, centroid);
      //printf("center: %lf %lf  ",centroid[0],centroid[1]);
      /* build a list of all the particles _near_ a cell */
      int Nlist=0;
      void action(int FOO, int b) {
		if (Nlist >= listbuf) {
		  listbuf *= 2;
		  list = realloc(list,sizeof(int)*listbuf);
		}
		list[Nlist] = b;
		Nlist++;
      }
      SpatialHash_ScanPt(&sh, centroid, action);
      //printf("nlist: %d\n", Nlist); 
      /* And prune the list by performing the check */
      int len_ker_inside = kernel_inps_len(ke_inside, l_edge+Nlist);
      real_t ker_in_inside[len_ker_inside + dim*Nlist];
      // I just do it myself...
      for(i=0;i<len_ker_centroid;i++) ker_in_inside[i] = ker_in[i];
      for( i=0; i<Nlist;i++) {
		for( j=0; j<dim; j++) {
		  ker_in_inside[len_ker_centroid + dim*i + j] = x[ dim*list[i] + j];
		}
      }
      /* for(i=0;i<len_ker_centroid+dim*Nlist;i++) printf("%lf ",ker_in_inside[i]); printf("\n"); */
      real_t testeval[Nlist];
      for(i=0;i<Nlist;i++) testeval[i]=0.0;
      ke_inside->eval(l_edge+Nlist, ker_in_inside, testeval);
      //for(i=0;i<Nlist;i++) printf("%lf ", testeval[i]);printf("\n"); 
      
      /* Push the result to the graph */
      int Naccept=0;
      hypervertex_t newedge[l_edge+Nlist];
      for(i=0;i<l_edge;i++) newedge[i] = edge[i];
      for(i=0;i<Nlist;i++) {
		if (testeval[i]>0 && found[list[i]]==0) {
		  newedge[l_edge+Naccept] = (PV? PV[list[i]] : list[i]);
		  found[list[i]]=1;
		  Naccept++;
		}
      }

      //printf("Naccept %d\n",Naccept);
      if(Naccept > 0) {
		Hypergraph_Push_Edge(hgnew,  l_edge+Naccept, newedge);
      }
    } // end edge loop
  } // end subgraph loop
  free(list);

  /* Clean up */
  CFData_Destroy(&hs);
  SpatialHash_destroy(&sh);
  Target_Destroy(att+0);
  free(found);
}
