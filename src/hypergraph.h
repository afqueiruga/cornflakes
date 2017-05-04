#ifndef __H_HYPERGRAPH__
#define __H_HYPERGRAPH__

/*
 * This is a serial library for generic hypergraphs. 
 * The edges of the hypergraph are meant to represent a computation,
 * and the vertices are mappings to DOFs of data to be collected from
 * the global data structures to the local computation kernel.
 * 
 * The is _NO_ geometric information associated with the graphs, they 
 * are a superset of meshes. Vertices of the graphs do not need to be 
 * a nodal-vertex-point-thingy, and hyperedges don't need to an 
 * edge/face/cell/element/w/e
 * 
 * This data structure is meant to stay on a single processor. The 
 * distribution will be handled by a second c module that converts to
 * METIS/SCOTCH style data structures adjancy information, and then 
 * redistributes the union of all the hypergraphs.
 */


/*
 * The data type of vertices.
 * An edge is just a pointer to vertices with an implied length
 * Vertices don't actually have to be anything; they should just be numbers.
 * The numbers can even start to refer to edges if you define 
 * the numbering scheme, e.g.
 *  0  1  2 ... Nv-1  Nv Nv+1 ... Nv+Ne-1
 * V0 V1 V2     VNv   E0  E1        ENe
 */
typedef int hypervertex_t;


/*
 * This structure represents a set of hyperedges of a fixed 
 */
typedef struct hyperedges_t {
  //int n_vertex; //How many vertices are referenced? The graph doesn't care
  int n_edge; //Number of edges
  int n_e_alloc; //How much memory has been allocated
  int l_edge; //How hyper is the edge? Number of verts per edge
  hypervertex_t * edges; //Pointer to the edges
} hyperedges_t;

/* 
 * A hypergraph is a set of hyperedge collections.
 * The hyperedges are sorted by their length for simplicity.
 * Does that make it simpler? It def. makes it cheaper...
 */
typedef struct Hypergraph {
  int n_types; // How many different collections are there?
  int n_alloc; // How many collections have we preallocated?
  hyperedges_t * he;
} hypergraph_t;

/*
 * Hypegraph routines
 */
void Hypergraph_Alloc(hypergraph_t * hg, int n_types);
void Hypergraph_Push_Edge(hypergraph_t * hg, int l_edge, hypervertex_t * verts);
void Hypergraph_Destroy(hypergraph_t * hg);

/*
 * Hyperedge routines
 */
void Hyperedges_Alloc(hyperedges_t * he, int l_edge, int alloc_init);
void Hyperedges_Push_Edge(hyperedges_t * he, int * verts);
hypervertex_t * Hyperedges_Get_Edge(hyperedges_t * he, int i);
void Hyperedges_Destroy(hyperedges_t * he);
#endif
