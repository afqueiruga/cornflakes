%module mylibrary
%{
#define SWIG_FILE_WITH_INIT
#include "hypergraph.h"
#include "SpatialHash.h"
#include "Graphers.h"
#include "kernel.h"
#include "assemble.h"
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (int DIM1, int* IN_ARRAY1) {(int nvert, int * verts)};
%apply (int DIM1, int DIM2, real_t* IN_ARRAY2) {(int Npart, int dim, real_t * x)};
%include "carrays.i"
%array_class(int,intArray)
%array_class(dof_t,dofArray)
%include "cpointer.i"
%pointer_class(int, intp)

%include "hypergraph.h"
%include "SpatialHash.h"
%include "Graphers.h"
%include "kernel.h"
%include "assemble.h"


%exception Hypergraph_Push_Edge_np {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}
%exception Hypergraph_Get_Edge_np {
    $action
    if (PyErr_Occurred()) SWIG_fail;
}

%inline %{
  void Hypergraph_Push_Edge_np(hypergraph_t * hg, int nvert, int * verts) {
    if (nvert < hg->l_edge) {
        PyErr_Format(PyExc_ValueError,
                     "Array of insufficent length (verts %d < edge length %d)",
                     nvert, hg->l_edge);
        return;
    }
    Hypergraph_Push_Edge(hg,verts);
  }
  void Hypergraph_Get_Edge_np(hypergraph_t * hg, int i,
			      int * DIM1,int** ARGOUTVIEW_ARRAY1) {
    if(i<0 || i>=hg->n_edge) {
        PyErr_Format(PyExc_ValueError,
                     "Bad index: require 0<%d<%d",
                     i, hg->n_edge);
        return;
    }
    *DIM1 = hg->l_edge;
    *ARGOUTVIEW_ARRAY1 = Hypergraph_Get_Edge(hg,i);
  }

  void Hypergraph_Get_View_np(hypergraph_t * hg,
			       int * DIM1, int *DIM2, int** ARGOUTVIEW_ARRAY2) {
    *DIM1 = hg->n_edge;
    *DIM2 = hg->l_edge;
    *ARGOUTVIEW_ARRAY2 = hg->edges;
  }
  //void Hypergraph_Get_Edge

  void assemble_vector_np(int DIM1, real_t * INPLACE_ARRAY1,
			  kernel_t * ke, hypergraph_t * hg,
			  int * INPLACE_ARRAY_FLAT, int DIM_FLAT,
			  PyObject * datalist) {
    if(PyList_Check(datalist)) {
      int ndata = PyList_Size(datalist);
      
      real_t *data_ptrs[ndata];
      PyArrayObject * nda[ndata];
      int isnewobj[ndata];
      int i;
      for(i=0;i<ndata;i++) {
	PyObject * obj = NULL;
	nda[i] = NULL;
	obj = PyList_GetItem(datalist, i);
	isnewobj[i];
	nda[i] = obj_to_array_contiguous_allow_conversion(obj,NPY_DOUBLE,isnewobj+i);
	printf("%d ", array_size(nda[i],0));
	data_ptrs[i] = array_data(nda[i]);
      }
      printf("%d\n",ndata);

      // Call the assembler
      assemble_vector(INPLACE_ARRAY1, ke,hg, INPLACE_ARRAY_FLAT, data_ptrs);
      
      // Tell python we don't want those arrays anymore, if we actually made any
      for(i=0;i<ndata;i++) {
	if(isnewobj[i]) Py_DECREF( nda[i] );
      }
    }
  }
%}
