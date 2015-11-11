%module cornflakes_library
%{
#define SWIG_FILE_WITH_INIT
#include "hypergraph.h"
#include "spatialhash.h"
#include "graphers.h"
#include "kernel.h"
#include "assemble.h"
%}

%include "numpy.i"
%init %{
import_array();
%}

%apply (int DIM1, int* IN_ARRAY1) {(int nvert, hypervertex_t * verts),
                                   (int l_edge, hypervertex_t * verts)};
%apply (int DIM1, int DIM2, real_t* IN_ARRAY2) {(int Npart, int dim, real_t * x)};
%apply (int DIM1, int* INPLACE_ARRAY1) {(int dim_II, int * array_II),
                                        (int dim_JJ, int * array_JJ)};
%apply (int DIM1, real_t* INPLACE_ARRAY1) {(int dim_KK, real_t * array_KK)};
%include "carrays.i"
%array_class(int,intArray)
%array_class(dof_t,dofArray)
%include "cpointer.i"
%pointer_class(int, intp)

%include "hypergraph.h"
%include "spatialhash.h"
%include "graphers.h"
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
  /*
   * Extra wrappers for Hypergraph and Hyperedges
   */
  void Hyperedges_Push_Edge_np(hyperedges_t * he, int nvert, int * verts) {
    if (nvert < he->l_edge) {
        PyErr_Format(PyExc_ValueError,
                     "Array of insufficent length (verts %d < edge length %d)",
                     nvert, he->l_edge);
        return;
    }
    Hyperedges_Push_Edge(he,verts);
  }
  void Hyperedges_Get_Edge_np(hyperedges_t * he, int i,
			      int * DIM1,int** ARGOUTVIEW_ARRAY1) {
    if(i<0 || i>=he->n_edge) {
        PyErr_Format(PyExc_ValueError,
                     "Bad index: require 0<%d<%d",
                     i, he->n_edge);
        return;
    }
    *DIM1 = he->l_edge;
    *ARGOUTVIEW_ARRAY1 = Hyperedges_Get_Edge(he,i);
  }

  void Hyperedges_Get_View_np(hyperedges_t * he,
			       int * DIM1, int *DIM2, int** ARGOUTVIEW_ARRAY2) {
    *DIM1 = he->n_edge;
    *DIM2 = he->l_edge;
    *ARGOUTVIEW_ARRAY2 = he->edges;
  }

  void Hypergraph_Get_Edge_View_np(hypergraph_t * hg, int i,
				int * DIM1, int *DIM2, int** ARGOUTVIEW_ARRAY2) {
    Hyperedges_Get_View_np(hg->he+i, DIM1,DIM2,ARGOUTVIEW_ARRAY2);
  }

  /*
   * Assembly wrappers
   */
  void assemble_targets_np(PyObject * targetlist,
			   kernel_t * ke, hyperedges_t * he,
			   int * INPLACE_ARRAY_FLAT, int DIM_FLAT,
			   PyObject * datalist)
  {
    int i=0, isnewobj=0;
    PyObject *obj, *subobj;
    PyArrayObject * arrobj;
    /* Step 0: Check data types and build a list of possible object allocations */
    if(!PyList_Check(targetlist)) return;
    if(!PyList_Check(datalist)) return;
    int ntarget = PyList_Size(targetlist);
    int ndata = PyList_Size(datalist);
    int nnewobj = 0;
    PyArrayObject * newobjs[3*ntarget + ndata];
    
    /* Step 1: Build the target list */
    assemble_target_t att[ntarget];
    for(i=0;i<ntarget;i++) {
      //printf("Target %d ",i);
      obj = PyList_GetItem(targetlist, i);
      if(PyTuple_Check(obj)) {
	//printf("Is a tuple\n");
	att[i].rank = 2;
	// 0: Get KK
	subobj = PyTuple_GetItem(obj,0);
	isnewobj = 0;
	arrobj = obj_to_array_contiguous_allow_conversion(subobj,NPY_DOUBLE,&isnewobj);
	att[i].V = array_data(arrobj);
	if(isnewobj) { newobjs[nnewobj] = arrobj; nnewobj++; }
	// 1: Get II
	subobj = PyTuple_GetItem(obj,1);
	isnewobj = 0;
	arrobj = obj_to_array_contiguous_allow_conversion(subobj,NPY_INT,&isnewobj);
	att[i].II = array_data(arrobj);
	if(isnewobj) { newobjs[nnewobj] = arrobj; nnewobj++; }
	// 2: Get JJ
	subobj = PyTuple_GetItem(obj,2);
	isnewobj = 0;
	arrobj = obj_to_array_contiguous_allow_conversion(subobj,NPY_INT,&isnewobj);
	att[i].JJ = array_data(arrobj);
	if(isnewobj) { newobjs[nnewobj] = arrobj; nnewobj++; }
      }
      else { //  obj is (BETTER BE) a ndarray
	//printf("Is an ndarray\n");
	isnewobj = 0;
	arrobj = obj_to_array_contiguous_allow_conversion(obj,NPY_DOUBLE,&isnewobj);
	att[i].V = array_data(arrobj);
	if(isnewobj) { newobjs[nnewobj] = arrobj; nnewobj++; }
	if(array_size(arrobj,0)>1) {
	  att[i].rank=1;
	} else {
	  att[i].rank=0;
	}
	att[i].II = NULL;
	att[i].JJ = NULL;
      }
    }

    //printf("Here\n");
    /* Step 2: Collect the data ptrs */
    real_t * data_ptrs[ndata];
    for(i=0;i<ndata;i++) {
      obj = PyList_GetItem(datalist,i );
      isnewobj = 0;
      arrobj = obj_to_array_contiguous_allow_conversion(obj,NPY_DOUBLE,&isnewobj);
      data_ptrs[i] = array_data(arrobj);
      if(isnewobj) {
	newobjs[nnewobj] = arrobj;
	nnewobj++;
      }
    }
    //printf("Assembling\n");
    /* Step 3: Assemble */
    assemble_targets(ntarget, att,
		     ke, he,
		     INPLACE_ARRAY_FLAT,
		     data_ptrs);

    //printf("Recreasing references of %d objects\n",nnewobj);
    /* Step 4: Decrease reference counts */
    for(i=0;i<nnewobj;i++) {
      Py_DECREF(newobjs[i]);
    }
    
    /* Step 5: Profit! */
  }

  void assemble_vector_np(int DIM1, real_t * INPLACE_ARRAY1,
			  kernel_t * ke, hyperedges_t * he,
			  int * INPLACE_ARRAY_FLAT, int DIM_FLAT,
			  PyObject * datalist)
  {
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
	isnewobj[i] = 0;
	nda[i] = obj_to_array_contiguous_allow_conversion(obj,NPY_DOUBLE,isnewobj+i);
	data_ptrs[i] = array_data(nda[i]);
      }
      

      // Call the assembler
      assemble_vector(INPLACE_ARRAY1, ke,he, INPLACE_ARRAY_FLAT, data_ptrs);
      
      // Tell python we don't want those arrays anymore, if we actually made any
      for(i=0;i<ndata;i++) {
	if(isnewobj[i]) Py_DECREF( nda[i] );
      }
    }
  }

  void assemble_matrix_np(int dim_II, int * array_II,
			  int dim_JJ, int * array_JJ,
			  int dim_KK, real_t * array_KK,
			  kernel_t * ke, hyperedges_t * he,
			  int * INPLACE_ARRAY_FLAT, int DIM_FLAT,
			  PyObject * datalist)
  {
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
	isnewobj[i] = 0;
	nda[i] = obj_to_array_contiguous_allow_conversion(obj,NPY_DOUBLE,isnewobj+i);
	data_ptrs[i] = array_data(nda[i]);
      }
      

      // Call the assembler
      assemble_matrix(array_II,array_JJ,array_KK, ke,he, INPLACE_ARRAY_FLAT, data_ptrs);
      
      // Tell python we don't want those arrays anymore, if we actually made any
      for(i=0;i<ndata;i++) {
	if(isnewobj[i]) Py_DECREF( nda[i] );
      }
    }
  }

  void assemble_vector_matrix_np(int DIM1, real_t * INPLACE_ARRAY1,
				 int dim_II, int * array_II,
				 int dim_JJ, int * array_JJ,
				 int dim_KK, real_t * array_KK,
				 kernel_t * ke, hyperedges_t * he,
				 int * INPLACE_ARRAY_FLAT, int DIM_FLAT,
				 PyObject * datalist)
  {
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
	isnewobj[i] = 0;
	nda[i] = obj_to_array_contiguous_allow_conversion(obj,NPY_DOUBLE,isnewobj+i);
	data_ptrs[i] = array_data(nda[i]);
      }
      

      // Call the assembler
      assemble_vector_matrix(INPLACE_ARRAY1,array_II,array_JJ,array_KK, ke,he, INPLACE_ARRAY_FLAT, data_ptrs);
      
      // Tell python we don't want those arrays anymore, if we actually made any
      for(i=0;i<ndata;i++) {
	if(isnewobj[i]) Py_DECREF( nda[i] );
      }
    }
  }
%}