%module cornflakes_library
%{
#define SWIG_FILE_WITH_INIT
#include "hypergraph.h"
#include "spatialhash.h"
#include "graphers.h"
#include "kernel.h"
#include "assemble.h"
#include "dofmap.h"
#include "util.h"
  
  //#include "kernels/sample_peri.h"
  //#include "kernels/sample_state.h"
  //#include "kernels/darcy_state.h"
  //#include "kernels/darcy_afq_state.h"
  //#include "kernels/darcy_support_afq_state.h"
  //#include "kernels/darcy_afq_CG.h"
%}

%include "numpy.i"
%init %{
  import_array();
  //import_managed();
%}

%apply (int DIM1, int* IN_ARRAY1) {(int nvert, hypervertex_t * verts),
                                   (int l_edge, hypervertex_t * verts)};
%apply (int DIM1, int DIM2, real_t* IN_ARRAY2) {(int Npart, int dim, real_t * x),
     (int npart, int dim, real_t * x),
     (int nx1, int dx1, real_t * x1),
     (int nx2, int dx2, real_t * x2),
     (int nu1, int du1, real_t * u1)
     };
%apply (int DIM1, int  DIM2, real_t * INPLACE_ARRAY2) {
  (int nu2, int du2, real_t * u2)
    };
%apply (int DIM1, int* INPLACE_ARRAY1) {(int dim_II, int * array_II),
                                        (int dim_JJ, int * array_JJ)};
%apply (int DIM1, real_t* INPLACE_ARRAY1) {(int dim_KK, real_t * array_KK)};
%apply (int* ARGOUT_ARRAY1, int * DIM1) {(int * dofs, int *ndofs)};
%apply (int DIM1, int DIM2, int * IN_ARRAY2) { (int Nentry, int stride, int * table) };

%include "carrays.i"
%array_class(int,intArray)
%array_class(k_map_t,k_mapArray)
%array_class(inp_t,inpArray)
%array_class(outp_t,outpArray)
%array_class(hyperedges_t,hyperedgesArray)

%include "cpointer.i"
%pointer_class(int, intp)

%include "hypergraph.h"
%include "spatialhash.h"
%include "graphers.h"
%include "kernel.h"
%include "assemble.h"
%include "dofmap.h"
%include "util.h"

 //%include "kernels/sample_peri.h"
 //%include "kernels/sample_state.h"
 //%include "kernels/darcy_state.h"
 //%include "kernels/darcy_afq_state.h"
 //%include "kernels/darcy_support_afq_state.h"
 //%include "kernels/darcy_afq_CG.h"

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
   * Extrawrappers for the Dofmap
   */
  void Dofmap_Get_np(dofmap_t * dm, hypervertex_t V,
		     int* DIM1, int** ARGOUTVIEW_ARRAY1) {
    int len;
    int *pay;
    len = Dofmap_Max_Len(dm);
    pay = (int*)malloc( sizeof(int)* (len));
    Dofmap_Get(dm,V, pay, DIM1);
    *ARGOUTVIEW_ARRAY1 = pay;
    // MEMLEAK!
  }
  void Dofmap_Get_List_np(dofmap_t * dm, int nvert, hypervertex_t * verts,
		     int* DIM1, int** ARGOUTVIEW_ARRAY1) {
    int len;
    int *pay;
    len = nvert * Dofmap_Max_Len(dm);
    pay = (int*)malloc(sizeof(int)* (len));
    Dofmap_Get_List(dm, nvert,verts, pay,DIM1);
    *ARGOUTVIEW_ARRAY1 = pay;
  }
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
   * Other wrappers
   */
  void Interpolate_np(int nu1, int du1, real_t * u1,
		      int nx1, int dx1, real_t * x1,
		      int nu2, int du2, real_t * u2,
		      int nx2, int dx2, real_t * x2,
		      real_t rad)
  {
    Interpolate( u1, x1, nu1,
		 u2, x2, nu2,
		 du1, dx1, rad);
  }
  
  /*
   * Assembly wrappers
   */

  void assemble_targets_np(PyObject * targetlist,
			   kernel_t * ke, hypergraph_t * hg,
			   PyObject * dofmaplist,
			   PyObject * datalist)
  {
    int i=0, isnewobj=0;
    PyObject *obj, *subobj;
    PyArrayObject * arrobj;

    if(!PyList_Check(targetlist)) return;
    if(!PyList_Check(datalist)) return;
    if(!PyList_Check(dofmaplist)) return;

    int ntarget = PyList_Size(targetlist);
    int ndata = PyList_Size(datalist);
    int ndofmap = PyList_Size(dofmaplist);

    
    dofmap_t * dofmaps[ndofmap];
    assemble_target_t att[ntarget];
    real_t * data_ptrs[ndata];
    
    int nnewobj = 0;
    PyArrayObject * newobjs[3*ntarget + ndata];
    /* Step 1: Build the target list */
    for(i=0;i<ntarget;i++) {
      obj = PyList_GetItem(targetlist, i);
      if(PyTuple_Check(obj)) {
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
      // Set the iterators just 'cuz
      att[i].IIiter = att[i].II;
      att[i].JJiter = att[i].JJ;
      att[i].Viter  = att[i].V;
    }
    
    /* Step 2: Collect the data ptrs */
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

    /* Step 3: Create the dofmap list */
    for(i=0;i<ndofmap;i++) {
      obj = PyList_GetItem(dofmaplist,i);
      // BEEN FIXED
      const int rest = SWIG_ConvertPtr(obj, (void**)(dofmaps+i),SWIGTYPE_p_dofmap_t, 0);
      //if (!SWIG_IsOK(res)) {
      //	SWIG_exception_fail(SWIG_ArgError(res), "error in dofmaptlist");	
      //}
    }
    /* Step 4: assemble! */
    assemble_targets(ke, hg,
		     dofmaps, data_ptrs,
		     att);

    /* Step 5: Decrease reference counts */
    for(i=0;i<nnewobj;i++) {
      Py_DECREF(newobjs[i]);
    }
    
  }

  void Tie_Cells_and_Particles_np(hypergraph_t * hgnew,
				  hypergraph_t * mesh,
				  kernel_t * ke_circum,
				  kernel_t * ke_centroid,
				  kernel_t * ke_inside,
				  PyObject * dofmaplist,
				  PyObject * datalist,
				  int Npart, int dim, real_t * x,
				  int nvert, hypervertex_t * verts
				  )
  {
    int i=0, isnewobj=0, nnewobj=0;
    PyObject *obj;
    PyArrayObject *arrobj;
    if(!PyList_Check(datalist)) return;
    if(!PyList_Check(dofmaplist)) return;
    int ndata = PyList_Size(datalist);
    int ndofmap = PyList_Size(dofmaplist);
    
    dofmap_t * dofmaps[ndofmap];
    real_t * data_ptrs[ndata];

    /* Step 2: Collect the data ptrs */
    PyArrayObject * newobjs[ndata];
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

    /* Step 3: Create the dofmap list */
    for(i=0;i<ndofmap;i++) {
      obj = PyList_GetItem(dofmaplist,i);
      // BEEN FIXED
      const int rest = SWIG_ConvertPtr(obj, (void**)(dofmaps+i),SWIGTYPE_p_dofmap_t, 0);
      //if (!SWIG_IsOK(res)) {
      //	SWIG_exception_fail(SWIG_ArgError(res), "error in dofmaptlist");	
      //}
    }

    if(nvert != Npart) {
      verts = NULL;
    }
    Tie_Cells_and_Particles(hgnew,mesh,
			    ke_circum,ke_centroid,ke_inside,
			    dofmaps,data_ptrs,
			    Npart,dim,x,
			    verts);
    /* Step 5: Decrease reference counts */
    for(i=0;i<nnewobj;i++) {
      Py_DECREF(newobjs[i]);
    }
  }
%}
