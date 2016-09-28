%module cornflakes_library
%{
#define SWIG_FILE_WITH_INIT
#include "cornflakes.h"
%}

%include "numpy.i"
%init %{
  import_array();
  //import_managed();
%}

%apply (int DIM1, int* IN_ARRAY1) {
         (int nvert, hypervertex_t * verts),
         (int l_edge, hypervertex_t * verts)
     };
%apply (int* IN_ARRAY1, int DIM1) {
         (index_t *BCs, index_t NBC)
     };
%apply (int DIM1, int DIM2, real_t* IN_ARRAY2) {
         (int Npart, int dim, real_t * x),
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

%apply ( int* DIM1, real_t** ARGOUTVIEW_ARRAY1 ) {
  (int* NA, real_t** VA),
    (int* Vnnz, real_t** VA) 
};

%apply ( int** ARGOUTVIEW_ARRAY1, int* DIM1 ) {
    (int** JA, int* Jnnz),
      (int** IA, int* Ni) // for some reason calling the array I makes a bug in swig. Hygiene?
      };

%include "carrays.i"
%array_class(int,intArray)
%array_class(k_map_t,k_mapArray)
%array_class(inp_t,inpArray)
%array_class(outp_t,outpArray)
%array_class(hyperedges_t,hyperedgesArray)
%array_class(target_t, targetArray)

%include "cpointer.i"
%pointer_class(int, intp)

%include "hypergraph.h"
%include "spatialhash.h"
%include "graphers.h"
%include "kernel.h"
%include "assemble.h"
%include "dofmap.h"
%include "util.h"
%include "target.h"
%include "cfdata_default.h"
%include "cfdata_bc.h"
%include "cfmat_default.h"
%include "cfmat_csr.h"
%include "cfmat_bc.h"

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
   * Wrappers for CFData and CFMat viewing
   */
  void CFData_Default_View_np(cfdata_t * self,  int* NA, real_t** VA) {
    *VA = CFData_Default_Data(self);
    *NA = self->N;
  }
  void CFMat_CSR_View_np(cfmat_t * self,
			 
			 int** IA, int* Ni,
			 int** JA, int* Jnnz,
			 int* Vnnz, real_t** VA)
  {
    
    *Ni  = self->N+1;
    *IA    = CFMat_CSR_Data(self)->IA;
    
    *Jnnz = CFMat_CSR_Data(self)->nnz;
    *JA    = CFMat_CSR_Data(self)->JA;
    
    *Vnnz = CFMat_CSR_Data(self)->nnz;
    *VA    = CFMat_CSR_Data(self)->V;
    
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
    target_t  att[ntarget];
    cfdata_t  data[ndata];
    cfdata_t* data_ptrs[ndata];
    int nnewobj = 0;
    PyArrayObject * newobjs[3*ntarget + ndata];
    /* Step 1: Build the target list */
    for(i=0;i<ntarget;i++) {
      obj = PyList_GetItem(targetlist, i);
      if(PyTuple_Check(obj)) {
	real_t *VV;
	int *II,*JJ;
	// 0: Get KK
	subobj = PyTuple_GetItem(obj,0);
	isnewobj = 0;
	arrobj = obj_to_array_contiguous_allow_conversion(subobj,NPY_DOUBLE,&isnewobj);
        VV = array_data(arrobj);
	if(isnewobj) { newobjs[nnewobj] = arrobj; nnewobj++; }
	// 1: Get II
	subobj = PyTuple_GetItem(obj,1);
	isnewobj = 0;
	arrobj = obj_to_array_contiguous_allow_conversion(subobj,NPY_INT,&isnewobj);
	II = array_data(arrobj);
	if(isnewobj) { newobjs[nnewobj] = arrobj; nnewobj++; }
	// 2: Get JJ
	subobj = PyTuple_GetItem(obj,2);
	isnewobj = 0;
	arrobj = obj_to_array_contiguous_allow_conversion(subobj,NPY_INT,&isnewobj);
	JJ = array_data(arrobj);

	cfmat_t * mat = malloc(sizeof(cfmat_t));
	CFMat_Default_From_Array(mat,-1/*don't know N!*/, VV,II,JJ);
	Target_New_From_Ptr(att+i,2,mat);
	if(isnewobj) { newobjs[nnewobj] = arrobj; nnewobj++; }
      }
      else { //  obj is (BETTER BE) a ndarray
	isnewobj = 0;
	arrobj = obj_to_array_contiguous_allow_conversion(obj,NPY_DOUBLE,&isnewobj);
	if(isnewobj) { newobjs[nnewobj] = arrobj; nnewobj++; }
	cfdata_t * dat = malloc(sizeof(cfdata_t));
	CFData_Default_New_From_Ptr(dat, array_size(arrobj,0), array_data(arrobj) );
	Target_New_From_Ptr(att+i,1,dat);


      }

    }
    
    /* Step 2: Collect the data ptrs */
    for(i=0;i<ndata;i++) {
      obj = PyList_GetItem(datalist,i );
      isnewobj = 0;
      arrobj = obj_to_array_contiguous_allow_conversion(obj,NPY_DOUBLE,&isnewobj);
      CFData_Default_New_From_Ptr(data+i, array_size(arrobj,0),  array_data(arrobj));
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

    /* Fill up the ptr array */
    for(i=0;i<ndata;i++) data_ptrs[i] = data+i;
    /* Step 4: assemble! */
    assemble(ke, hg, dofmaps, data_ptrs, att);

    /* Step 5: Decrease reference counts */
    for(i=0;i<nnewobj;i++) {
      Py_DECREF(newobjs[i]);
    }
    /* Step 6: Free the target data structures */
    for(i=0;i<ntarget;i++) {
      //Target_Destroy(att+i);
      if(att[i].rank==2) {
	free(att[i].K);
      } else {
	free(att[i].R);
      }
    }
  }
  void fill_sparsity_np(kernel_t * ke, hypergraph_t * hg,
			PyObject * dofmaplist,
			PyObject * targetlist)
  {
    int i;
    PyObject *obj;
    if(!PyList_Check(dofmaplist)) return;
    if(!PyList_Check(targetlist)) return;
    /* Step 1: Create the dofmap list */
    int ndofmap = PyList_Size(dofmaplist);
    dofmap_t * dofmaps[ndofmap];
    for(i=0;i<ndofmap;i++) {
      obj = PyList_GetItem(dofmaplist,i);
      const int rest = SWIG_ConvertPtr(obj, (void**)(dofmaps+i),SWIGTYPE_p_dofmap_t, 0);
    }

    /* Step 2: Create the target array */
    int ntarget = PyList_Size(targetlist);
    target_t att[ntarget];
    for(i=0;i<ntarget;i++) {
      obj = PyList_GetItem(targetlist,i);
      target_t * t;
      const int rest = SWIG_ConvertPtr(obj, (void**)(&t),SWIGTYPE_p_target_t, 0);
      att[i].rank = t->rank;
      if(att[i].rank==2)
	att[i].K = t->K;
      else
	att[i].R = t->R;
    }

    /* Step 3: Make the call */
    fill_sparsity(ke,hg, dofmaps, att);
  }
  void assemble_np(kernel_t * ke, hypergraph_t * hg,
		   PyObject * dofmaplist,
		   PyObject * datalist,
		   target_t * att)
  {
    
  }
  void filter_np(
		 kernel_t * ke, hypergraph_t * hg,
		 PyObject * dofmaplist,
		 PyObject * datalist,
		 hypergraph_t * htrue, hypergraph_t * hfalse)
  {
    int i=0, isnewobj=0;
    PyObject *obj, *subobj;
    PyArrayObject * arrobj;

    if(!PyList_Check(datalist)) return;
    if(!PyList_Check(dofmaplist)) return;

    int ndata = PyList_Size(datalist);
    int ndofmap = PyList_Size(dofmaplist);

    
    dofmap_t * dofmaps[ndofmap];
    cfdata_t  data[ndata];
    cfdata_t* data_ptrs[ndata];
    int nnewobj = 0;
    PyArrayObject * newobjs[ndata];
    
    /* Step 1: Build the target list */
    
    /* Step 2: Collect the data ptrs */
    for(i=0;i<ndata;i++) {
      obj = PyList_GetItem(datalist,i );
      isnewobj = 0;
      arrobj = obj_to_array_contiguous_allow_conversion(obj,NPY_DOUBLE,&isnewobj);
      CFData_Default_New_From_Ptr(data+i, array_size(arrobj,0),  array_data(arrobj));
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

    /* Fill up the ptr array */
    for(i=0;i<ndata;i++) data_ptrs[i] = data+i;
    /* Step 4: assemble! */
    filter(ke, hg, dofmaps, data_ptrs, htrue, hfalse);

    /* Step 5: Decrease reference counts */
    for(i=0;i<nnewobj;i++) {
      Py_DECREF(newobjs[i]);
    }
    /* Step 6: Free the target data structures */
    
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
    cfdata_t data[ndata];
    cfdata_t *data_ptrs[ndata];

    /* Step 2: Collect the data ptrs */
    PyArrayObject * newobjs[ndata];
    for(i=0;i<ndata;i++) {
      obj = PyList_GetItem(datalist,i );
      isnewobj = 0;
      arrobj = obj_to_array_contiguous_allow_conversion(obj,NPY_DOUBLE,&isnewobj);
      CFData_Default_New_From_Ptr(data+i, array_size(arrobj,0),array_data(arrobj));
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
    /* Fill up the ptr array */
    for(i=0;i<ndata;i++) data_ptrs[i] = data+i;
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


  /* Call the C version. No vargs supported, and only one hg will be read */
  PyObject* load_gmsh_np(int gdim_set,
		    hypergraph_t * hg,
		    char * fname)
  {
    real_t *xc;
    int nc;
    load_gmsh(&xc, &nc, gdim_set,
	      &hg, fname);

    npy_intp dims[2] = { nc, gdim_set };
    PyObject* x_wrap = PyArray_SimpleNewFromData(2,dims, NPY_DOUBLE, xc);
    PyObject * x_np = PyArray_SimpleNew(2,dims, NPY_DOUBLE);
    PyArray_CopyInto((PyArrayObject*)x_np,(PyArrayObject*)x_wrap);
    Py_DECREF(x_wrap);
    free(xc);
    return x_np;
  }
%}
