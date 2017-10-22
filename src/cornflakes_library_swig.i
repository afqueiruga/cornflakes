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
	 (int Nparty,  int dimy, real_t * y),
	 (int nx1, int dx1, real_t * x1),
	 (int nx2, int dx2, real_t * x2),
	 (int nu1, int du1, real_t * u1)
     };
%apply (int DIM1, real_t *IN_ARRAY1) {
    (int N, real_t * payload),
    (int Ncfbc, real_t* Acfbc),
    (int Norig, real_t* Aorig)
    }
%apply (int DIM1, int  DIM2, real_t * INPLACE_ARRAY2) {
  (int nu2, int du2, real_t * u2)
    };
%apply (int DIM1, int* INPLACE_ARRAY1) {(int dim_II, int * array_II),
                                        (int dim_JJ, int * array_JJ)};
%apply (int DIM1, real_t* INPLACE_ARRAY1) {
  (int dim_KK, real_t * array_KK)
    };
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
PyObject * make_np_copy_i(int len, int * arr) {
    npy_intp dims[1] = { len };
    PyObject * cpy   = PyArray_SimpleNewFromData(1,dims, NPY_INT, arr);
    PyObject * npret = PyArray_SimpleNew        (1,dims, NPY_INT);
    PyArray_CopyInto((PyArrayObject*)npret,(PyArrayObject*)cpy);
    Py_DECREF(cpy);
    return npret;
  }
%}

/*
 * IndexMap object transcription
 */
%inline %{
   /* IndexMap operations on numpy arrays */
  void IndexMap_Set_Values_np(indexmap_t * self,
			int Ncfbc, real_t *Acfbc,
			int Norig, real_t* Aorig )
  {
    cfdata_t cfbc, orig;
    CFData_Default_New_From_Ptr(&cfbc,Ncfbc,Acfbc);
    CFData_Default_New_From_Ptr(&orig,Norig,Aorig);
    IndexMap_Set_Values(self,&cfbc,&orig);
  }
  void IndexMap_Get_Values_np(indexmap_t * self,
			int Ncfbc, real_t *Acfbc,
			int Norig, real_t* Aorig )
  {
    cfdata_t cfbc, orig;
    CFData_Default_New_From_Ptr(&cfbc,Ncfbc,Acfbc);
    CFData_Default_New_From_Ptr(&orig,Norig,Aorig);
    IndexMap_Get_Values(self,&cfbc,&orig);
  }
  
  void IndexMap_Push_np(indexmap_t * self,
			int Ncfbc, real_t *Acfbc,
			int Norig, real_t* Aorig )
  {
    cfdata_t cfbc, orig;
    CFData_Default_New_From_Ptr(&cfbc,Ncfbc,Acfbc);
    CFData_Default_New_From_Ptr(&orig,Norig,Aorig);
    IndexMap_Push(self,&cfbc,&orig);
  }
  void IndexMap_Pull_np(indexmap_t * self,
			int Ncfbc, real_t *Acfbc,
			int Norig, real_t* Aorig )
  {
    cfdata_t cfbc, orig;
    CFData_Default_New_From_Ptr(&cfbc,Ncfbc,Acfbc);
    CFData_Default_New_From_Ptr(&orig,Norig,Aorig);
    IndexMap_Pull(self,&cfbc,&orig);
  }
%}
%extend IndexMap {
    IndexMap(index_t istart, index_t iend,
	     index_t *BCs, index_t NBC) {
      indexmap_t * imap = malloc(sizeof(indexmap_t));
      IndexMap_New(imap, istart,iend,BCs,NBC);
      return imap;
    }
    ~IndexMap() {
      IndexMap_Destroy($self);
      free($self);
    }
    index_t Get(index_t i);
    
    void Set_Values(cfdata_t * cfbc,
		    cfdata_t * orig);
    void Get_Values(cfdata_t * cfbc,
		    cfdata_t * orig);
    void Push(cfdata_t * cfbc, cfdata_t * orig);
    void Pull(cfdata_t * cfbc, cfdata_t * orig);

    void Set_Values_np(int Ncfbc, real_t *Acfbc,
		       int Norig, real_t* Aorig );
    void Get_Values_np(int Ncfbc, real_t *Acfbc,
		       int Norig, real_t* Aorig );
    void Push_np(int Ncfbc, real_t *Acfbc,
		 int Norig, real_t* Aorig );
    void Pull_np(int Ncfbc, real_t *Acfbc,
		 int Norig, real_t* Aorig );
};
%pythoncode %{
IndexMap = indexmap_t
%}


/*
 * CFData object transcription
 */
%inline %{
  void CFData_Default_View_np(cfdata_t * self,  int* NA, real_t** VA) {
    *VA = CFData_Default_Data(self);
    *NA = self->N;
  }
  cfdata_t * CFData_BC(cfdata_t * R, indexmap_t * imap) {
    cfdata_t * dat = malloc(sizeof(cfdata_t));
    CFData_BC_New(dat,R,imap);
    return dat;
  }
  cfdata_t * CFData_From_Ptr(int N, real_t * payload) {
    cfdata_t * dat = malloc(sizeof(cfdata_t));
    CFData_Default_New_From_Ptr(dat,N,payload);
    return dat;
  }
  
%}
%extend CFData {
  CFData(int N) {
    // Default Constructor
    cfdata_t * dat = malloc(sizeof(cfdata_t));
    CFData_Default_New(dat,N);
    return dat;
  }
  ~CFData() {
    printf("I'm getting destroyed with size %d",$self->N);
    CFData_Destroy($self);
    free($self);
  }
  void Get_Values( int ndof,int * dofs, real_t * vals);
  void Set_Values( int ndof, int *dofs, real_t * vals);
  void Scatter( real_t * src);
  void Copy( cfdata_t * src);
  real_t * Place( int n, int * dofs, real_t * vals);
  
  void Wipe();
  void Finalize();
  void Get_Ptr( real_t **ptr);
  void Release_Ptr( real_t **ptr);
  void Print();

  void Default_View_np(int* NA, real_t** VA);
};
%pythoncode %{
CFData = cfdata_t
CFData.np = CFData.Default_View_np
%}


/*
 * CFMat object transcription
 */
%inline %{
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
  cfmat_t * CFMat_BC(cfmat_t * K, cfdata_t * R, cfdata_t * u,
					 indexmap_t * map) {
    cfmat_t * mat = malloc(sizeof(cfmat_t));
    CFMat_BC_New(mat, K,R,u,map);
    return mat;
  }
  %}
%extend CFMat {
  CFMat(int N) {
    cfmat_t * mat = malloc(sizeof(cfmat_t));
    CFMat_CSR_New(mat,N);
    return mat;
  }
  ~CFMat() {
    CFMat_Destroy($self);
    free($self);
  }
  
  void Add_Sparsity(int n, int *dofs);
  void Finalize_Sparsity();
  real_t * Place(int n, int * dofs, real_t * vals);
  void Set_Value(int i, int j, real_t v);
  void Wipe();
  void Finalize();

  void CSR_View_np(int** IA, int* Ni,
				   int** JA, int* Jnnz,
				   int* Vnnz, real_t** VA);
 };
%pythoncode %{
CFMat = cfmat_t
def CFMat_np(mat):
  import scipy.sparse
  " Wrap as a scipy csr matrix "
  I,J,V = CFMat_CSR_View_np(mat)
  return scipy.sparse.csr_matrix( (V,J,I) , shape=(mat.N, mat.N) )
CFMat.np=CFMat_np
%}


/*
 * Dofmap object transcription
 */
%inline %{
  PyObject * Dofmap_Get_np(dofmap_t * dm, hypervertex_t V) {
    int len;
    int *pay, pdim;
    len = Dofmap_Max_Len(dm);
    pay = (int*)malloc( sizeof(int)* (len));
    Dofmap_Get(dm,V, pay, &pdim);
    
    // MEMLEAK!
    PyObject * npret = make_np_copy_i(len,pay);
    free(pay);
    return npret;
  }
  PyObject * Dofmap_Get_List_np(dofmap_t * dm, int nvert, hypervertex_t * verts) {
    int len;
    int *pay, pdim;
    len = nvert * Dofmap_Max_Len(dm);
    pay = (int*)malloc(sizeof(int)* (len));
    
    Dofmap_Get_List(dm, nvert,verts, pay, &pdim);

    PyObject * npret = make_np_copy_i(len,pay);
    free(pay);
    return npret;
  }
  dofmap_t * new_Dofmap_Tabled(int Nentry, int stride, int * table, int offset ) {
    dofmap_t * dmap = malloc(sizeof(dofmap_t));
    Dofmap_Tabled(dmap, Nentry,stride,table,offset);
    return dmap;
  }
  %}
%extend Dofmap {
  Dofmap(int stride, int offset) {
    dofmap_t * dmap = malloc(sizeof(dofmap_t));
    Dofmap_Strided(dmap, stride,offset);
    return dmap;
  }
  ~Dofmap() {
    Dofmap_Destroy($self);
    free($self);
  }
  
  void Get_List(int nvert, hypervertex_t * verts,int * dofs, int * ndofs);
  void Get(hypervertex_t V, int * dofs, int *ndofs);
  
  PyObject * Get_List_np(int nvert, hypervertex_t * verts);
  PyObject * Get_np(hypervertex_t V);
  
  int  Max_Len();
  void Destroy();
 };
%pythoncode %{
Dofmap = dofmap_t
Dofmap.Get = Dofmap.Get_np
Dofmap.Get_List = Dofmap.Get_List_np

def Dofmap_Strided(stride,offset=0):
  return Dofmap(stride,offset)

def Dofmap_From_Vertices(stride, vertices, offset=0):
  start = int(vertices.min())
  stop =  int(vertices.max())+1
  table = np.zeros(stop-start, dtype=np.intc)
  for i,l in enumerate(vertices):
    table[ l-start ] = i + offset
  return new_Dofmap_Tabled(table.reshape((table.size/stride,stride)), -start )
%}


/*
 * Hypergraph object transcription
 */
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
  %}
/* TODO: Hypergraph needs to be redone anyways
   %extend Hypergraph {
   Hypergraph(int alloc_init) {
   hypergraph_t * hg = malloc(sizeof(hypergraph_t));
   Hypergraph_Alloc(hg,alloc_init);
   return hg;
   }
   ~Hypergraph() {
   Hypergraph_Destroy($self);
   free($self);
   }
   Push_Edge(int l_edge, hypervertex_t * verts);
  
   };
   %pythoncode %{
   Hypergraph = hypergraph_t
   %}
*/

/*
 * Extra Wrappers
 */
%inline %{  
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
  void Interpolate_Closest_np(int nu1, int du1, real_t * u1,
							  int nx1, int dx1, real_t * x1,
							  int nu2, int du2, real_t * u2,
							  int nx2, int dx2, real_t * x2,
							  real_t rad)
  {
    Interpolate_Closest( u1, x1, nu1,
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
      const int rest = SWIG_ConvertPtr(obj, (void**)(dofmaps+i),SWIGTYPE_p_Dofmap, 0);
      //if (!SWIG_IsOK(res)) {
      //    SWIG_exception_fail(SWIG_ArgError(res), "error in dofmaptlist");    
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
      const int rest = SWIG_ConvertPtr(obj, (void**)(dofmaps+i),SWIGTYPE_p_Dofmap, 0);
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
				   PyObject * targetlist)
  {
    int i, isnewobj;
    PyObject *obj;
    PyArrayObject *arrobj;
    
    if(!PyList_Check(dofmaplist)) return;
    if(!PyList_Check(targetlist)) return;
    if(!PyList_Check(datalist))   return;
    
    /* Step 1: Create the dofmap list */
    int ndofmap = PyList_Size(dofmaplist);
    dofmap_t * dofmaps[ndofmap];
    for(i=0;i<ndofmap;i++) {
      obj = PyList_GetItem(dofmaplist,i);
      const int rest = SWIG_ConvertPtr(obj, (void**)(dofmaps+i),SWIGTYPE_p_Dofmap, 0);
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

    /* Step 3: Collect the data ptrs */
    int ndata = PyList_Size(datalist);
    cfdata_t  data[ndata];
    cfdata_t* data_ptrs[ndata];
    int nnewobj = 0;
    PyArrayObject * newobjs[ndata];
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
    for(i=0;i<ndata;i++) data_ptrs[i] = data+i;
    
    /* Step 4: Make the call */
    assemble(ke,hg, dofmaps,data_ptrs, att); 

    /* Step 5: Decrease reference counts */
    for(i=0;i<nnewobj;i++) {
      Py_DECREF(newobjs[i]);
    }
  }


  void assemble2_np(kernel_t * ke, hypergraph_t * hg,
		    PyObject * datadict,
		    PyObject * outpdict)
  {
    if(!PyDict_Check(datadict)) return;
    if(!PyDict_Check(outpdict)) return;

    /* We may need to allocate new numpy arrays */
    int isnewobj, n_newobj=0;
    PyArrayObject *newobjs[ke->ninp];;

    /* Extract the input signature from the data dictionary */
    cfdata_t data[ke->ninp];
    cfdata_t *data_ptrs[ke->ninp];
    dofmap_t * idofmaps[ke->ninp];
    for(int i=0; i<ke->ninp; i++) {
      PyObject *pair, *obj_dat, *obj_dm;
      PyArrayObject *arrobj;
      pair = PyDict_GetItemString(datadict, ke->inp[i].name);
	  
      /* Get the CFData */
      obj_dat = PySequence_GetItem(pair,0);
      isnewobj = 0;
      arrobj = obj_to_array_contiguous_allow_conversion(obj_dat,NPY_DOUBLE,&isnewobj);
      if(isnewobj) {
		newobjs[n_newobj] = arrobj;
		n_newobj++;
      }
      CFData_Default_New_From_Ptr(data+i, array_size(arrobj,0), array_data(arrobj));
      data_ptrs[i] = data+i;
      /* Get the dofmap */
      obj_dm  = PySequence_GetItem(pair,1);
      const int rest = SWIG_ConvertPtr(obj_dm, (void**)(idofmaps+i),SWIGTYPE_p_Dofmap, 0);
      Py_DECREF(obj_dat);
      Py_DECREF(obj_dm);
    }

    /* Extract the output signature from the data dicctionary.
       The python layer of Assemble is reponsible for initializing empty
       output targets and calling fill_sparsity if needed. */
    //target_t targets[ke->noutp];
    void * targets[ke->noutp];
    dofmap_t * odofmaps[ke->noutp*KERNEL_OUT_MAP_MAX]; // Yup, making a giant argument.
    for(int i=0; i<ke->noutp; i++) {
      PyObject *pair, *obj_targ, *seq_dms, *obj_dm;
      pair = PyDict_GetItemString(outpdict, ke->outp[i].name);
      /* Get the target */
      obj_targ = PySequence_GetItem(pair,0);
      if(ke->outp[i].rank==2) {
	const int rest = SWIG_ConvertPtr(obj_targ, (void**)(targets+i), SWIGTYPE_p_CFMat, 0);
      } else {
	const int rest = SWIG_ConvertPtr(obj_targ, (void**)(targets+i), SWIGTYPE_p_CFData, 0);
      }
      //targets[i].rank = t->rank;
      //if(targets[i].rank==2) targets[i].K = t->K;
      //else targets[i].R = t->R;
      /* Get the list of dofmaps */
      // seq_dms = PySequence_GetItem(pair,1);
      for(int j=0; j<ke->outp[i].nmap; j++) {
	obj_dm = PySequence_GetItem(pair,1 + j);
	const int rest = SWIG_ConvertPtr(obj_dm,
					 (void**)(odofmaps+i*KERNEL_OUT_MAP_MAX+j),
					 SWIGTYPE_p_Dofmap, 0);
	Py_DECREF(obj_dm);
      }
      Py_DECREF(obj_targ);
    }

    /* Make the call */
    assemble2(ke,hg,  data_ptrs, idofmaps, targets,odofmaps);
    
    /* Decrease reference counts. n_newobjs hasn't been observed to be >0 yet */
    for(int i=0;i<n_newobj;i++) {
      Py_DECREF(newobjs[i]);
    }
  }

  void filter2_np(kernel_t * ke, hypergraph_t * hg,
		    PyObject * datadict,
		    hypergraph_t * htrue, hypergraph_t * hfalse)
  {
	if(!PyDict_Check(datadict)) return;

    /* We may need to allocate new numpy arrays */
    int isnewobj, n_newobj=0;
    PyArrayObject *newobjs[ke->ninp];;

    /* Extract the input signature from the data dictionary */
    cfdata_t data[ke->ninp];
    cfdata_t *data_ptrs[ke->ninp];
    dofmap_t * idofmaps[ke->ninp];
    for(int i=0; i<ke->ninp; i++) {
      PyObject *pair, *obj_dat, *obj_dm;
      PyArrayObject *arrobj;
      pair = PyDict_GetItemString(datadict, ke->inp[i].name);
      /* Get the CFData */
      obj_dat = PySequence_GetItem(pair,0);
      isnewobj = 0;
      arrobj = obj_to_array_contiguous_allow_conversion(obj_dat,NPY_DOUBLE,&isnewobj);
      if(isnewobj) {
		newobjs[n_newobj] = arrobj;
		n_newobj++;
      }
      CFData_Default_New_From_Ptr(data+i, array_size(arrobj,0), array_data(arrobj));
      data_ptrs[i] = data+i;
      /* Get the dofmap */
      obj_dm  = PySequence_GetItem(pair,1);
      const int rest = SWIG_ConvertPtr(obj_dm, (void**)(idofmaps+i),SWIGTYPE_p_Dofmap, 0);
    }

    /* Make the call */
    printf("FINISH ME!\n");
    //filter2(ke,hg,  data_ptrs, idofmaps, targets,odofmaps);
    
    /* Decrease reference counts. n_newobjs hasn't been observed to be >0 yet */
    for(int i=0;i<n_newobj;i++) {
      Py_DECREF(newobjs[i]);
    }
  }  
  void fill_sparsity2_np(kernel_t * ke, hypergraph_t * hg,
						 PyObject * datadict,
						 PyObject * outpdict)
  {
    if(!PyDict_Check(datadict)) return;
    if(!PyDict_Check(outpdict)) return;
	
    /* We may need to allocate new numpy arrays */
    int isnewobj, n_newobj=0;
    PyArrayObject *newobjs[ke->ninp];;
	
    /* Extract the input signature from the data dictionary */
    cfdata_t data[ke->ninp];
    cfdata_t *data_ptrs[ke->ninp];
    dofmap_t * idofmaps[ke->ninp];
    for(int i=0; i<ke->ninp; i++) {
      PyObject *pair, *obj_dat, *obj_dm;
      PyArrayObject *arrobj;
      pair = PyDict_GetItemString(datadict, ke->inp[i].name);
      /* Get the CFData */
      obj_dat = PySequence_GetItem(pair,0);
      isnewobj = 0;
      arrobj = obj_to_array_contiguous_allow_conversion(obj_dat,NPY_DOUBLE,&isnewobj);
      if(isnewobj) {
		newobjs[n_newobj] = arrobj;
		n_newobj++;
      }
      CFData_Default_New_From_Ptr(data+i, array_size(arrobj,0), array_data(arrobj));
      data_ptrs[i] = data+i;
      /* Get the dofmap */
      obj_dm  = PySequence_GetItem(pair,1);
      const int rest = SWIG_ConvertPtr(obj_dm, (void**)(idofmaps+i),SWIGTYPE_p_Dofmap, 0);
      Py_DECREF(obj_dm);
      Py_DECREF(obj_dat);
    }

    /* Extract the output signature from the data dicctionary.
       The python layer of Assemble is reponsible for initializing empty
       output targets and calling fill_sparsity if needed. */
    //target_t targets[ke->noutp];
    void * targets[ke->noutp];
    dofmap_t * odofmaps[ke->noutp*KERNEL_OUT_MAP_MAX]; // Yup, making a giant argument.
    for(int i=0; i<ke->noutp; i++) {
      PyObject *pair, *obj_targ, *seq_dms, *obj_dm;
      pair = PyDict_GetItemString(outpdict, ke->outp[i].name);
      /* Get the target */
      obj_targ = PySequence_GetItem(pair,0);
      if(obj_targ==NULL) printf("WTF??\n");
      printf("%lx\n",obj_targ);
      if(ke->outp[i].rank==2) {
	const int rest = SWIG_ConvertPtr(obj_targ, (void**)(targets+i), SWIGTYPE_p_CFMat, 0);
	printf("K was %d\n",rest);
      } else {
	const int rest = SWIG_ConvertPtr(obj_targ, (void**)(targets+i), SWIGTYPE_p_CFData, 0);
	printf("R was %d\n",rest);
      }
      //targets[i].rank = t->rank;
      //if(targets[i].rank==2) targets[i].K = t->K;
      //else targets[i].R = t->R;
      /* Get the list of dofmaps */
      // seq_dms = PySequence_GetItem(pair,1);
      for(int j=0; j<ke->outp[i].nmap; j++) {
	obj_dm = PySequence_GetItem(pair,1 + j);
	if(obj_dm==NULL) printf("WTF??\n");
	const int rest = SWIG_ConvertPtr(obj_dm,
					 (void**)(odofmaps+i*KERNEL_OUT_MAP_MAX+j),
					 SWIGTYPE_p_Dofmap, 0);
	Py_DECREF(obj_dm);
      }
      printf("%lx : %d\n",obj_targ, obj_targ->ob_refcnt);
      Py_DECREF(obj_targ);
    }

    /* Make the call */
    fill_sparsity2(ke,hg, targets,odofmaps);
    
    /* Decrease reference counts. n_newobjs hasn't been observed to be >0 yet */
    if(n_newobj>0) printf("There were %d new objects in fill_sparsity2\n",n_newobj);
    for(int i=0;i<n_newobj;i++) {
      Py_DECREF(newobjs[i]);
    }
  }
  
  void filter_np(kernel_t * ke, hypergraph_t * hg,
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
      const int rest = SWIG_ConvertPtr(obj, (void**)(dofmaps+i),SWIGTYPE_p_Dofmap, 0);
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

  
  void Build_Proximity_Graph_Variable_np( hypergraph_t * hg,
					  int Npart, int dim, real_t * x,
					  int DIM1, real_t * IN_ARRAY1)
  {
    Build_Proximity_Graph_Variable(hg, Npart, dim, x, IN_ARRAY1);
  }
  void Build_Proximity_Graph_Given_Length_np(hypergraph_t * hg,
					     int Npart, int dim, real_t * x,
					     int N_desired, real_t cutoff,
					     int DIM1, real_t * INPLACE_ARRAY1)
  {
    Build_Proximity_Graph_Given_Length(hg,Npart,dim,x,N_desired,cutoff, INPLACE_ARRAY1);
  }
  void Build_Proximity_Graph_2Sets_Variable_np( hypergraph_t * hg,
						int Npart, int dim, real_t * x,
						int Nparty, int dimy, real_t * y,
						int DIM1, real_t * IN_ARRAY1)
  {
    Build_Proximity_Graph_2Sets_Variable(hg,
					 Npart, dim, x,
					 Nparty,dimy,y,
					 IN_ARRAY1);
  }
  void Build_Proximity_Graph_2Sets_Given_Length_np(hypergraph_t * hg,
						   int Npart, int dim, real_t * x,
						   int Nparty, int dimy, real_t * y,
						   int N_desired, real_t cutoff,
						   int DIM1, real_t * INPLACE_ARRAY1)
  {
    Build_Proximity_Graph_2Sets_Given_Length(hg,
					     Npart,dim,x,
					     Nparty,dimy,y,
					     N_desired,cutoff, INPLACE_ARRAY1);
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
      const int rest = SWIG_ConvertPtr(obj, (void**)(dofmaps+i),SWIGTYPE_p_Dofmap, 0);
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


  PyObject* Remove_Duplicate_Particles_np(int Npart, int dim, real_t * x,

				  real_t cutoff, real_t binsize) {
    int Naccept;
    real_t * y;

    y = malloc( dim*Npart * sizeof(real_t) );

    Remove_Duplicate_Particles(Npart,dim,x,
			     &Naccept,y,
			     cutoff,binsize);
    
    
    npy_intp dims_trim[2] = {Naccept, dim};
    // This should work, right? y is just malloc'd a little longer
    PyObject * y_long = PyArray_SimpleNewFromData(2,dims_trim, NPY_DOUBLE, y); 
    PyObject * y_trim = PyArray_SimpleNew(2,dims_trim, NPY_DOUBLE);
    PyArray_CopyInto((PyArrayObject*)y_trim,(PyArrayObject*)y_long);
    Py_DECREF(y_long);
    free(y);

    return y_trim;
  }
%}
