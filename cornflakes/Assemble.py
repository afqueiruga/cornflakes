import numpy as np
import scipy.sparse
from . import cornflakes_library as cflib
from .Hypergraph import Hypergraph


IndexMap = cflib.IndexMap
CFData = cflib.CFData
CFData_BC = cflib.CFData_BC
CFData_From_Ptr = cflib.CFData_From_Ptr
CFMat = cflib.CFMat
CFMat_BC = cflib.CFMat_BC

def Collect(ke, edge, data):
    if not isinstance(data, dict):
        from itertools import chain
        data = dict(chain(*[f.items() for f in data]))
    # Check the input dictionary
    inps = cflib.inpArray_frompointer(ke.inp)
    for i in range(ke.ninp):
        name = inps[i].name
        try:
            data[name]
        except KeyError:
            print("cornflakes runtime error: kernel ", ke.name,": You're missing key ", name, " in your data dict!")
            raise KeyError('kernel assembly error')        
    return cflib.collect_np(ke,edge,data)
    
def Fill_Sparsity(ke, H, data, cftargets):
    cflib.fill_sparsity_np(ke,H.hg, data, cftargets)
    outps = cflib.outpArray_frompointer(ke.outp)
    onames = [ outps[j].name for j in range(ke.noutp) ]
    # for o in onames:
		# try:
		# cftargets[o][0].Finalize_Sparsity()
            # except AttributeError:
                # pass

def Assemble(ke,H, data, cftargets, ndof=0, wipe=True):
    """
    Assemble a kernel across a graph with specified input data.

    This is the meat-and-potatos of cornflakes. 
    Example Usage:
    --------------

    1) Passing preallocated cfmat/cfdata objects:
    Assemble(Kernal, Hypergraph,
             {'u':( u, dm_u ), 'p':( p,dm_p ), 'params':(params, dm_params) },
             {'R':( cfdat_R, dm_u), 'K':( cfmat_K, dm_u ) })
    2) Multiple dofmaps to the outputs (e.g. monolithic systems):
    Assemble(Kernal, Hypergraph,
             {'u':( u, dm_u ), 'p':( p,dm_p ), 'params':(params, dm_params) },
             {'R':( cfmat_R, dm_u,dm_p), 'K':( cfmat_K, dm_u,dm_p  ) })
    3) Allocate for me:
    R,K = Assemble(Kernal, Hypergraph,
             {'u':( u, dm_u ), 'p':( p,dm_p ), 'params':(params, dm_params) },
             {'R':(dm_u), 'K':( dm_u ) },ndof=1000)
    
    Always returns the output objects in the order that the kernel defines them.

    """
    if not isinstance(data, dict):
        from itertools import chain
        data = dict(chain(*[f.items() for f in data]))
    # Check the input dictionary
    inps = cflib.inpArray_frompointer(ke.inp)
    for i in range(ke.ninp):
        name = inps[i].name
        try:
            data[name]
        except KeyError:
            raise RuntimeError("Cornflakes runtime error: kernel ", ke.name,\
                               ": You're missing key ", name, " in your data dict.")
            
    # Sanitize the output dictionary
    outps = cflib.outpArray_frompointer(ke.outp)
    onames = [ outps[j].name for j in range(ke.noutp) ]
    need_to_sparsify = False
    made_new = False
    for  j in range(ke.noutp):
        name = outps[j].name
        try:
            cftargets[name]
        except KeyError:
            raise RuntimeError("Cornflakes runtime error: kernel ", ke.name,\
                               ": You're missing key ", name, " in your target dict")
        # Do we need to make it for it?
        if not hasattr(cftargets[name][0],'Place'):
            if ndof<=0:
                raise RuntimeError('Cornflakes runtime error: You want Assemble to allocate new outputs, but ndof was not set.')
            made_new = True
            if outps[j].rank==2:
                cft = CFMat(ndof)
            else:
                cft = CFData(ndof)
            cftargets[name] = [ cft ] + list(cftargets[name])
            need_to_sparsify = True
    if need_to_sparsify:
        cflib.fill_sparsity_np(ke,H.hg, data, cftargets)
        for o in onames:
            try:
                cftargets[o][0].Finalize_Sparsity()
            except AttributeError:
                pass
    
    # Wipe them if needed
    if(wipe):
        for o in onames:
            cftargets[o][0].Wipe()
    # Call the C routines
    cflib.assemble_np(ke,H.hg, data, cftargets)

    if (wipe):
        for o in onames:
            cftargets[o][0].Finalize()
        if (made_new):
            l = [ cftargets[o][0].np().copy() for o in onames ]
            return l
        else:
            return [ cftargets[o][0].np() for o in onames ]


        

# Deprecated call signature
def Filter(ke,H, data):
    # Check sanity of arguments
    if not isinstance(data, dict):
        from itertools import chain
        data = dict(chain(*[f.items() for f in data]))
    # Check the input dictionary
    inps = cflib.inpArray_frompointer(ke.inp)
    for i in range(ke.ninp):
        name = inps[i].name
        try:
            data[name]
        except KeyError:
            raise RuntimeError("Cornflakes runtime error: kernel ", ke.name,\
                               ": You're missing key ", name, " in your data dict.")

    # Make the new graphs and call C. No clean up required
    htrue = Hypergraph()
    hfalse = Hypergraph()
    cflib.filter_np(ke,H.hg, data, htrue.hg, hfalse.hg)
    return htrue,hfalse

def Apply_BC(dofs,vals, K=None,R=None):
    if K is not None:
        for i in dofs:
            K.data[K.indptr[i]:K.indptr[i+1]] = 0.0
            K[i,i] = 1.0
    if R is not None:
        R[dofs]=vals

