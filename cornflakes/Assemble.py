import numpy as np
import scipy.sparse
import cornflakes_library as cflib
from Hypergraph import Hypergraph


IndexMap = cflib.IndexMap
CFData = cflib.CFData
CFData_BC = cflib.CFData_BC
CFData_From_Ptr = cflib.CFData_From_Ptr
CFMat = cflib.CFMat
CFMat_BC = cflib.CFMat_BC

def Fill_Sparsity2(ke, H, data, cftargets):
    cflib.fill_sparsity2_np(ke,H.hg, data, cftargets)
    outps = cflib.outpArray_frompointer(ke.outp)
    onames = [ outps[j].name for j in xrange(ke.noutp) ]
    # for o in onames:
		# try:
		# cftargets[o][0].Finalize_Sparsity()
            # except AttributeError:
                # pass

def Assemble2(ke,H, data, cftargets, wipe=True, ndof=0):
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
    xAssemble(Kernal, Hypergraph,
    x         {'u':( u, dm_u ), 'p':( p,dm_p ), 'params':(params, dm_params) },
    x         {'R':( cfmat_R, (dm_u,dm_p)), 'K':( cfmat_K, (dm_u,dm_p)  ) })
    xor
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
    # Sanitize the output dictionary
	# from IPython import embed ; embed()
    outps = cflib.outpArray_frompointer(ke.outp)
    onames = [ outps[j].name for j in xrange(ke.noutp) ]
    need_to_sparsify = False
    made_new = False
    for  j in xrange(ke.noutp):
        name = outps[j].name
        try:
            cftargets[name]
        except KeyError:
            print "cornflakes runtime error: kernel ', ke.name,': You're missing key ", name, " in your target dict!"
            raise KeyError('kernel assembly error')
        # Do we need to make it for it?
        if not hasattr(cftargets[name][0],'Place'):
            made_new = True
            if outps[j].rank==2:
                cft = CFMat(ndof)
            else:
                cft = CFData(ndof)
            cftargets[name] = [ cft ] + list(cftargets[name])
            need_to_sparsify = True
    # from IPython import embed ; embed()
    if need_to_sparsify:
        cflib.fill_sparsity2_np(ke,H.hg, data, cftargets)
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
    cflib.assemble2_np(ke,H.hg, data, cftargets)

    if (wipe):
        for o in onames:
            cftargets[o][0].Finalize()
        if (made_new):
            #print "Returning a copy for ", ke.name
            # l = [ np.copy(cftargets[o][0].np()) for o in onames ]
            l = [ cftargets[o][0].np().copy() for o in onames ]
            return l
        else:
            #print "Returning the mask for ", ke.name
            return [ cftargets[o][0].np() for o in onames ]


# Deprecated call signature
def Filter(ke,H, dofmaps,data):
    htrue = Hypergraph()
    hfalse = Hypergraph()
    cflib.filter_np(ke,H.hg, dofmaps, data, htrue.hg, hfalse.hg)
    return htrue,hfalse

def Apply_BC(dofs,vals, K=None,R=None):
    if K is not None:
        for i in dofs:
            K.data[K.indptr[i]:K.indptr[i+1]] = 0.0
            K[i,i] = 1.0
    if R is not None:
        R[dofs]=vals

