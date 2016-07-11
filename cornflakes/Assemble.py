import numpy as np
import scipy.sparse
import cornflakes_library as cflib
from Hypergraph import Hypergraph


def Assemble_Targets(ke,H, dofmaps,data, ndof):
    he = cflib.hyperedgesArray_frompointer(H.hg.he)
    n_edges = np.sum([ he[i].n_edge for i in xrange(H.hg.n_types) ])
    

    outps = cflib.outpArray_frompointer(ke.outp)
    forms = []
    for j in xrange(ke.noutp):
        op = outps[j]
        if op.rank==0:
            forms.append(np.zeros(1,dtype=np.double))
        elif op.rank==1:
            forms.append(np.zeros(ndof,dtype=np.double))
        else:
            matsize = 0
            for i in xrange(H.hg.n_types):
                oplen = cflib.kernel_outp_len(ke,op,he[i].l_edge)
                matsize += he[i].n_edge*oplen
            forms.append((np.zeros(matsize,dtype=np.double),
                          np.zeros(matsize,dtype=np.intc),
                          np.zeros(matsize,dtype=np.intc)))
    
    cflib.assemble_targets_np(forms, ke,H.hg, [ d.dm for d in dofmaps], data)
    
    for j in xrange(ke.noutp):
        if outps[j].rank==2:
            Kcoo = scipy.sparse.coo_matrix((forms[j][0],(forms[j][1],forms[j][2])), (ndof,ndof))
            K = Kcoo.tocsr()
            forms[j]=K
    return forms

def Filter(ke,H, dofmaps,data):
    htrue = Hypergraph()
    cflib.filter_np(ke,H.hg, [d.dm for d in dofmaps], data, htrue.hg)
    return htrue

def Apply_BC(dofs,vals, K=None,R=None):
    if K!=None:
        for i in dofs:
            K.data[K.indptr[i]:K.indptr[i+1]] = 0.0
            K[i,i] = 1.0
    if R!=None:
        R[dofs]=vals
