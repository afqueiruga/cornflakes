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

def _sanitize_targets(cftargets):
    # Need to figure out what type cftargets is?
    try: # Is it that stupid data structure up there?
        att = cftargets.targets
    except AttributeError:
        att = cftargets # It better be a list
    if type(att[0]) is not cflib.target_t: # Do I need to make list of target_t's?
        att2 = []
        for cf in att:
            targ = cflib.target_t()
            try:
                cf.mat
                rank = 2
            except AttributeError:           
                rank = 1
            cflib.Target_New_From_Ptr(targ,rank, cf.top() )
            att2.append(targ)
        att = att2
    return att


class CFTargets():
    def __init__(self, ke,H, dofmaps, ndof, bcs=None,bcvals=None):
        # First, if there are BCs, create the indexmap
        if False: #bcs is not None:
            self.imap = cflib.indexmap_t()
            cflib.IndexMap_New(self.imap, 0,ndof, bcs)
            ubc = cflib.cfdata_t()
            cflib.CFData_Default_New_From_Ptr(ubc,  bcvals.ravel())
        else:
            self.imap = None
            ubc = None
        # Initialize the data structures
        self.cfobjs = []
        self.targets = []
        outps = cflib.outpArray_frompointer(ke.outp)
        Rtarg = None
        # This loop requires that the R's come before the K's
        for j in xrange(ke.noutp):
            op = outps[j]
            targ = cflib.target_t()
            if op.rank==2:
                K = CFMat(ndof) #,Rtarg,ubc,self.imap)
                cflib.Target_New_From_Ptr(targ,2,K) #.top())
                self.cfobjs.append(K)
            else:
                R = CFData(ndof) #,self.imap)
                #if not Rtarg:
                #    Rtarg = R
                cflib.Target_New_From_Ptr(targ,1,R) #.top())
                self.cfobjs.append(R)
            self.targets.append(targ)
            
        # Fill in the sparsities
        cflib.fill_sparsity_np(ke,H.hg,  dofmaps, self.targets)
        
        # Finalize the sparsities
        for k in self.cfobjs:
            try:
                k.Finalize_Sparsity()
            except AttributeError:
                pass
        
    def np(self):
        return [ f.np() for f in self.cfobjs ]
    
    def Finalize(self):
        for f in self.cfobjs:
            f.Finalize()
    def Push(self, cfbc, orig):
        if self.imap:
            cflib.IndexMap_Set_Values_np(self.imap,cfbc,orig)
    def Pull(self, cfbc, orig):
        if self.imap:
            cflib.IndexMap_Get_Values_np(self.imap,cfbc,orig)
    def Wipe(self):
        for f in self.cfobjs:
            f.Wipe()

def Fill_Sparsity(ke, H, dofmaps, cftargets):
    att = _sanitize_targets(cftargets)    
    cflib.fill_sparsity_np(ke, H.hg, dofmaps, att)
    
def Assemble(ke,H, dofmaps,data, cftargets=None, wipe=True,ndof=0):
    if cftargets!=None:
        att = _sanitize_targets(cftargets)
        ret_np = False
    else:
        ret_np = True
        cftargets = CFTargets(ke,H,dofmaps,ndof)
        wipe = True
        att = cftargets.targets
    # call wipe
    if(wipe):
        cftargets.Wipe()
    cflib.assemble_np(ke,H.hg, dofmaps, data, att)
    if(wipe):
        cftargets.Finalize()
    if ret_np:
        return [ _.copy() for _ in cftargets.np() ]
    else:
        return cftargets.cfobjs

def Filter(ke,H, dofmaps,data):
    htrue = Hypergraph()
    hfalse = Hypergraph()
    cflib.filter_np(ke,H.hg, dofmaps, data, htrue.hg, hfalse.hg)
    return htrue,hfalse

def Apply_BC(dofs,vals, K=None,R=None):
    if K!=None:
        for i in dofs:
            K.data[K.indptr[i]:K.indptr[i+1]] = 0.0
            K[i,i] = 1.0
    if R!=None:
        R[dofs]=vals



#
# Deprecated
#
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
    
    cflib.assemble_targets_np(forms, ke,H.hg, dofmaps, data)
    
    for j in xrange(ke.noutp):
        if outps[j].rank==2:
            Kcoo = scipy.sparse.coo_matrix((forms[j][0],(forms[j][1],forms[j][2])), (ndof,ndof))
            K = Kcoo.tocsr()
            forms[j]=K
    return forms
