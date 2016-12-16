import numpy as np
import scipy.sparse
import cornflakes_library as cflib
from Hypergraph import Hypergraph

class IndexMap():
    def __init__(self, start,end,bcs):
        self.imap = cflib.indexmap_t()
        cflib.IndexMap_New(self.imap, start,end,bcs)
    def __del__(self):
        cflib.IndexMap_Destroy(self.imap)

# The python version wraps in the BCs? IDK How I feel about this......
class CFData():
    def __init__(self, ndof, imap=None, fromptr=None):
        self.dat = cflib.cfdata_t()
        if imap is not None:
            self.dat_bc = cflib.cfdata_t()
            cflib.CFData_BC_New(self.dat_bc, self.dat, imap)
            if fromptr is not None:
                cflib.CFData_Default_New_From_Ptr(self.dat,fromptr)
            else:
                cflib.CFData_Default_New(self.dat,imap.Nsys)
                
        else:
            self.dat_bc = None
            if fromptr is not None:
                cflib.CFData_Default_New_From_Ptr(self.dat,fromptr)
            else:
                cflib.CFData_Default_New(self.dat, ndof)
    def __del__(self):
        print "CFDATA WAS DESTROYED"
        cflib.CFData_Destroy(self.dat)
        if(self.dat_bc):
            cflib.CFData_Destroy(self.dat_bc)
    def top(self):
        " Return the CFData that need to be placed into "
        if(self.dat_bc):
            return self.dat_bc
        else:
            return self.dat
    def Finalize(self):
        cflib.CFData_Finalize(self.dat)
    def Wipe(self):
        cflib.CFData_Wipe(self.dat)
    def np(self):
        " Wrap as a numpy array "
        return cflib.CFData_Default_View_np(self.dat)
    
class CFMat():
    def __init__(self, ndof, Rtarg=None,ubc=None,imap=None):
        self.mat = cflib.cfmat_t()

        if imap:
            self.ubc = ubc
            self.mat_bc = cflib.cfmat_t()
            cflib.CFMat_BC_New(self.mat_bc, self.mat, Rtarg.dat,self.ubc,imap)
            cflib.CFMat_CSR_New(self.mat,imap.Nsys)
        else:
            self.mat_bc = None
            cflib.CFMat_CSR_New(self.mat,ndof)
    def __del__(self):
        print "CFMAT WAS DESTROYED"
        cflib.CFMat_Destroy(self.mat)
        if(self.mat_bc):
            cflib.CFMat_Destroy(self.mat_bc)
    def top(self):
        " Return the cfmat that needs to be placed into "
        if(self.mat_bc):
            return self.mat_bc
        else:
            return self.mat
    def Finalize_Sparsity(self):
        cflib.CFMat_Finalize_Sparsity(self.mat)
    def Finalize(self):
        cflib.CFMat_Finalize(self.mat)
    def Wipe(self):
        cflib.CFMat_Wipe(self.mat)
    def np(self):
        " Wrap as a scipy csr matrix "
        I,J,V = cflib.CFMat_CSR_View_np(self.mat)
        return scipy.sparse.csr_matrix( (V,J,I) , shape=(self.mat.N, self.mat.N) )
    
class CFTargets():
    def __init__(self, ke,H, dofmaps, ndof, bcs=None,bcvals=None):
        # First, if there are BCs, create the indexmap
        if bcs is not None:
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
                K = CFMat(ndof,Rtarg,ubc,self.imap)
                cflib.Target_New_From_Ptr(targ,2,K.top())
                self.cfobjs.append(K)
            else:
                R = CFData(ndof,self.imap)
                if not Rtarg:
                    Rtarg = R
                cflib.Target_New_From_Ptr(targ,1,R.top())
                self.cfobjs.append(R)
            self.targets.append(targ)
            
        # Fill in the sparsities
        cflib.fill_sparsity_np(ke,H.hg, [ d.dm for d in dofmaps ], self.targets)
        
        # Finalize the sparsities
        for k in self.cfobjs:
            try:
                k.Finalize_Sparsity()
            except AttributeError:
                pass
        
    def __del__(self):
        if self.imap:
            cflib.IndexMap_Destroy(self.imap)
        # The targets don't need to be finalized since they don't own anything
        
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

def Fill_Sparsity(ke, H, dofmaps, cftargets):
    att = _sanitize_targets(cftargets)
    
    cflib.fill_sparsity_np(ke, H.hg, [d.dm for d in dofmaps], att)
    
def Assemble(ke,H, dofmaps,data, cftargets, wipe=True):
    att = _sanitize_targets(cftargets)
    
    # call wipe
    if(wipe):
        cftargets.Wipe()
    cflib.assemble_np(ke,H.hg, [ d.dm for d in dofmaps], data, att)
    if(wipe):
        cftargets.Finalize()
    # call finalize
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
    hfalse = Hypergraph()
    cflib.filter_np(ke,H.hg, [d.dm for d in dofmaps], data, htrue.hg, hfalse.hg)
    return htrue,hfalse

def Apply_BC(dofs,vals, K=None,R=None):
    if K!=None:
        for i in dofs:
            K.data[K.indptr[i]:K.indptr[i+1]] = 0.0
            K[i,i] = 1.0
    if R!=None:
        R[dofs]=vals
