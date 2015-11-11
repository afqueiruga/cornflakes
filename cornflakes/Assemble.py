import numpy as np
import scipy.sparse
import cornflakes_library as clfib

def Assemble_Targets(ranks, ke,hg, dofmap,ndof, data):
    #ranks=(1,2,2)
    len_loc_out = cflib.kernel_outp_len(ke,hg.l_edge)
    matsize = len_loc_out*len_loc_out*hg.n_edge
    allocator = {
        0:(lambda : np.zeros(1,dtype=np.double)),
        1:(lambda : np.zeros(ndof,dtype=np.double)),
        2:(lambda : (np.zeros(matsize,dtype=np.double),
                     np.zeros(matsize,dtype=np.intc),
                     np.zeros(matsize,dtype=np.intc)))
    }
    forms = [ allocator[r]() for r in ranks ]
    from IPython import embed
    #print "ayuuup"
    cflib.assemble_targets_np(forms,
                                  ke,hg,dofmap,
                                  data)
    #print "anoope"
    for i,r in enumerate(ranks):
        if r==2:
            Kcoo = scipy.sparse.coo_matrix((forms[i][0],(forms[i][1],forms[i][2])), (ndof,ndof))
            K = Kcoo.tocsr()
            forms[i]=K
    return forms

def Apply_BC(dofs,vals, K=None,R=None):
    if K!=None:
        for i in dofs:
            K[i,:] = 0.0
            K[i,i] = 1.0
    if R!=None:
        R[dofs]=vals



def Assemble_Matrix(ke,hg, dofmap,ndof, data):
    len_loc_out = cflib.kernel_outp_len(ke,hg.l_edge)
    II = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.intc)
    JJ = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.intc)
    KK = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.double)
    cflib.assemble_matrix_np(II,JJ,KK,
                              ke,hg,dofmap,
                              data)
    Kcoo = scipy.sparse.coo_matrix((KK,(II,JJ)),(ndof,ndof))
    
    K = Kcoo.tocsr()
    return K

def Assemble_Vector_Matrix(ke,hg, dofmap,ndof, data):
    R = np.zeros( ndof, dtype=np.double)

    len_loc_out = cflib.kernel_outp_len(ke,hg.l_edge)
    II = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.intc)
    JJ = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.intc)
    KK = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.double)
    cflib.assemble_vector_matrix_np(R,
                                 II,JJ,KK,
                                 ke,hg,dofmap,
                                 data)
    
    Kcoo = scipy.sparse.coo_matrix((KK,(II,JJ)),(ndof,ndof))
    
    K = Kcoo.tocsr()

    return R,K

