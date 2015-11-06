import numpy as np
import scipy.sparse
import mylibrary

def Assemble_Matrix(ke,hg, dofmap,ndof, data):
    len_loc_out = mylibrary.kernel_outp_len(ke,hg.l_edge)
    II = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.intc)
    JJ = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.intc)
    KK = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.double)
    mylibrary.assemble_matrix_np(II,JJ,KK,
                              ke,hg,dofmap,
                              data)
    Kcoo = scipy.sparse.coo_matrix((KK,(II,JJ)),(ndof,ndof))
    
    K = Kcoo.tocsr()
    return K

def Assemble_Vector_Matrix(ke,hg, dofmap,ndof, data):
    R = np.zeros( ndof, dtype=np.double)

    len_loc_out = mylibrary.kernel_outp_len(ke,hg.l_edge)
    II = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.intc)
    JJ = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.intc)
    KK = np.zeros( len_loc_out*len_loc_out*hg.n_edge , dtype=np.double)
    mylibrary.assemble_vector_matrix_np(R,
                                 II,JJ,KK,
                                 ke,hg,dofmap,
                                 data)
    
    Kcoo = scipy.sparse.coo_matrix((KK,(II,JJ)),(ndof,ndof))
    
    K = Kcoo.tocsr()

    return R,K

def Apply_BC(dofs,vals, K=None,R=None):
    if K!=None:
        for i in dofs:
            K[i,:] = 0.0
            K[i,i] = 1.0
    if R!=None:
        R[dofs]=vals
