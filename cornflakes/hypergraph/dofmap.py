import mylibrary

def make_dofmap(h, ke, Nnode):
    # Make the dofmap for where the kernel output will be
    # assembled into. It will look like
    # G1 ... E1 E3 ... N1 N2 N3
    # for now...
    # TODO: Do I need to give edges global indices??
    
    v = h.view()
    outlen = 0
    dofs = mylibrary.dofArray_frompointer(ke.outp)
    
    for i in xrange(ke.noutp):
        if dofs[i].loc == mylibrary.LOC_NODE:
            outlen+=dofs[i].len*h.hg.l_edge
        else:
            outlen+=dofs[i].len
    
    outmap = np.zeros((v.shape[0], outlen),dtype=np.intc)
    # Loop over the edges
    for i,e in enumerate(v):
        # Loop over the dofs
        off = 0
        for d in xrange(ke.noutp):
            if dofs[d].loc == mylibrary.LOC_NODE:
                for j,n in enumerate(e):
                    outmap[i,off+dofs[d].len*j:dofs[d].len*(j+1)] = range(n*dofs[d].len,(n+1)*dofs[d].len)
                off+= len(e)*dofs[d].len
            else:
                pass
    return outmap
