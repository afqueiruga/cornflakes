import numpy as np
import mylibrary as ml

def write_graph(fname, H, X):
    celltypekey = {
        1:1,
        2:3,
        3:5,
        4:9,
        #4:10,
        8:12}
    elems = H.view()
    fh = open(fname,"w")
    fh.write("# vtk DataFile Version 2.0\nGraph connectivity\nASCII\n")
    fh.write("DATASET UNSTRUCTURED_GRID\n")
    fh.write("POINTS {0} double\n".format(X.shape[0]))
    print X.shape[1]
    if X.shape[1]==1:
        for pt in X:
            fh.write("{0} 0.0 0.0\n".format(*pt))
    elif X.shape[1]==2:
        for pt in X:
            fh.write("{0} {1} 0.0\n".format(*pt))
    else:
        for pt in X:
            fh.write("{0} {1} {2}\n".format(*pt))
    print elems.shape
    fh.write("\nCELLS {0} {1}\n".format(elems.shape[0],elems.shape[0]*(1+elems.shape[1]))) # I assume they're all the same
    for el in elems:
        fh.write("{0} ".format(len(el))+" ".join([str(x) for x in el])+"\n")
    fh.write("\nCELL_TYPES {0}\n".format(elems.shape[0]))
    for el in elems:
        fh.write("{0}\n".format(celltypekey[elems.shape[1]]))
    fh.close()
    
