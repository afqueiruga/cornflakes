from Hypergraph import Hypergraph
import cornflakes_library as cflib

import os, errno
import numpy as np
#
# IO Routines
#
def make_path_and_open(fname,*args,**kwargs):
    path=os.path.dirname(fname)
    if path:
        try:
            os.makedirs(path)
        except OSError as exc: # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else: raise
    fh = open(fname,*args,**kwargs)
    return fh

def write_graph(fname, H, X, nodefields=None,edgefields=None):
    """Writes to a vtk file."""
    celltypekey = {
        1:1,
        2:3,
        3:5,
        4:9 if X.shape[1]==2 else 10,
        #4:10,
        8:12}
    vecformatdict = {
        1:"{0} 0 0\n",
        2:"{0} {1} 0\n",
        3:"{0} {1} {2}\n"
        }
    tenformatdict = {
        1:"{0} 0 0\n0 0 0\n0 0 0\n",
        2:"{0} {1} 0\n{2} {3} 0\n0 0 0\n",
        3:"{0} {1} {2}\n{3} {4} {5}\n{6} {7} {8}\n"
        }
    # yes, I know, eyerollemoji
    variformatdict = {
        1:"{0} 0 0\n0 0 0\n0 0 0\n",
        2:"{0} {1} 0\n0 0 0\n0 0 0\n",
        3:"{0} {1} {2}\n0 0 0\n0 0 0\n",
        4:"{0} {1} {2}\n{3} 0 0\n0 0 0\n",
        5:"{0} {1} {2}\n{3} {4} 0\n0 0 0\n",
        6:"{0} {1} {2}\n{3} {4} {5}\n0 0 0\n",
        7:"{0} {1} {2}\n{3} {4} {5}\n{6} 0 0\n",
        8:"{0} {1} {2}\n{3} {4} {5}\n{6} {7} 0\n",
        9:"{0} {1} {2}\n{3} {4} {5}\n{6} {7} {8}\n"
        }
        
    elems = H.view()[0]
    vecfmt = vecformatdict[X.shape[1]]
    tenfmt = tenformatdict[X.shape[1]]
    fh = make_path_and_open(fname,"w")
    fh.write("# vtk DataFile Version 2.0\nGraph connectivity\nASCII\n")
    fh.write("DATASET UNSTRUCTURED_GRID\n")
    
    fh.write("POINTS {0} double\n".format(X.shape[0]))
    for pt in X:
        fh.write(vecfmt.format(*pt))
    #if X.shape[1]==1:
    #    for pt in X:
    #        fh.write("{0} 0.0 0.0\n".format(*pt))
    #elif X.shape[1]==2:
    #    for pt in X:
    #        fh.write("{0} {1} 0.0\n".format(*pt))
    #else:
    #    for pt in X:
    #        fh.write("{0} {1} {2}\n".format(*pt))
    
    fh.write("\nCELLS {0} {1}\n".format(elems.shape[0],elems.shape[0]*(1+elems.shape[1]))) # I assume they're all the same
    for el in elems:
        fh.write("{0} ".format(len(el))+" ".join([str(x) for x in el])+"\n")
    fh.write("\nCELL_TYPES {0}\n".format(elems.shape[0]))
    for el in elems:
        fh.write("{0}\n".format(celltypekey[elems.shape[1]]))

    # Macro to write a data block
    def PUTFIELD(n,f):
        if len(f.shape)==1 or f.shape[1]==1:
            fh.write("SCALARS {0} double\n".format(n))
            fh.write("LOOKUP_TABLE default\n")
            for l in f:
                fh.write(str(l)+"\n")
        elif f.shape[1]==X.shape[1]:
            fh.write("VECTORS {0} double\n".format(n))
            for l in f:
                fh.write(vecfmt.format(*l))
        else:
            fh.write("TENSORS {0} double\n".format(n))
            fmt=variformatdict[len(f[0])]
            for l in f:
                fh.write(fmt.format(*l))
    # Dump all of the node fields
    if nodefields:
        fh.write("POINT_DATA {0}\n".format(X.shape[0]))
        for n,f in nodefields:
            PUTFIELD(n,f)
    # Cell fields now
    if edgefields:
        fh.write("CELL_DATA {0}\n".format(elems.shape[0]))
        for n,f in edgefields:
            PUTFIELD(n,f)
            
    fh.close()

             
def write_silo(fname, H,X,cycle=0, time=0,nodefields=[], edgefields=[], PUTMESH=True, PUTCONN=True):
    from pyvisfile.silo import SiloFile, IntVector, DB_ZONETYPE_BEAM, DB_NODECENT, DB_ZONECENT, DBOPT_CYCLE, DBOPT_DTIME, DBOPT_TIME, DB_CLOBBER
    
    silo = SiloFile(fname, mode=DB_CLOBBER)

    pair_edges = H.view()[0]
    zonelist_name = "foo_zonelist"
    nodelist = IntVector()
    nodelist.extend( int(i) for i in pair_edges[:,0:2].flat)
    shapetypes = IntVector()
    shapetypes.append(DB_ZONETYPE_BEAM)
    shapesizes = IntVector()
    shapesizes.append(2)
    shapecounts = IntVector()
    shapecounts.append(len(pair_edges))
    #from IPython import embed
    #embed()
    if PUTCONN:
        silo.put_zonelist_2(zonelist_name, len(pair_edges), 2, nodelist,
                            0,0, shapetypes, shapesizes, shapecounts)
    if PUTMESH:
        silo.put_ucdmesh("foo", [],
                     np.asarray(X.T,order="C"), len(pair_edges),
                     zonelist_name, None)
    def putvar(n,fo,LOC):
        if len(f.shape)==1 or f.shape[1]==1:
            silo.put_ucdvar1("node_"+n,"foo",
                             np.asarray(f,order="C",dtype=np.double),
                             LOC, {DBOPT_CYCLE:cycle,DBOPT_DTIME:float(time),DBOPT_TIME:float(time)})
        elif f.shape[1]==2:
            silo.put_ucdvar("node_"+n,"foo", [n+"x",n+"y"],
                            np.asarray(f.T,order="C",dtype=np.double),
                            LOC, {DBOPT_CYCLE:cycle,DBOPT_DTIME:float(time),DBOPT_TIME:float(time)})
        else:
            silo.put_ucdvar("node_"+n,"foo", [n+"x",n+"y",n+"z"],
                            np.asarray(f.T,order="C",dtype=np.double),
                            LOC,  {DBOPT_CYCLE:cycle,DBOPT_DTIME:float(time),DBOPT_TIME:float(time)})
   
    for n,f in nodefields:
        putvar(n,f,DB_NODECENT)
    for n,f in edgefields:
        putvar(n,f,DB_ZONECENT)

    silo.close()

    

def write_silo_meshfile(fname, H,X):
    from pyvisfile.silo import SiloFile, IntVector, DB_ZONETYPE_BEAM, DB_NODECENT, DB_ZONECENT, DBOPT_CYCLE, DBOPT_DTIME, DBOPT_TIME, DB_CLOBBER
    import numpy as np
    silo = SiloFile(fname, mode=DB_CLOBBER)

    pair_edges = H.view()[0]
    zonelist_name = "foo_zonelist"
    nodelist = IntVector()
    nodelist.extend( int(i) for i in pair_edges[:,0:2].flat)
    shapetypes = IntVector()
    shapetypes.append(DB_ZONETYPE_BEAM)
    shapesizes = IntVector()
    shapesizes.append(2)
    shapecounts = IntVector()
    shapecounts.append(len(pair_edges))
    silo.put_zonelist_2(zonelist_name, len(pair_edges), 2, nodelist,
                            0,0, shapetypes, shapesizes, shapecounts)
    silo.put_ucdmesh("foo", [],
                     np.asarray(X.T,order="C"), len(pair_edges),
                     zonelist_name, None)
    silo.close()
    
def write_silo_datfile(fname,mname,cycle=0, time=0, nodefields=[], edgefields=[]):
    from pyvisfile.silo import SiloFile, IntVector, DB_ZONETYPE_BEAM,\
        DB_NODECENT, DB_ZONECENT, DBOPT_CYCLE, DBOPT_DTIME, DBOPT_TIME, DB_CLOBBER
    from pyvisfile.silo import DBObjectType as DBOBjectType
    import numpy as np
    silo = SiloFile(fname, mode=DB_CLOBBER)
    silo.put_multimesh('foo', [(mname+":foo",DBOBjectType.DB_UCDMESH)])
    def putvar(n,fo,LOC):
        if len(f.shape)==1 or f.shape[1]==1:
            silo.put_ucdvar1("node_"+n,"foo",
                             np.asarray(f.T,order="C",dtype=np.double),
                             LOC, {DBOPT_CYCLE:cycle,DBOPT_DTIME:float(time),DBOPT_TIME:float(time)})
        elif f.shape[1]==2:
            silo.put_ucdvar("node_"+n,"foo", [n+"x",n+"y"],
                            np.asarray(f.T,order="C",dtype=np.double),
                            LOC, {DBOPT_CYCLE:cycle,DBOPT_DTIME:float(time),DBOPT_TIME:float(time)})
        else:
            silo.put_ucdvar("node_"+n,"foo", [n+"x",n+"y",n+"z"],
                            np.asarray(f.T,order="C",dtype=np.double),
                            LOC,  {DBOPT_CYCLE:cycle,DBOPT_DTIME:float(time),DBOPT_TIME:float(time)})
   
    for n,f in nodefields:
        putvar(n,f,DB_NODECENT)
    for n,f in edgefields:
        putvar(n,f,DB_ZONECENT)
    silo.close()
    


def write_gmsh_file(fname, H,X):
    nodeformatdict = {
        1:"{0} {1} 0.0 0.0\n",
        2:"{0} {1} {2} 0.0\n",
        3:"{0} {1} {2} {3}\n"
    }
    vecfmt = nodeformatdict[X.shape[1]]
    
    fh = make_path_and_open(fname,"w")
    fh.write("$MeshFormat\n2.2 0 8\n&EndMeshFormat\n")

    fh.write("$Nodes\n")
    fh.write("{0}\n".format(X.shape[0]))
    for i,l in enumerate(X):
      fh.write(vecfmt.format(1+i,*l))
    fh.write("$EndNodes\n")

    etype = 3 # TODO:
    # 1 = 2-node lines
    # 2 = 3-node triangle
    # 3 = 4-node quad
    # 4 = 4-node tet
    # 6 = 8-node hex
    elems = H.view()[0]
    fh.write("$Elements\n")
    fh.write("{0}\n".format(elems.shape[0]))
    for i,e in enumerate(elems):
        fh.write("{0} {1} 2 0 0".format(i+1,etype))
        for v in e: fh.write(" {0}".format(v+1))
        fh.write("\n")
    fh.write("$EndElements\n")
    fh.close()



def Load_gmsh(fname,gdim=3):
    """
    Load a gmsh mesh. Requires version 2 ascii format.
    
    Specify gdim to weed-out irrelevant elements.
    """
    # check if file exists and that it's the right version
    try:
        with open(fname,'r') as f:
            l = f.readline()
            if not "$MeshFormat" in l:
                raise RuntimeError("cornflakes: Not a gmsh file.")
            l = f.readline()
            if not l[0]=='2':
                raise RuntimeError("cornflakes: Not version 2 of gmsh mesh.")
    except Exception as e:
        raise e
    # Swig wrapper initialized the hypergraph, but the C routine assumes it 
    # wasn't initialized.
    H = Hypergraph()
    cflib.Hypergraph_Destroy(H.hg)
    X = cflib.load_gmsh_np(gdim, H.hg, fname)
    return X,H
    
    
def write_cloud(fname, X, nodefields):
    """Write a pointcloud to a vtk file without a graph."""
    Hcloud = Hypergraph()
    for l in xrange(X.shape[0]): Hcloud.Push_Edge([l])
    write_graph(fname, Hcloud, X,nodefields)
