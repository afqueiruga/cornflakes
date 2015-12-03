#
# IO Routines
#

def write_graph(fname, H, X, nodefields=None,edgefields=None):
    celltypekey = {
        1:1,
        2:3,
        3:5,
        4:9,
        #4:10,
        8:12}
    vecformatdict = {
        1:"{0} 0.0 0.0\n",
        2:"{0} {1} 0.0\n",
        3:"{0} {1} {2}\n"
        }
    elems = H.view()
    vecfmt = vecformatdict[X.shape[1]]
    
    fh = open(fname,"w")
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

    if nodefields:
        fh.write("POINT_DATA {0}\n".format(X.shape[0]))
        for n,f in nodefields.iteritems():
            fh.write("VECTORS {0} double\n".format(n))
            #fh.write("LOOKUP_TABLE default\n")
            for l in f:
                fh.write(vecfmt.format(*l))
    if edgefields:
        fh.write("CELL_DATA {0}\n".format(H.hg.n_edge))
        for n,f in edgefields.iteritems():
            fh.write("SCALARS {0} double\n".format(n))
            fh.write("LOOKUP_TABLE default\n")
            for l in f:
                fh.write("{0}\n".format(l))
    fh.close()

             
def write_silo(fname, H,X,cycle=0, time=0,nodefields=[], edgefields=[], PUTMESH=True, PUTCONN=True):
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
    
