from dftpart.eda import mayer
from pyscf import lib

with lib.with_omp_threads(4):
    test1 = mayer.MayerEDA()
    test1.gjf = 'c3.gjf' # provide geom, charge, multiplicity
    test1.method = ['b3lyp','6-31gs'] 
    test1.output = 'c3'
    test1.verbose = 9
    #test1.build()
    #test1.showinter=True
    #test1.lso = 'c8b/c8b.lso'
    totE, conv, RR2 = test1.kernel()
        

