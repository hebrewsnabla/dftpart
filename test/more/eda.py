from dftpart.eda import mayer
from pyscf import lib

with lib.with_omp_threads(4):
    #for mol in ['c2h6', 'c2h4', 'c2h5oh']:
    for mol in [ 'c2h4']:
        test1 = mayer.MayerEDA()
        test1.gjf = mol+'.gjf' # provide geom, charge, multiplicity
        test1.method = ['b3lyp','6-31gs'] 
        test1.output = mol
        test1.verbose = 9
        #test1.build()
        #test1.showinter=True
        #test1.lso = 'c8b/c8b.lso'
        totE, conv, RR2 = test1.kernel()
            

