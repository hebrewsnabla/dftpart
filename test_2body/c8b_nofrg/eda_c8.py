from dftpart.eda import scfeda
from pyscf import lib

with lib.with_omp_threads(4):
    test1 = scfeda.EDA()
    test1.gjf = 'c8b_nofrg.gjf'
    test1.method = ['pbe0','6-31gss']
    test1.output = 'c8b_nofrg'
    test1.verbose = 9
    #test1.build()
    #test1.showinter=True
    #test1.lso = 'c8b/c8b.lso'
    test1.kernel()
        

