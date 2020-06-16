from dftpart import scfeda
from pyscf import lib

par = 8
with lib.with_omp_threads(par):
    test1 = scfeda.EDA()
    test1.gjf = 'w32.gjf'
    test1.method = ['hf','6-31gss']
    test1.output = 'w32'
    test1.verbose = 9
    #test1.build()
    test1.showinter=True
    test1.lso = 'w32/w32.lso'
    test1.kernel()
        

