from dftpart import gfea3
from pyscf import lib

par = 8
with lib.with_omp_threads(par):
    fr0 = gfea3.GFEA()
    fr0.inputstyle = 'frg'
    fr0.method = ['hf', '6-31gs','cart','charge']
    fr0.gjfname = 'w19'
    fr0.output = 'w19'
    fr0.verbose = 9
    fr0.showinter = True
    #fr0.do_deriv = True
    fr0.kernel()

