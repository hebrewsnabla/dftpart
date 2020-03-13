from pyscf import gto, dft, lib
from QCKit import gjf_kit
import time

gjf = 'c8.gjf' 
for par in [24,16,12,8,4,2,1]:
    t1 = time.time()
    with lib.with_omp_threads(par):
        mol = gto.Mole()
        mol.basis = 'cc-pvdz'
        mol.atom, coords, charges, mol.charge, mol.spin = gjf_kit.gjf_parser(gjf)
        mf = dft.RKS(mol)
        mf.xc = 'm06-2x'
        mf.kernel()
    t2 = time.time()
    print(par,t2-t1)
