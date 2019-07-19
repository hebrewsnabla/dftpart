from pyscf import gto,dft
import numpy as np
#import gen_grid_sep
#import numint_sep
import xceda
import os
import time

name = 'C12H26'

#baslist = ['6-31+gss','6-31++gss','def2dzvp','def2tzvp','cc-pvdz','cc-pvtz','aug-cc-pvtz']
#xclist = [] 

mol = gto.Mole()
with open('c12h26.xyz','r') as f:
    mol.atom = f.readlines()
mol.basis = '6-31gss'
mol.output = name+'-pyscf.log'
mol.verbose = 9
mol.max_memory = 20000
mol.build()

t1 = time.time()
ks = dft.RKS(mol)
ks.xc = 'b3lyp'
#ks.grids = gen_grid_sep.Grids(ks.mol)
ks.kernel()

t2 = time.time()
#dm = ks.make_rdm1()
os.system("echo '' > %s" % (name+'-xceda.log'))

atom_exc, dexc = xceda.analysis(ks)
#os.system("echo \"E_xc of C: %.10f\n\" >> %s" % (np.mean(atom_exc[0:38:4]), name+'-xceda.log'))
#os.system("echo \"E_xc of H: %.10f\n\" >> %s" % (np.mean(atom_exc[6:18]), name+'-xceda.log'))
t3 = time.time()
print(t2-t1, t3-t2)

