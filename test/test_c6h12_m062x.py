from pyscf import gto,dft
import numpy as np
#import gen_grid_sep
#import numint_sep
import xceda
import os
import time

name = 'C6H12_m062x'

#baslist = ['6-31+gss','6-31++gss','def2dzvp','def2tzvp','cc-pvdz','cc-pvtz','aug-cc-pvtz']
#xclist = [] 

mol = gto.Mole()
mol.atom = ''' C                 -2.86831901   -0.79141989    0.00000000
 C                 -1.35321301   -0.79141989    0.00000000
 C                 -0.80128201    0.61965811    0.00000000
 C                 -1.35094501    1.42419511    1.16066100
 C                 -2.86607001    1.42485611    1.16017200
 C                 -3.41887001    0.01423511    1.15887600
 H                  0.31731699    0.58566711    0.06271400
 H                 -0.98066401   -1.33705389    0.90656200
 H                 -0.97791901   -1.34127889   -0.90191000
 H                 -3.24100701   -0.35839089   -0.96538500
 H                 -3.24391601   -1.84563989    0.06350200
 H                 -0.97841001    0.98935011    2.12527200
 H                 -0.97492201    2.47833211    1.09866600
 H                 -3.24139001    1.97428411    2.06228600
 H                 -3.23792501    1.97145711    0.25384900
 H                 -3.15286601   -0.49148689    2.12419100
 H                 -4.53741801    0.04929611    1.09393800
 H                 -1.06990001    1.12529111   -0.96454600'''
mol.basis = '6-31gss'
mol.output = name+'-pyscf.log'
mol.verbose = 9
mol.max_memory = 20000
mol.build()

t1 = time.time()
ks = dft.RKS(mol)
ks.xc = 'm062x'
#ks.grids = gen_grid_sep.Grids(ks.mol)
ks.kernel()

t2 = time.time()
#dm = ks.make_rdm1()
os.system("echo '' > %s" % (name+'-xceda.log'))

atom_exc, dexc = xceda.analysis(ks)
os.system("echo \"E_xc of C: %.10f\n\" >> %s" % (np.mean(atom_exc[0:6]), name+'-xceda.log'))
os.system("echo \"E_xc of H: %.10f\n\" >> %s" % (np.mean(atom_exc[6:18]), name+'-xceda.log'))
t3 = time.time()
print(t2-t1, t3-t2)

