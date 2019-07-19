from pyscf import gto,dft
import numpy as np
#import gen_grid_sep
#import numint_sep
import xceda
import os

name = 'C2H6_m062x'

#baslist = ['6-31+gss','6-31++gss','def2dzvp','def2tzvp','cc-pvdz','cc-pvtz','aug-cc-pvtz']
#xclist = [] 

mol = gto.Mole()
mol.atom = '''C1      0.0000  0.0000  0.7614
C2      0.0000  0.0000  -0.7614
H3      0.0000  1.0155  1.1574
H4      -0.8794 -0.5077 1.1574
H5      0.8794  -0.5077 1.1574
H6      0.0000  -1.0155 -1.1574
H7      -0.8794 0.5077  -1.1574
H8      0.8794  0.5077  -1.1574'''
mol.basis = '6-31gss'
mol.output = name+'-pyscf.log'
mol.verbose = 9
mol.build()

ks = dft.RKS(mol)
ks.xc = 'm062x'
#ks.grids = gen_grid_sep.Grids(ks.mol)
ks.kernel()

#dm = ks.make_rdm1()
os.system("echo '' > %s" % (name+'-xceda.log'))

atom_exc, dexc = xceda.analysis(ks)

#os.system("echo \"E_xc of C: %.10f\n\" >> %s" % (np.mean(atom_exc[0:2]), name+'-xceda.log'))
#os.system("echo \"E_xc of H: %.10f\n\" >> %s" % (np.mean(atom_exc[2:8]), name+'-xceda.log'))

