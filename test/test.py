from pyscf import gto,dft
import numpy as np
#import gen_grid_sep
#import numint_sep
import xceda

mol = gto.Mole()
mol.atom = '''C1	0.0000	0.0000	0.0000
H2	0.6282	0.6282	0.6282
H3	-0.6282	-0.6282	0.6282
H4	-0.6282	0.6282	-0.6282
H5	0.6282	-0.6282	-0.6282'''
mol.basis = '6-31gss'
mol.output = 'test.log'
mol.verbose = 9
mol.build()

ks = dft.RKS(mol)
ks.xc = 'b3lyp'
#ks.grids = gen_grid_sep.Grids(ks.mol)
ks.kernel()

dm = ks.make_rdm1()

xceda.get_atmexc(ks)

