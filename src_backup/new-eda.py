'''
Example for path/to/pyscf/pyscf/eda

Created Oct/15/2018
'''

from pyscf import gto, scf, dft, symm, qmmm,tools
from pyscf.eda import edanew
from pyscf.scf import _vhf
#from frame import ecoul
#from frame_old import preri
#from frame_run import preri
from pyscf.gto import moleintor
from h1e_new import h1e
#from h1e import h1e

import numpy as np
import sys
import time
import copy

starttime = time.time()
xyznam = sys.argv[1]
mol = gto.Mole()
with open(xyznam) as f:
   geom = f.read()
mol.atom = geom
mol.cart=False
#mol.basis = 'cc-pvqz'
#mol.basis = '6-311g(d,p)'
mol.basis = '6-31g(d,p)'
#mol.basis = 'sto-3g'
mol.verbose = 4
#mol.symmetry = 1
mol.output = xyznam[0:-4] + '-pyscf.log'
mol.build()
atom_number = len(mol._atom)

atom_coord = np.zeros((atom_number,3))
for i in range(atom_number):
    atom_coord[i] = mol.atom_coord(i)
atomlist1 = np.ones(atom_number,dtype=int)

realatomlabel = list(range(atom_number))

spinlist1 = [0,1,1]*12
#spinlist1 = [0,1,1,1,0,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1]
mf = scf.RHF(mol)
mf.kernel()

#g = mf.nuc_grad_method()
#g.kernel()

dm = mf.make_rdm1()
nao = len(dm)

atm_ = mol._atm
bas_ = mol._bas
env_ = mol._env
print("mol._atm=",atm_)
#print(nao)
vj,vk = scf.hf.get_jk(mol,dm)
endtime = time.time()
print('pyscf-time=',(endtime-starttime))
#print(vj,vk)
# -----------------------------------------------------------------------------------
a = time.time()
natm_ = []
for n in range(len(atm_)):
    back = copy.copy(atm_)
    for i in range(len(atm_)):
        if (i!=n):
            back[i,0] = 0
    natm_.append(back)
natm_ = np.array(natm_)
#print("natm_",natm_)

#f12 = open("1.mat","w")
Hcore = []
for i in range(len(natm_)):
  hcore = moleintor.getints("int1e_nuc_sph",natm_[i],mol._bas, mol._env)
  #if i==2:
  #  tools.dump_mat.dump_rec(f12,hcore,ncol=16,digits=4)
  Hcore.append(hcore)
#f12.close()
b = time.time()
print("Time_calc_Hcore=",b-a)
#tem = moleintor.getints("int1e_nuc_sph",mol._atm,mol._bas, mol._env)
#etem = np.einsum('ij,ji',tem,dm).real
#print("h1e_tot=",etem)

a = time.time()
atom_energy = []

kin = mol.intor_symmetric('int1e_kin')
basis = []
for i in realatomlabel:
    dm_new, e_nuc, e1, basis_range = edanew.get_dm(mf, dm, [i], atomlist1, spinlist1, 'hf')
    basis.append(basis_range)
b = time.time()
print("Time_basis=",b - a)
#for i in realatomlabel:
#    tem = []
#    for j in realatomlabel:
#        if (j!=i):
#            tem.append(basis[j])
#    outbasis.append(tem)

bas2atm = np.zeros(nao)
for i in range(nao):
    for j in range(atom_number):
        if (i in basis[j]):
            bas2atm[i] = j
# -----------------------------------------------------------------------------------
# -------------------------------------------------------- kinetic item 
c = time.time()
atom_kinE = []
#for i in range(atom_number):
#    kinE = np.sum(np.dot(dm[basis[i]],kin[:,basis[i]]))
#    atom_kinE.append(kinE)
atom_kinE = np.zeros(atom_number)
for i in range(nao):
    for j in range(nao):
        kinE = dm[i,j] * kin[i,j]
        a = int(bas2atm[i])
        b = int(bas2atm[j])
        atom_kinE[a] = atom_kinE[a] + kinE
        atom_kinE[b] = atom_kinE[b] + kinE
atom_kinE = atom_kinE/2
print("atom_kinE",atom_kinE,np.sum(atom_kinE))
d = time.time()
print("Time_kinE=",d-c)
# -------------------------------------------------------- kinetic item 
# -------------------------------------------------------- Za*Zb/Rab item
#print("mol.atom_charges(),atom_coord",mol.atom_charges(),atom_coord)
a = time.time()
totnuc = mf.energy_nuc()
atom_nucE = []
charges = mol.atom_charges(); coords = atom_coord
for i in range(atom_number):
    q1 = charges[i]
    r1 = coords[i]
    nucE = 0.0
    for j in range(atom_number):
        if (j!=i):
          q2 = charges[j]
          r2 = coords[j]
          r = np.linalg.norm(r2-r1,ord=2)
          nucE += q1 * q2/r
    atom_nucE.append(nucE)
atom_nucE = np.array(atom_nucE)/2
print("atom_nucE=",atom_nucE,np.sum(atom_nucE),totnuc)
b = time.time()
print("Time_nucE=",b-a)
# -------------------------------------------------------- Za*Zb/Rab item
# -------------------------------------------------------- nuc with elec item
a = time.time()
e1 = 0.0
atom_h1E = np.zeros(atom_number)
# -------------------------------------------------------- 
#for i in range(nao):
#    for j in range(nao):
#        for l in range(atom_number):
#            tem = dm[i,j] * Hcore[l][j,i] 
#            atom_h1E[l] = atom_h1E[l] + tem*1/3
#            atom_h1E[int(bas2atm[i])] = atom_h1E[int(bas2atm[i])] + tem*1/3
#            atom_h1E[int(bas2atm[j])] = atom_h1E[int(bas2atm[j])] + tem*1/3
atom_h1E = h1e(dm,bas2atm,Hcore,atom_number,nao)
print("atom_h1E=",atom_h1E)
# --------------------------------------------------------- use fortran to replace it
b = time.time()
print("Time_h1E=",b-a)
# -------------------------------------------------------- nuc with elec item
a = time.time()
for i in range(atom_number):
    e_coul = np.einsum('ij,ji',dm[basis[i]],(vj[:,basis[i]]-0.5*vk[:,basis[i]]))*.5
    atom_energy.append(e_coul)
b = time.time()
print("TIME_eri",b-a)
Atom_energy = atom_kinE + atom_h1E + atom_nucE + atom_energy 
print("Atom_energy=",Atom_energy,np.sum(Atom_energy))
endtime = time.time()
print('timeTot=',(endtime-starttime))
print('Hello wolrd')

