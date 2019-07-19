'''
Example for path/to/pyscf/pyscf/eda

Created Oct/15/2018
'''

from pyscf import gto, scf, dft, symm, qmmm
from pyscf.eda import edanew
from pyscf.scf import _vhf
from frame_small6 import preri
from pyscf.gto import moleintor
from h1e_old import h1e_old

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
#mol.cart=True
#mol.basis = 'cc-pvqz'
#mol.basis = '6-311g(d,p)'
#mol.basis = '6-31g(d,p)'
mol.basis = 'sto-3g'
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

#spinlist1 = [0,1,1]*30
#spinlist1 = [1,1]
#spinlist1 = [0,1,1,1,1]
spinlist1 = [0,1,1,1,0,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1]
mf = scf.RHF(mol)
mf.kernel()
tatm_= mol._atm
natm_ = []
for n in range(len(tatm_)):
    back = copy.copy(tatm_)
    for i in range(len(tatm_)):
        if (i!=n):
            back[i,0] = 0   
    natm_.append(back)
natm_ = np.array(natm_)
atm_ = mol._atm.T
atml = np.shape(mol._atm)[0]
bas_ = mol._bas.T 
basl = np.shape(mol._bas)[0]
env_ = mol._env
envl = np.shape(mol._env)[0]

print("mol._atm=",atm_.T,atml) 
print("mol._bas=",bas_.T,basl)
print("mol._env=",env_.T,envl)
#g = mf.nuc_grad_method()
#g.kernel()
pyscf_time = time.time()
print("pyscf_time=",pyscf_time-starttime)

dm = mf.make_rdm1()
nao = len(dm)
print("nao=",nao)
nbas = mol._bas.shape[0]
print("nbas=",nbas)
def include(a,b):
    if len(list(set(a) & set(b)))==0:
        return True
    else:
        return False

dml = []
E_nuc = []
E1 = []
basis = []
shls = []
singleatom = []

for i in realatomlabel:
    dm_new, e_nuc, e1, basis_range = edanew.get_dm(mf, dm, [i], atomlist1, spinlist1, 'hf')
    #singleatom.append(basis_range)
    basis.append(basis_range)
    basis_range = [i + 1 for i in basis_range]
    singleatom.append(basis_range)

bas2atm = np.zeros(nao)
for i in range(nao):
    for j in range(atom_number):
        if (i in basis[j]):
            bas2atm[i] = j
# -------------------------------------------------------- kinetic item 
atom_kinE = []
kin = mol.intor_symmetric('int1e_kin')
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
# -------------------------------------------------------- kinetic item 
# -------------------------------------------------------- Za*Zb/Rab item
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
# -------------------------------------------------------- nuc with elec item
Hcore = []
for i in range(len(natm_)):
    hcore = moleintor.getints("int1e_nuc_sph",natm_[i],mol._bas, mol._env)
    Hcore.append(hcore)
e1 = 0.0
atom_h1E = np.zeros(atom_number)
#for i in range(nao):
#    for j in range(nao):
#        for l in range(atom_number):
#            tem = dm[i,j] * Hcore[l][j,i]
#            a = int(bas2atm[i])
#            b = int(bas2atm[j])
#            buf = []
#            buf.append(a)
#            buf.append(b)
#            buf.append(l)
#            m = list(set(list(buf)))
#            body = len(m)
#            if body==1: 
#               atom_h1E[m[0]] = atom_h1E[m[0]] + tem
#            if body==2:
#               atom_h1E[m[0]] = atom_h1E[m[0]] + tem*1/2
#               atom_h1E[m[1]] = atom_h1E[m[1]] + tem*1/2
#            if body==3: 
#               atom_h1E[m[0]] = atom_h1E[m[0]] + tem*1/3
#               atom_h1E[m[1]] = atom_h1E[m[1]] + tem*1/3
#               atom_h1E[m[2]] = atom_h1E[m[2]] + tem*1/3
atom_h1E = h1e_old(dm,bas2atm,Hcore,atom_number,nao)               
atom_h1E = atom_h1E[0:atom_number]
print("atom_h1E=",atom_h1E,np.sum(atom_h1E))
# -------------------------------------------------------- nuc with elec item
e_coul = []
vhf = []
num = []; num2 = []; num3= []; num4=[]
for i in range(len(singleatom)):
    num.append(len(singleatom[i]))
num1 = (max(num))

for i in range(len(singleatom)):
    if len(singleatom[i])<num1:
        singleatom[i] = singleatom[i] + [0]*(num1-len(singleatom[i]))
singleitem = len(singleatom)
# -----------------------------------------------------------------------------------------  two elec item 
atom_energy = []
t1 = time.time()
atom_energy = preri(atm_,atml,bas_,basl,env_,envl,nao,nbas,dm,singleatom,singleitem,num1)
t2 = time.time()
print("T2-T1=",t2-t1)
atom_energy = np.array(atom_energy)[0:atom_number]
print("atom_eriE=",atom_energy)
#------------------------------------------------------------------------------------------  two elec item
Atom_energy = atom_kinE + atom_h1E + atom_nucE + atom_energy
print("Atom_energy=",Atom_energy,np.sum(Atom_energy))
endtime = time.time()
print("Run total time",endtime-starttime)
