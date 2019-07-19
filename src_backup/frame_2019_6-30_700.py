'''
Example for path/to/pyscf/pyscf/eda

Created Oct/15/2018
'''

from pyscf import gto, scf, dft, symm, qmmm
import edanew
from pyscf.scf import _vhf
#from frame_small2 import preri
from jkeda2 import preri
#from frame_small5 import preri
from pyscf.gto import moleintor

import numpy as np
import sys
import time

starttime = time.time()
xyznam = sys.argv[1]
mol = gto.Mole()
with open(xyznam) as f:
   geom = f.read()
mol.atom = geom
#mol.cart=True
#mol.basis = 'cc-pvtz'
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

#spinlist1 = [0,1,1]*10
#spinlist1 = [1,1]
#spinlist1 = [0,1,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,1]
#spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1]
spinlist1 = [0,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1]
mf = scf.RHF(mol)
mf.kernel()
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
slice = []
basis = []
shls = []
singleatom = []
twoatom = []
threeatom = []
fouratom = []

for i in realatomlabel:
    dm_new, e_nuc, e1, basis_range = edanew.get_dm(mf, dm, [i], atomlist1, spinlist1, 'hf')
    dml.append(dm_new)
    frag = gto.Mole()
    a = basis_range[0] + 1
    b = basis_range[-1] + 1
    slice.append([a,b])
    #singleatom.append(basis_range)
    e_nuc = frag.energy_nuc([mol.atom_charges()[i]],atom_coord[i])
    E_nuc.append(e_nuc)
    E1.append(e1)
    basis_range = [i + 1 for i in basis_range]
    singleatom.append(basis_range)

for i in realatomlabel:
    for j in realatomlabel:
        if j > i:
            dm_new, e_nuc, e1, basis_range = edanew.get_dm(mf, dm, [i,j], atomlist1, spinlist1, 'hf')
            #twoatom.append(basis_range)
            dml.append(dm_new)
            frag = gto.Mole()
            e_nuc = frag.energy_nuc([mol.atom_charges()[i],mol.atom_charges()[j]],atom_coord[[i,j]])
            E_nuc.append(e_nuc)
            E1.append(e1)
            basis_range = [i + 1 for i in basis_range]
            twoatom.append(basis_range)

for i in realatomlabel:
    for j in realatomlabel:
        for k in realatomlabel:
            if j>i and k>j:
                dm_new, e_nuc, e1, basis_range = edanew.get_dm(mf, dm,  [i,j,k], atomlist1, spinlist1, 'hf')
                #threeatom.append(basis_range)
                dml.append(dm_new)
                frag = gto.Mole()
                slice.append([a,b])
                e_nuc = frag.energy_nuc([mol.atom_charges()[i],mol.atom_charges()[j],mol.atom_charges()[k]],atom_coord[[i,j,k]])
                E_nuc.append(e_nuc)
                E1.append(e1)
                basis_range = [i + 1 for i in basis_range]
                threeatom.append(basis_range)

for i in realatomlabel:
    for j in realatomlabel:
        for k in realatomlabel:
            for l in realatomlabel:
                if j>i and k>j and l>k:
                    dm_new, e_nuc, e1, basis_range = edanew.get_dm(mf, dm, [i,j,k,l], atomlist1, spinlist1, 'hf')
                    #fouratom.append(basis_range)
                    dml.append(dm_new)
                    frag = gto.Mole()
                    a = basis_range[0] + 1 
                    b = basis_range[-1] + 1 
                    slice.append([a,b])
                    e_nuc = frag.energy_nuc([mol.atom_charges()[i],mol.atom_charges()[j],mol.atom_charges()[k],mol.atom_charges()[l]],atom_coord[[i,j,k,l]])
                    E_nuc.append(e_nuc)
                    e1 = 0.0
                    E1.append(e1)
                    basis_range = [i + 1 for i in basis_range]
                    fouratom.append(basis_range)
b = time.time()
print("combination-time=",b-pyscf_time)
e_coul = []
vhf = []

#eri_s8 = mol.intor("int2e_sph",aosym='s8')
#eri = mol.intor("int2e_sph")


#print(max(singleatom))
num = []; num2 = []; num3= []; num4=[]
for i in range(len(singleatom)):
    num.append(len(singleatom[i]))
num1 = (max(num))
for i in range(len(twoatom)):
    num2.append(len(twoatom[i]))
num2 = (max(num2))
#for i in range(len(threeatom)):
#    num3.append(len(threeatom[i]))
#num3 = (max(num3))
#for i in range(len(fouratom)):
#    num4.append(len(fouratom[i]))
#num4 = (max(num4))
#print(num1,num2,num3,num4)

for i in range(len(singleatom)):
    if len(singleatom[i])<num1:
        singleatom[i] = singleatom[i] + [0]*(num1-len(singleatom[i]))
for i in range(len(twoatom)):
    if len(twoatom[i])<num2:
        twoatom[i] = twoatom[i] + [0]*(num2-len(twoatom[i]))
#for i in range(len(threeatom)):
#    if len(threeatom[i])<num3:
#        threeatom[i] = threeatom[i] + [0]*(num3-len(threeatom[i]))
#for i in range(len(fouratom)):
#    if len(fouratom[i])<num4:
#        fouratom[i] = fouratom[i] + [0]*(num4-len(fouratom[i]))
singleitem = len(singleatom)
twoitem = len(twoatom)
#threeitem = len(threeatom)
#fouritem = len(fouratom)
#print(singleitem,twoitem,threeitem,fouritem)

#--------------------------------calc e_coul step by step
#for i in range(nbas):
#    tem = moleintor.getints('int2e_cart',mol._atm, mol._bas, mol._env,(i,i+1,0,nbas,0,nbas,0,nbas))
#    a = tem[0].reshape((1,nao,nao,nao))
#    print(np.shape(a))
#print("singleatom=",singleatom)
#print("singleitem=",singleitem)
#print("num1=",num1)
#exit(1)
# ---------------------------- restructure slice 
#e_coul = ecoul(eri,nao,dm,list(range(1,nao+1)),nao)
# ---------------------------------------------------
atom_energy = []
slice = 1
t1 =time.time()
#tem = moleintor.getints('int2e_cart',mol._atm, mol._bas, mol._env,shls_slice = None,aosym='s8')
#print("TEM=",tem)
#exit(1)
#gto.getints_by_shell('int2e_sph', (0,1,0,1), mol._atm, mol._bas, mol._env, comp=3)


########################################################################
#tem = moleintor.getints('int2e_sph',mol._atm, mol._bas, mol._env,(0,nbas,0,nbas,0,nbas,0,nbas))
#print("sad",np.shape(tem))
#for i in range(nao):
#    print("tem=",tem[i][1][2][3])
t1 = time.time()
atom_energy = preri(atm_,atml,bas_,basl,env_,envl,nao,nbas,dm,singleatom,singleitem,num1)
t2 = time.time()
print("T2-T1=",t2-t1)
print("Run total time",t2-starttime)
#print("atom_energy=",atom_energy)
#atom_energy = preri(atm_,atml,bas_,basl,env_,envl)
#exit(1)
#########################################################################
#
## ---------------------------------------------------------------------- SINGLE RUN
#t1 = time.time()
#tem = moleintor.getints('int2e_sph',mol._atm, mol._bas, mol._env,(0,1,0,nbas,0,nbas,0,nbas))
#t2 = time.time()
#a = tem[0].reshape((1,nao,nao,nao))
#atom_energy.append(preri(a,slice,nao,dm,singleatom,singleitem,num1))
#t3 = time.time()
#print("Single cycle=",t2-t1,' ',t3-t2,' ',t3-t1)
#exit(1)
## ----------------------------------------------------------------------
#for i in range(nbas):
#    tem = moleintor.getints('int2e_sph',mol._atm, mol._bas, mol._env,(i,i+1,0,nbas,0,nbas,0,nbas))
#    #print("SHAPE1=",np.shape(tem))
#    atom = np.shape(tem)[0]
#    #print("SHAPE=",atom)
#    a = tem[0].reshape((1,nao,nao,nao))
#    if (atom==1):
#        # ---------------------------- call fortran subroutine
#        #atom_energy = preri(a,nao,dm,singleatom,singleitem,num1)
#         atom_energy.append(preri(a,slice,nao,dm,singleatom,singleitem,num1))
#         slice+=1
#         #print("atom")
#        #print("Atom_energy=",atom_energy)
#        # ---------------------------- end call fortran
#    else:
#       for j in range(atom):
#            a = tem[j].reshape((1,nao,nao,nao))
#           # ---------------------------- call fortran subroutine
#           #atom_energy = preri(a,nao,dm,singleatom,singleitem,num1)
#            atom_energy.append(preri(a,slice,nao,dm,singleatom,singleitem,num1))
#            slice+=1
#           #print("Atom_energy=",atom_energy)
#           # ---------------------------- end call fortran
#       #print(tem[]
#t2 = time.time()
#print("Time=",t2-t1)
#atom_energy = np.array(atom_energy)
##print("Atom_energy=",atom_energy)
print("atom_ej=",atom_energy[0])
print("atom_ek=",atom_energy[1])
atom_energy = np.array(atom_energy[0])[0:atom_number] + np.array(atom_energy[1])[0:atom_number]
print("Atom_energy=",atom_energy)
#atom_energy = np.array(atom_energy).sum()
#print("Atom_energy=",atom_energy)
#-------------------------------------------------------
#print(len(E_coul_list),len(E1),len(E_nuc))
#Elec = np.array(E_coul_list) + np.array(E1) + np.array(E_nuc)
Elec = np.array(E1) + np.array(E_nuc)

frag = list(range(atom_number)) 
command = []
for i in frag:
    command.append((i,))

for i in frag:
    for j in frag:
        if j>i:
            command.append((i,j))

for i in frag:
    for j in frag:
        for k in frag:
            if k>j and j>i:
                command.append((i,j,k))

for i in frag:
    for j in frag:
        for k in frag:
            for l in frag:
                if l>k and k>j and j>i:
                    command.append((i,j,k,l))
energy = {}
E_coul = {}
i = 0
for item in command:
    energy[item] = Elec[i]
    #E_coul[item] = E_coul_list[i]
    i += 1

def inter_energy(energy):
    blank = []
    for item in command:
        if len(item)==1:  # single-body energy
            blank.append(energy[item])
        elif len(item)==2:  # two-body interaction
            #print(item)
            for element in item:
                #print(energy[item],energy[(item[0],)],energy[(item[1],)])
                blank.append(energy[item] - energy[(item[0],)] - energy[(item[1],)])
                break
        elif len(item)==3:  # three-body interaction 
            for element in item:
                #print(item)
                EIJK = energy[item]
                EIJ = 0.0
                for i in item:
                    for j in item:
                        if j>i:
                            Eij = energy[(i,j)]
                            #print(Eij)
                            EIJ+= Eij
                        EI = 0.0
                for i in item:
                    Ei = energy[(i,)]
                    #print(Ei)
                    EI += Ei 
                #print(Eijk - EIJ + EI)
                blank.append(EIJK - EIJ + EI)
                break
        elif len(item)==4:  # four-body interaction 
            for element in item:
                EIJKL = energy[item]
                EIJK = 0.0 
                for i in item:
                    for j in item:
                        for k in item:
                            if k>j and j>i:
                                Eijk = energy[(i,j,k)]
                                EIJK+= Eijk
                EIJ = 0.0
                for i in item:
                    for j in item:
                        if j>i:
                            Eij = energy[(i,j)]
                            EIJ+= Eij
                EI = 0.0
                for i in item:
                    Ei = energy[(i,)]
                    EI+= Ei
                blank.append(EIJKL - EIJK + EIJ -EI)
                break
    return blank
# -------------------------------------------------------------------------------------
blank = inter_energy(energy)
#blank1 = inter_energy(E_coul)

Inter_energy = {}
i = 0
for item in command:
    if len(item)!=4:
        Inter_energy[item] = blank[i]
        i+=1
    else:
        Inter_energy[item] = 0.0
        i+=1
sum = 0.0
for i in range(len(command)):
    #print(i,Inter_energy[command[i]])
    sum+= Inter_energy[command[i]]
#print(sum)
atomenergy = []
for central_number in frag:
    central_number = int(central_number)
    central_energy = 0.0
    for i in range(len(blank)):
        if len(command[i])==1 and central_number in command[i]:
            central_energy+= Inter_energy[command[i]]
            #print(command[i])
        elif len(command[i])==2 and central_number in command[i]:
            central_energy+= 1/2 * Inter_energy[command[i]]
            #print(command[i])
        elif len(command[i])==3 and central_number in command[i]:
            central_energy+= 1/3 * Inter_energy[command[i]]
            #print(command[i])
        elif len(command[i])==4 and central_number in command[i]:
            central_energy+= 1/4 * Inter_energy[command[i]]
            #print(command[i])
    atomenergy.append(central_energy)

atomenergy = np.array(atomenergy) + atom_energy
scfenergy = np.sum(atomenergy)
outfile = xyznam[0:-3] + "atom"
with open(outfile,"w") as f:
    f.write("SCF Energy=%16.6f\n" %(scfenergy))
    for i in range(len(atomenergy)):
        f.write("%i%16.6f\n" %(i+1,atomenergy[i]))
endtime = time.time()
print('timeTot=',(endtime-starttime))
#print(command)
#print(len(command))
print('Hello wolrd')

