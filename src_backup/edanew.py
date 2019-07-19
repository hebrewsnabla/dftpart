'''
DFT-EDA for octane

Created May/11/2018
'''
from pyscf import gto, scf, lib, dft, symm, qmmm
import numpy as np
import copy
import time

#print(np.shape(dm))

#N = len(atomlist) # total numbers of fragments

def build_monofrag(mol, fragnumber, atomlist, spinlist):
    i = fragnumber
    #print(atomlist)
    startatom = sum(atomlist[0:i])
    endatom = startatom + atomlist[i]
    startbasis = mol.aoslice_by_atom()[startatom, 2]
    endbasis = mol.aoslice_by_atom()[endatom - 1, 3]
    basis_rangei = list(range(startbasis,endbasis))
    if isinstance(mol.atom, str):
        geomlines = mol.atom.split('\n')
        geomi = ''
        for line in geomlines[startatom:endatom]:
            geomi += (line + '\n')
    elif isinstance(mol.atom, list):
        geomlines = mol.atom
        geomi = []
        for line in geomlines[startatom:endatom]:
            geomi.append(list(line))
    else:
        print("wrong form of geom")
    spini = spinlist[i]
    return [geomi, spini, basis_rangei]

def build_frag(mol, fragnumbers = [], atomlist = [], spinlist = []):
    frag = gto.Mole()
    #frag.verbose = 4
    frag.basis = mol.basis
    basis_range = []
    #geom_frag = ''
    spin_frag = 0
    if isinstance(build_monofrag(mol, 0, atomlist, spinlist)[0], str):
        geom_frag = ''
    else:
        geom_frag = []
    for i in fragnumbers:
        geom_frag += build_monofrag(mol, i, atomlist, spinlist)[0]       
        spin_frag += build_monofrag(mol, i, atomlist, spinlist)[1]
        basis_range += build_monofrag(mol, i, atomlist, spinlist)[2]
    frag.atom = geom_frag
    frag.spin = spin_frag
    #frag.build()
    if len(fragnumbers)!=4:
        frag.build()
        h1e = scf.hf.get_hcore(frag)
    else:
        h1e = 0.0
    
    E_nuc = frag.energy_nuc()
    #endtime2 = datetime.datetime.now()
    return [h1e, E_nuc, basis_range]
#print(energy_nuc)

def get_dm(mf, dm, fragnumbers = [], atomlist = [], spinlist = [], method = 'hf'):
    mol = mf.mol
    h1e,E_nuc,basis_range = build_frag(mol, fragnumbers, atomlist, spinlist)
    #print('basis_range=',basis_range)
    #frag.cart = mol.cart
    #basis_range = build_frag(mol, fragnumbers, atomlist, spinlist)[2]
   #
   # dm_tot = np.array(mf.make_rdm1())
    dm_tot = dm
    #dm_new = np.zeros((np.shape(dm_tot)[0],np.shape(dm_tot)[1]))
    #dm_new[np.ix_(basis_range,basis_range)] = dm_tot[np.ix_(basis_range,basis_range)]
    dm_new = dm_tot[np.ix_(basis_range,basis_range)] 
    #print(dm_tot[np.ix_(basis_range,basis_range)])
    #if len(fragnumbers)!=4:
    #    e1 = np.einsum('ij,ji',h1e,dm_tot[np.ix_(basis_range,basis_range)]).real
    #else:
    #    e1 = 0.0
    if len(fragnumbers)!=4:
        e1 = np.einsum('ij,ji',h1e,dm_new).real
    else:
        e1 = 0.0
    return dm_new, E_nuc, e1, basis_range
    #delta = 0

    #E_elec = fragf.energy_elec(dm = dm)[0]
    #E_nuc = build_frag(mol, fragnumbers, atomlist, spinlist)[1]
    #E_tot = E_elec + E_nuc
    #return [E_elec, E_tot, delta]
    

def inter_energy(mf, fragnumbers = [], atomlist = [], spinlist = [], method = 'dft'):
    #xc = mf.xc
    '''
    oe = len(fragnumbers)%2
    for j in range(i):
        if set(command[j]).issubset(set(command[i])):
            #print(i,j)
            #print(set(command[j]),set(command[i]))
            #print(result_list[j])
            if len(command[j])%2 == oe:
                result += result_list[j]
            else:
                result -= result_list[j]
    result_output.append(result)
    '''
    N = len(fragnumbers)
    E_inter, E_h1e = get_energy(mf, fragnumbers, atomlist, spinlist, method)[1:3]
    if N == 2:
        for i in fragnumbers:
            a,b = get_energy(mf, [i], atomlist, spinlist, method)[1:3]
            E_inter = E_inter - a 
            E_h1e = E_h1e - b
    elif N == 3:
        for i in fragnumbers:
            for j in fragnumbers:
                if j > i:
                    a, b = get_energy(mf, [i,j], atomlist, spinlist, method)[1:3]
                    E_inter = E_inter - a
                    E_h1e = E_h1e - b
            a , b = get_energy(mf, [i], atomlist, spinlist, method)[1:3]
            E_inter = E_inter + a
            E_h1e = E_h1e + b
    elif N == 4:
        for i in fragnumbers:
            for j in fragnumbers:
                for k in fragnumbers:
                    if j > i and k > j:
                        a,b = get_energy(mf, [i,j,k], atomlist, spinlist, method)[1:3]
                        E_inter = E_inter - a
                        E_h1e = E_h1e - b    
                if j > i:
                    a, b= get_energy(mf, [i,j], atomlist, spinlist, method)[1:3]
                    E_inter = E_inter + a
                    E_h1e = E_h1e + b    
            a,b = get_energy(mf, [i], atomlist, spinlist, method)[1:3]
            E_inter = E_inter - a
            E_h1e = E_h1e - b
    elif N != 1:
        print("Invalid number of frags")
    
    return E_inter,E_h1e

        
def eda_full(mf, atomlist = [], spinlist = [], method = 'dft',coords = [], charges = []):
    result_list = []
    energy_list = []
    h1e_list = []
    command = []
    N = len(atomlist)
    for i in range(N):
        command.append([i])
    for i in range(N):
        for j in range(i+1,N):
            command.append([i,j])
    for i in range(N):
        for j in range(i+1,N):
            for k in range(j+1,N):
                command.append([i,j,k])
    for i in range(N):
        for j in range(i+1,N):
            for k in range(j+1,N):
                for l in range(k+1,N):
                    command.append([i,j,k,l])
    #print(command)
    
    for item in command:
        a,b = inter_energy(mf, item, atomlist, spinlist, method, coords, charges)
        energy_list.append(get_energy(mf, item, atomlist, spinlist, method,coords,charges)[1])
        result_list.append(a)
        h1e_list.append(b)
    print('result_list=',result_list)
    fulllog = mf.mol.output + "_eda_full"
    g = open(fulllog,'w')
    command_str_list_ = []
    command_str_list = []
    for item in command:
        s = ""
        s_ = ""
        for i in item:
            s += "%s," % str(i+1)
            s_ += "%s" % str(i+1)
        s = s.rstrip(',')
        command_str_list_.append(s_)
        command_str_list.append(s)
    for i in range(len(command)):
        g.write("E_%s = % .10g \n" % (command_str_list[i], result_list[i]))
    for i in range(len(command)):
        g.write("e_%s = % .10g \n" % (command_str_list[i], h1e_list[i]))
    for i in range(len(command)):
        g.write("E_%s = % .10g \n" % (command_str_list_[i], energy_list[i]))
    g.close() 



