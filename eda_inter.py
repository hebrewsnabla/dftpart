#from scfeda import p2f
from new_eda import preri
from kit import logger
#import pymp
import numpy as np

def p2f(atm2bas_p):
    atm2bas_f = []
    for item in atm2bas_p:
        item_f = []
        for j in item:
            item_f.append(j+1)
        atm2bas_f.append(item_f)
    return atm2bas_f


def jk_inter(eda, atm2bas_p, jk='jk'):

    atm_ = eda.mol._atm.T
    atml = np.shape(eda.mol._atm)[0]
    bas_ = eda.mol._bas.T
    basl = np.shape(eda.mol._bas)[0]
    env_ = eda.mol._env
    envl = np.shape(eda.mol._env)[0]
    dm = eda.dm
    #logger.log(eda.stdout,"dm",dm[-5:])
    mol = eda.mol
    nao = len(dm)
    nbas = eda.mol._bas.shape[0]
    #if eda.verbose > 5:
    #logger.mlog(eda.stdout,"atm_",atm_)
    #logger.mlog(eda.stdout,"bas_",bas_)
    #logger.mlog(eda.stdout,"env_",env_)
    logger.slog(eda.stdout,"nao=%d",nao)    # num of cGTO
    logger.slog(eda.stdout,"nbas=%d",nbas)  # num of shells

    atom_ejk, e1, e2, e3, e4 = preri(atm_,atml,bas_,basl,env_,envl, eda.cart, nao,nbas, dm, 
        eda.bas2atm_f, eda.bas2frg, mol.natm, eda.totnum_frag)
    #atom_ejk = preri(atm_,atml,bas_,basl,env_,envl,nao,nbas,dm,atm2bas_f,singleitem,num1)
    #atom_energy = np.array(atom_energy)
    #print(atom_ejk)
    ejksum = atom_ejk.sum()
    atom_ejk = np.array(atom_ejk)[0:eda.mol.natm]
    logger.log(eda.stdout_inter,"Atom_ejk=",atom_ejk)
    ejk1 = np.array(e1)
    ejk2 = np.array(e2)
    #ek1 = np.array(ek1)
    #ek2 = np.array(ek2)
    ejk3 = np.array(e3)
    ejk3 = simp3(ejk3, eda.nfrag)
    ejk4 = np.array(e4)
    ejk4 = simp4(ejk4, eda.nfrag)
    #ejk1 = ej1+ek1
    #ejk2 = np.triu(ej2+ek2)
    logger.log(eda.stdout_inter,"ejk1=",ejk1)
    logger.log(eda.stdout_inter,"ejk2=",ejk2)
    logger.mlog(eda.stdout_inter,"ejk3=",ejk3)
    logger.mlog(eda.stdout_inter,"ejk4=",ejk4)
    
    interejksum = ejk1.sum() + ejk2.sum() + sum(ejk3.values()) + sum(ejk4.values())
    logger.mlog(eda.stdout_inter,"interejksum=",interejksum)
    logger.mlog(eda.stdout_inter,"ejksum=",ejksum)
    #if eda.anal:
    #    tot_aej = atom_ej.sum()
    #    ej_err = tot_aej - np.einsum('ij,ji',dm,eda.mf.get_j(mol,dm))
    #    logger.mlog(eda.stdout, "err_ej: ",(ej_err))
    return ejk1,ejk2 , ejk3, ejk4

def get_atm2sub(natm, atomlist):
    atm2sub = {}
    for frag in atomlist:
        frag_num = atomlist.index(frag) +1
        for i in frag:
            atm2sub[i] = frag_num
    return atm2sub

def simp3(e3, nfrag):
    e3simp = {}
    for f in range(nfrag+2):
        for g in range(f+1,nfrag+2):
            for h in range(g+1, nfrag+2):
                fgh = "%d,%d,%d"%(f,g,h)
                if abs(e3[f,g,h]) > 1e-12:
                    e3simp[fgh] = e3[f,g,h]
    return e3simp

def simp4(e4, nfrag):
    e4simp = {}
    for f in range(nfrag+2):
        for g in range(f+1,nfrag+2):
            for h in range(g+1, nfrag+2):
                for i in range(h+1, nfrag+2):
                    fghi = "%d,%d,%d,%d"%(f,g,h,i)
                    if abs(e4[f,g,h,i]) > 1e-12:
                        e4simp[fghi] = e4[f,g,h,i]
    return e4simp

def get_RR_inter(e1, e2):
    natm = e1.shape[0]
    e_1 = {}
    for i in range(natm):
        e_1[i+1] = e1[i]
    e_2 = {}
    for i in range(natm):
        for j in range(natm):
            ij = tuple(sorted((i+1,j+1)))
            if ij not in e_2:
                e_2[ij] = e2[i,j]
            else:
                e_2[ij] += e2[i,j]
"""
    e_3 = {}
    for i in range(natm):
        for j in range(natm):
            for k in range(natm):
                ijk = tuple([i+1,j+1,k+1].sorted())
                if ijk not in e_3:
                    e_3[ijk] = e3[i,j,k]
                else:
                    e_3[ijk] += e3[i,j,k]
    #e_4 = {}
    with pymp.Parallel(10) as p:
        e4 = pymp.shared.array(e4)
        e_4 = pymp.shared.dict()
        for i in p.range(natm):
            for j in range(natm):
                for k in range(natm):
                    for l in range(natm):
                        ijkl = tuple([i+1,j+1,k+1,l+1].sorted())
                    if ijkl not in e_3:
                        e_3[ijkl] = e3[i,j,k,l]
                    else:
                        e_3[ijkl] += e3[i,j,k,l]
"""                    
    
