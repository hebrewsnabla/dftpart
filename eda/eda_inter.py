#from scfeda import p2f
from .new_eda import preri
from ..kit import logger, misc
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


def jk_inter(eda, atm2bas_p, jk='jk2'):

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

    jk_sts = ['j', 'jk1', 'jk2']
    jkst = jk_sts.index(jk)
    atom_ej, atom_ek, e1, e2, e3, e4, ek1, ek2, ek3, ek4 = preri(atm_,atml,bas_,basl,env_,envl, eda.cart, nao,nbas, dm, 
        eda.bas2atm_f, eda.bas2frg, mol.natm, eda.nfrag, jkst)
    #atom_ejk = preri(atm_,atml,bas_,basl,env_,envl,nao,nbas,dm,atm2bas_f,singleitem,num1)
    #atom_energy = np.array(atom_energy)
    #print(atom_ejk)
    if jk == 'jk1':
        atom_ejk = atom_ej 
        ejksum = atom_ejk.sum()
        atom_ejk = np.array(atom_ejk)[0:eda.mol.natm]
        logger.log(eda.stdout_inter,"Atom_ejk=",atom_ejk)
    else:
        ejsum = atom_ej.sum()
        atom_ej = np.array(atom_ej)[0:eda.mol.natm]
        logger.log(eda.stdout_inter,"Atom_ej=",atom_ej)
        if jk == 'jk2':
            eksum = atom_ek.sum()
            atom_ek = np.array(atom_ek)[0:eda.mol.natm]
            logger.log(eda.stdout_inter,"Atom_ek=",atom_ek)

    ejk1 = np.array(e1)
    ejk2 = np.array(e2)
    ejk3 = np.array(e3)
    ejk3 = simp3(ejk3, eda.nfrag, eda.frag2layer)
    ejk4 = np.array(e4)
    ejk4 = simp4(ejk4, eda.nfrag, eda.frag2layer)
    #ejk1 = ej1+ek1
    #ejk2 = np.triu(ej2+ek2)
    if jk == 'jk2':
        ek1 = np.array(ek1)
        ek2 = np.array(ek2)
        ek3 = np.array(ek3)
        ek3 = simp3(ek3, eda.nfrag, eda.frag2layer)
        ek4 = np.array(ek4)
        ek4 = simp4(ek4, eda.nfrag, eda.frag2layer)
    if jk == 'jk1':
        logger.log(eda.stdout_inter,"ejk1=",ejk1)
        logger.log(eda.stdout_inter,"ejk2=",ejk2)
        logger.ilog(eda.stdout_inter,"ejk3=",ejk3)
        logger.ilog(eda.stdout_inter,"ejk4=",ejk4)
    else:
        logger.log(eda.stdout_inter,"ej1=",ejk1)
        logger.log(eda.stdout_inter,"ej2=",ejk2)
        logger.ilog(eda.stdout_inter,"ej3=",ejk3)
        logger.ilog(eda.stdout_inter,"ej4=",ejk4)
        if jk == 'jk2':
            logger.log(eda.stdout_inter,"ek1=",ek1)
            logger.log(eda.stdout_inter,"ek2=",ek2)
            logger.ilog(eda.stdout_inter,"ek3=",ek3)
            logger.ilog(eda.stdout_inter,"ek4=",ek4)
    
    interejksum = ejk1.sum() + ejk2.sum() + sum(ejk3.energies()) + sum(ejk4.energies())
    intereksum = ek1.sum() + ek2.sum() + sum(ek3.energies()) + sum(ek4.energies())
    if jk == 'jk1':
        logger.mlog(eda.stdout_inter,"interejksum=",interejksum)
        logger.mlog(eda.stdout_inter,"ejksum=",ejksum)
    else:
        logger.mlog(eda.stdout_inter,"interejsum=",interejksum)
        logger.mlog(eda.stdout_inter,"ejsum=",ejsum)
        if jk == 'jk2':
            logger.mlog(eda.stdout_inter,"intereksum=",intereksum)
            logger.mlog(eda.stdout_inter,"eksum=",eksum)


    #if eda.anal:
    #    tot_aej = atom_ej.sum()
    #    ej_err = tot_aej - np.einsum('ij,ji',dm,eda.mf.get_j(mol,dm))
    #    logger.mlog(eda.stdout, "err_ej: ",(ej_err))
    return ejk1,ejk2 , ejk3, ejk4, ek1, ek2, ek3, ek4

def get_atm2sub(natm, atomlist):
    atm2sub = {}
    for frag in atomlist:
        frag_num = atomlist.index(frag) +1
        for i in frag:
            atm2sub[i] = frag_num
    return atm2sub

def simp3(e3, nfrag, f2layer):
    e3simp = misc.EDict()
    for f in range(nfrag+2):
        for g in range(f+1,nfrag+2):
            for h in range(g+1, nfrag+2):
                fgh = (f+1,g+1,h+1)
                if abs(e3[f,g,h]) > 1e-12:
                    layers = (f2layer[f+1], f2layer[g+1], f2layer[h+1])
                    if 'cap' in layers:
                        continue
                    e3simp[fgh] = [e3[f,g,h], layers]
    return e3simp

def simp4(e4, nfrag, f2layer):
    e4simp = misc.EDict()
    for f in range(nfrag+2):
        for g in range(f+1,nfrag+2):
            for h in range(g+1, nfrag+2):
                for i in range(h+1, nfrag+2):
                    fghi = (f+1,g+1,h+1,i+1)
                    if abs(e4[f,g,h,i]) > 1e-12:
                        layers = (f2layer[f+1], f2layer[g+1], f2layer[h+1], f2layer[i+1])
                        if 'cap' in layers:
                            continue
                        e4simp[fghi] = [e4[f,g,h,i], layers]
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
    
