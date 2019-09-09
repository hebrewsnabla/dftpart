from gfea3 import p2f
from new_eda import preri
from QCKit import logger
import pymp

def jk_inter(eda, atm2bas_p, jk='jk'):
    atm2bas_f = p2f(atm2bas_p)
    num = []
    for i in range(len(atm2bas_f)):
        num.append(len(atm2bas_f[i]))
    num1 = max(num)

    for i in range(len(atm2bas_f)):
        if len(atm2bas_f[i])<num1:
            atm2bas_f[i] = atm2bas_f[i] + [0]*(num1-len(atm2bas_f[i]))
    singleitem = len(atm2bas_f)

    atm_ = mol._atm.T
    atml = np.shape(mol._atm)[0]
    bas_ = mol._bas.T
    basl = np.shape(mol._bas)[0]
    env_ = mol._env
    envl = np.shape(mol._env)[0]
    nao = len(dm)
    nbas = mol._bas.shape[0]
    #if eda.verbose > 5:
    #logger.mlog(eda.stdout,"atm_",atm_)
    #logger.mlog(eda.stdout,"bas_",bas_)
    #logger.mlog(eda.stdout,"env_",env_)
    logger.slog(eda.stdout,"nao=%d",nao)
    logger.slog(eda.stdout,"nbas=%d",nbas)

    atom_ejk, ejk1,ejk2, ejk3, ejk4 = preri(atm_,atml,bas_,basl,env_,envl,nao,nbas,dm,atm2bas_f,singleitem,num1)
    #atom_ejk = preri(atm_,atml,bas_,basl,env_,envl,nao,nbas,dm,atm2bas_f,singleitem,num1)
    #atom_energy = np.array(atom_energy)
    #print(atom_ejk)
    atom_ejk = np.array(atom_ejk)[0:mol.natm]
    #logger.log(eda.stdout,"Atom_ejk=",atom_ejk)
    ejk1 = np.array(ejk1)
    ejk2 = np.array(ejk2)
    ejk3 = np.array(ejk3)
    ejk4 = np.array(ejk4)
    logger.log(eda.stdout_inter,"ejk1=",ejk1)
    logger.log(eda.stdout_inter,"ejk2=",ejk2)
    logger.log(eda.stdout_inter,"ejk3=",ejk3)
    logger.log(eda.stdout_inter,"ejk4=",ejk4)
    #if eda.anal:
    #    tot_aej = atom_ej.sum()
    #    ej_err = tot_aej - np.einsum('ij,ji',dm,eda.mf.get_j(mol,dm))
    #    logger.mlog(eda.stdout, "err_ej: ",(ej_err))
    return ejk1,ejk2, ejk3, ejk4

def get_atm2sub(natm, atomlist):
    atm2sub = {}
    for frag in atomlist:
        frag_num = atomlist.index(frag) +1
        for i in frag:
            atm2sub[i] = frag_num
    return atm2sub

def get_RR_inter(e1, e2, e3, e4):
    natm = e1.shape[0]
    e_1 = {}
    for i in range(natm):
        e_1[i+1] = e1[i]
    e_2 = {}
    for i in range(natm):
        for j in range(natm):
            ij = tuple([i+1,j+1].sorted())
            if ij not in e_2:
                e_2[ij] = e2[i,j]
            else:
                e_2[ij] += e2[i,j]
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
                    
    