from pyscf import scf
import numpy as np
from kit import logger, misc
import eda_inter

def get_Ejk(eda, jk='jk', jktype='py'):
    '''
    jk: j -> only calc E_j
        jk -> calc E_j, E_k
    jktype: def partition coeff for (ij|kl)
            atom-eq : find atom label a,b,c,d of i,j,k,l, and suppose atoms are equal
                      e.g. (ij|kl) -> from atom 1,2,3 -> coeff {1/3, 1/3, 1/3}
            bas-eq: each basis in (ij|kl) takes 1/4
    '''
    #e_coul = []
    mol = eda.mol
    dm = eda.dm
    atm2bas_p = eda.atm2bas_p
    atm2bas_f = eda.atm2bas_f

    if jktype=='py':
        vj,vk = scf.hf.get_jk(mol,dm)
        atom_ej = []
        atom_ek = []
        for i in range(mol.natm):
            ej = np.einsum('ij,ji',dm[atm2bas_p[i]],vj[:,atm2bas_p[i]])*.5
            atom_ej.append(ej)
            if jk=='jk':
                ek = np.einsum('ij,ji',dm[atm2bas_p[i]],-0.5*vk[:,atm2bas_p[i]])*.5
                atom_ek.append(ek)
        #with open(eda.output+'-eda.log','a') as f:
        logger.log(eda.stdout,"Atom_ej=",atom_ej)
        if jk=='jk':
            logger.log(eda.stdout,"Atom_ek=",atom_ek)
        if jk=='jk':
            atm_ejk = np.array(atom_ej), np.array(atom_ek)
        elif jk=='j':
            atm_ejk = np.array(atom_ej)
    elif jktype=='fort': 
        #atm2bas_f = misc.p2f(atm2bas_p)
        num = []
        for i in range(len(atm2bas_f)):
            num.append(len(atm2bas_f[i]))
        num1 = max(num)

        for i in range(len(atm2bas_f)):
            if len(atm2bas_f[i])<num1:
                atm2bas_f[i] = atm2bas_f[i] + [0]*(num1-len(atm2bas_f[i]))
        #singleitem = len(atm2bas_f)

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

        atom_ej, atom_ek = preri(atm_,atml,bas_,basl,env_,envl,nao,nbas,dm,atm2bas_f,singleitem,num1)
        #atom_energy = np.array(atom_energy)
        ##print("Atom_energy=",atom_energy)
        #with open(eda.output+'-eda.log','a') as f:
        if jk=='jk':
            atom_ej = np.array(atom_ej)[0:mol.natm]
            atom_ek = np.array(atom_ek)[0:mol.natm]
            logger.log(eda.stdout,"Atom_ej=",atom_ej)
            logger.log(eda.stdout,"Atom_ek=",atom_ek)
            if eda.anal:
                tot_aej = atom_ej.sum()
                tot_aek = atom_ek.sum()
                ej_err = tot_aej - 0.5*np.einsum('ij,ji',dm,eda.mf.get_j(mol,dm))
                ek_err = tot_aek - 0.5*(-0.5)*np.einsum('ij,ji',dm,eda.mf.get_k(mol,dm))
                logger.mlog(eda.stdout, "err_ej,ek: ",(ej_err,ek_err))
            return atom_ej, atom_ek
        elif jk=='j':
            atom_ej = np.array(atom_ej)[0:mol.natm]
            logger.log(eda.stdout,"Atom_ej=",atom_ej)
            if eda.anal:
                tot_aej = atom_ej.sum()
                ej_err = tot_aej - np.einsum('ij,ji',dm,eda.mf.get_j(mol,dm))
                logger.mlog(eda.stdout, "err_ej: ",(ej_err))
                return atom_ej

    if eda.showinter:
        #ejk1 = np.zeros(mol.natm)
        #ejk2 = np.zeros((mol.natm, mol.natm))
       
        #ejk3 = np.zeros((mol.natm, mol.natm, mol.natm))
        #ejk4 = np.zeros((mol.natm, mol.natm, mol.natm, mol.natm))
        ejk1,ejk2,ejk3,ejk4 = eda_inter.jk_inter(eda, atm2bas_p, 'jk') 
        atm_ejk = atm_ejk +( ejk1, ejk2 , ejk3, ejk4)
    return atm_ejk
