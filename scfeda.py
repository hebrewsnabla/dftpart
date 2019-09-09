'''

'''

from pyscf import gto, scf, dft, symm, qmmm
#import edanew
from pyscf.scf import _vhf
#from frame_small2 import preri
from new_eda import preri
from h1e_new import h1e, bgh1e, h1e_inter
#from h1e import h1e
import xceda
#from frame_small5 import preri
import bg
from pyscf.gto import moleintor
from pyscf.gto.mole import inter_distance
import numpy as np
import os,sys,subprocess
import time
import copy
#from gfea import gfea2_2
from QCKit import logger, dft_kit, gjf_kit

EDATHRESH = 1e-8

class EDA():
    r'''
    Attributes
        gjf : Gaussian gjf file
        method : e.g. ['m062x','6-31gss', ...]
                 '...' can be: charge -- GFEA's qmmm
                             qmmm -- common qmmm
                             cart
                             force
        output : A string like 'test'. 2 output file will be generated:
                 test-pyscf.log : PySCF output info
                 test-eda.log   : eda info

        built : whether build() is performed
        anal : check each part of SCF energy (J, K, nuc, etc.) equals the sum of their decomposition
    '''
    def __init__(self):
        self.mol = None

        self.mf = None
        self.dm = None
        #eda.nao = None
        #eda.nbas = None
        #self.atomlist = None
        #self.realatomlabel = None
        # must specify
        self.gjf = None
        self.method = None
        #self.chgstyle = 'gebf'
        #self.spinlist = None
        self.output = None
        self.verbose = 4
        self.built = False
        self.anal = False
        self.showinter = False
        ### GFEA ####
        self.molchgs = None # mol charges
        #self.envchgs = None # env frag charges
        self.bgchgs = None # background charges
        self.cen = None
        self.env = None # That's NOT mol.env_ !
    def build(self):
        return build(self, self.gjf, self.method)

    def kernel(self):
        if self.built == False:
            self.build()
        os.system("echo '' > "+self.output+'-eda.log')
        self.stdout = open(self.output+'-eda.log','a')
        if self.showinter:
            self.stdout_inter = open(self.output+'-inter.log','a')
        #with open(self.output+'-eda.log','a') as f:
        logger.mlog(self.stdout,"method,basis: ", self.method)
        t1 = time.time()
        atm2bas_f, atm2bas_p = get_atm2bas(self.mol)
        if self.showinter:
            atm_e1, e1_1, e1_2, e1_3 = get_E1(self,atm2bas_f)
        else:
            atm_e1 = get_E1(self,atm2bas_f)
        atm_enuc, enuc1, enuc2 = get_Enuc(self)
        t2 = time.time()
        #with open(self.output+'-eda.log','a') as f:
        logger.slog(self.stdout,"time for E1, E_nuc: %.5f\n", (t2-t1))
        if self.method[0] == 'hf':
            if self.showinter:
                atm_ejk, ejk1, ejk2, ejk3, ejk4 = get_Ejk(self, atm2bas_p,'jk')
                RR_inter = eda_inter.get_RR_inter(e1_1+enuc1+ejk1,e1_2+enuc2+ejk2,e1_3+ejk3,ejk4)
            else:
                atm_ej, atm_ek = get_Ejk(self, atm2bas_p,'jk')
                atm_ejk = atm_ej + atm_ek
            atm_E = atm_e1 + atm_enuc + atm_ejk
            t3 = time.time()
            #with open(self.output+'-eda.log','a') as f:
            logger.slog(self.stdout,"time for Ej, Ek: %.5f\n", (t3-t2))
        elif dft_kit.is_dft(self.method[0]):
            if self.showinter:
                atm_exc, atm_ej, ejxc1, ejxc2, ejxc3, ejxc4 = xceda.get_atmexc(self,atm2bas_p)
                RR_inter = eda_inter.get_RR_inter(e1_1+enuc1+ejxc1,e1_2+enuc2+ejxc2,e1_3+ejxc3,ejxc4)
            else:
                atm_exc, atm_ej = xceda.get_atmexc(self,atm2bas_p)
            atm_E = atm_e1 + atm_enuc + atm_ej + atm_exc
        if 'charge' in self.method:
            if self.showinter:
                bg_corrxn, bg_corrxn_fake, bg1,bg2,bg3 = get_bg_corrxn(self, atm2bas_f, 'charge')
                RC_inter = eda_inter.get_RC_inter(bg1,bg2,bg3)
            else:
                bg_corrxn, bg_corrxn_fake = get_bg_corrxn(self, atm2bas_f, 'charge')
            atm_E += bg_corrxn
        elif 'qmmm' in self.method:
            # show inter TBD
            bg_corrxn, bg_corrxn_fake = get_bg_corrxn(self, atm2bas_f, 'qmmm')
            atm_E += bg_corrxn
        #atm_ehf = atm_e1 + atm_enuc + atm_ej
        #atm_E = anal(self, atm_ek, atm_exc)
        totE = np.sum(atm_E)
        #outfile = xyznam[0:-3] + "atom"
        if '(' in self.output:
            op = self.output.replace('(','\(')
            self.output = op.replace(')','\)')
        scfE = subprocess.getstatusoutput("grep 'converged SCF energy' "+self.output+'-pyscf.log')[1]
        scfE = float(scfE.strip().split()[-1])
        #time.sleep(20)
        #with open(self.output+'-eda.log','a') as f:
        logger.slog(self.stdout,"tot Energy =%16.10f", (totE))
        logger.slog(self.stdout,"SCF Energy =%16.10f", (scfE))
        if ('charge' not in self.method) and ('qmmm' not in self.method):
            logger.slog(self.stdout,"Err of totE =%16.10f", (totE - scfE))
            conv = (totE - scfE) < EDATHRESH
        else:
            totE_fake = (atm_E - bg_corrxn + bg_corrxn_fake).sum()
            logger.slog(self.stdout,"fake tot Energy =%16.10f", (totE_fake))
            logger.slog(self.stdout,"Err of fake totE =%16.10f", (totE_fake - scfE))
            conv = (totE_fake - scfE) < EDATHRESH
        ok = ['FAIL','OK']
        logger.slog(self.stdout,ok[conv])

        for i in range(atm_E.shape[0]):
            logger.slog(self.stdout,"%s %i %16.10f", self.mol.atom_symbol(i),i+1,atm_E[i])
        inter_terms = {}
        if self.showinter: 
            inter_terms = RR_inter + RC_inter
            return atm_E, totE, conv, inter_terms
        else:
            return atm_E, totE, conv

def build(eda, gjf, method):

    starttime = time.time()
    #xyznam = sys.argv[1]
    mol = gto.Mole()
    ##with open(xyznam) as f:
    #   geom = f.read()
    mol.atom, coords, charges, molcharge, spin = gjf_kit.gjf_parser(gjf)
    if 'charge' in method:
        eda.bgchgs = (coords, charges)
        if eda.molchgs is None:
            logger.slog("Warning: center frag has no atomic charges")
    if 'qmmm' in method:
        eda.bgchgs = (coords, charges)
    #mol.cart=True
    mol.basis = method[1]
    #mol.symmetry = 1
    mol.output = eda.output + '-pyscf.log'
    mol.verbose = eda.verbose
    mol.build()
    eda.mol = mol

    if method[0] == 'hf':
        mf = scf.RHF(mol)
    elif dft_kit.is_dft(method[0]):
        mf = dft.RKS(mol)
        mf.xc = method[0]
        if 'ultrafine' in method:
            mf.grids.atom_grid = (99, 590)
    if ('charge' in method) or ('qmmm' in method):
        mf = qmmm.mm_charge(mf, coords, charges, unit='au')
    mf.kernel()
    eda.mf = mf
    if 'force' in method:
        g = mf.Gradients()
        eda.force = g.grad()
    pyscf_time = time.time()
    #print("pyscf_time=",pyscf_time-starttime)

    eda.dm = mf.make_rdm1()
    eda.built = True
    return eda

def get_atm2bas(mol):
    #dml = []
    #E_nuc = []
    #E1 = []
    #slice = []
    #basis = []
    #shls = []
    atm2bas_p = []
    atm2bas_f = []
    #twoatom = []
    #threeatom = []
    #fouratom = []

    #mol = eda.mol
    #realatomlabel = eda.realatomlabel
    #atomlist = eda.atomlist
    #spinlist = eda.spinlist
    for i in range(mol.natm):
        #dm_new, e_nuc, e1, basis_range = edanew.get_dm(mf, dm, [i], atomlist, spinlist, 'hf')
        #startatom = sum(atomlist[0:i])
        #endatom = startatom + atomlist[i]
        startbasis = mol.aoslice_by_atom()[i, 2]
        endbasis = mol.aoslice_by_atom()[i, 3]
        basis_range = list(range(startbasis,endbasis))
        basis_range_fort = [item + 1 for item in basis_range]
        basis_range_py = [item for item in basis_range]
        atm2bas_f.append(basis_range_fort)
        atm2bas_p.append(basis_range_py)

    return atm2bas_f, atm2bas_p

def get_bas2atm(atm2bas,nao,natm):
    bas2atm = np.zeros(nao)
    for i in range(1,nao+1):
        for j in range(natm):
            if (i in atm2bas[j]):
                bas2atm[i-1] = j
    return bas2atm


def get_bg_corrxn(eda, atm2bas, charge='charge'):
    '''
    For cen atom A and bg chg B
    bg correction (elec) = <i|B|j> - 0.5*Q_A*Q_B/R_AB
    bg correction (nuc)  = 0.5*Z_A* Q_B/R_AB
    bg corr fake = <i|B|j> + Z_A* Q_B/R_AB
    bg corr single qmmm = 0.5*<i|B|j> + 0.5*Z_A* Q_B/R_AB
    '''
    mol = eda.mol
    dm = eda.dm
    nao = len(dm)
    molchgs = eda.molchgs
    bgcoords, bgchgs = eda.bgchgs
    #cen = eda.cen
    #env = eda.env
    bg_corrxn = np.zeros(mol.natm)
    if eda.showinter:
        bg1 = np.zeros(mol.natm)
        bg2 = np.zeros((mol.natm, mol.natm))
        bg3 = np.zeros((mol.natm, mol.natm, mol.natm))
    vchg = bg.inter_elecbg(mol, dm, bgcoords, bgchgs)
    bas2atm = get_bas2atm(atm2bas,nao,mol.natm)
    if charge=='charge':
        atom_elecbg = 2*np.asarray(bgh1e(dm,bas2atm,vchg,mol.natm,nao))[0:mol.natm]
        atom_elecbg -= 0.5*bg.inter_bgbg(mol.atom_coords(), molchgs, bgcoords, bgchgs)
    elif charge=='qmmm':
        atom_elecbg = np.asarray(bgh1e(dm,bas2atm,vchg,mol.natm,nao))[0:mol.natm]
    logger.log(eda.stdout,"atom_elecbg=",atom_elecbg)

    if charge=='charge':
        atom_nucbg = bg.inter_nucbg(mol, bgcoords, bgchgs)
    elif charge=='qmmm':
        atom_nucbg = 0.5*bg.inter_nucbg(mol, bgcoords, bgchgs)
    logger.log(eda.stdout,"atom_nucbg=",atom_nucbg)

    bg_corrxn = atom_elecbg + atom_nucbg
    logger.log(eda.stdout,"atom_bg_correction=",bg_corrxn)
    bg_corrxn_fake = 2*np.asarray(bgh1e(dm,bas2atm,vchg,mol.natm,nao))[0:mol.natm] + bg.inter_nucbg(mol, bgcoords, bgchgs)
    if eda.showinter:
        return bg_corrxn, bg_corrxn_fake, bg1, bg2, bg3
    else:
        return bg_corrxn, bg_corrxn_fake

def get_E1(eda, atm2bas):
    #atom_kinE = []
    mol = eda.mol
    dm = eda.dm
    nao = len(dm)
    bas2atm = get_bas2atm(atm2bas,nao,mol.natm)
    kinmat = mol.intor_symmetric('int1e_kin')
    atom_kin = np.zeros(mol.natm)
    if eda.showinter:
        kin1 = np.zeros(mol.natm)
        kin2 = np.zeros((mol.natm, mol.natm))
    for i in range(nao):
        for j in range(nao):
            #print(dm,kin)
            kin = dm[i,j] * kinmat[i,j]
            a = int(bas2atm[i])
            b = int(bas2atm[j])
            atom_kin[a] = atom_kin[a] + kin
            atom_kin[b] = atom_kin[b] + kin
            if eda.showinter:
                if a==b:
                    kin1[a] = kin
                else:
                    kin2[a,b] = kin
    atom_kin = atom_kin/2
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_kinE=",atom_kin)
    if eda.showinter:
        logger.mlog(eda.stdout_inter,"kin1",kin1)
        logger.mlog(eda.stdout_inter,"kin2",kin2)

    fakeatm = []
    for n in range(mol.natm):
        back = copy.copy(mol._atm)
        for i in range(mol.natm):
            if (i!=n):
                back[i,0] = 0
        fakeatm.append(back)
    fakeatm = np.array(fakeatm)

    int1enuc = []
    for i in range(mol.natm):
        int1en = moleintor.getints("int1e_nuc_sph",fakeatm[i],mol._bas, mol._env)
        int1enuc.append(int1en)
    #e1 = 0.0
    #atom_h1E = np.zeros(atom_number)
    if eda.showinter:
        atom_1enuc, e1_1, e1_2, e1_3 = h1e_inter(dm,bas2atm,int1enuc,mol.natm,nao)
        atom_1enuc = np.asarray(atom_1enuc)[0:mol.natm]
        e1_1 = np.asarray(e1_1) + kin1
        e1_2 = np.asarray(e1_2) + kin2
        e1_3 = np.asarray(e1_3)
    else:
        atom_1enuc = np.asarray(h1e(dm,bas2atm,int1enuc,mol.natm,nao))[0:mol.natm]
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_1enuc=",atom_1enuc)
    if eda.showinter:
        logger.mlog(eda.stdout_inter,"e1_1",e1_1)
        logger.mlog(eda.stdout_inter,"e1_2",e1_2)
        logger.mlog(eda.stdout_inter,"e1_3",e1_3)
    atm_e1 = atom_kin + atom_1enuc
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_e1=",atm_e1)
    if eda.anal:
        tot_akin = atom_kin.sum()
        kin_err = tot_akin - np.einsum('ij,ji',dm,mol.intor_symmetric('int1e_kin'))
        tot_a1enuc = atom_1enuc.sum()
        a1enuc_err = tot_a1enuc - np.einsum('ij,ji',dm,moleintor.getints("int1e_nuc_sph",mol._atm,mol._bas, mol._env))
        logger.mlog(eda.stdout,"kin_err: ",kin_err)
        logger.mlog(eda.stdout,"1enuc_err: ",a1enuc_err)
    if eda.showinter: atm_e1 = atm_e1, e1_1, e1_2, e1_3
    return atm_e1

def get_Enuc(eda):
    mol = eda.mol
    charges = mol.atom_charges()
    if True: #eda.showinter:
        enuc1 = np.zeros(mol.natm)
        enuc2 = np.zeros((mol.natm, mol.natm))
    #coords = mol.atom_coords()
    rr = inter_distance(mol)
    #print(rr)
    rr[np.diag_indices_from(rr)] = 1e200
    #if CHECK_GEOM and numpy.any(rr < 1e-5):
    #    for atm_idx in numpy.argwhere(rr<1e-5):
    #    logger.warn(mol, 'Atoms %s have the same coordinates', atm_idx)
    #    raise RuntimeError('Ill geometry')
    atm_enucnuc = np.einsum('i,ij,j->i', charges, 1./rr, charges) * .5
    logger.log(eda.stdout,"atom_enucnuc=",atm_enucnuc)
    if eda.anal:
        tot_enucnuc = atm_enucnuc.sum()
        enucnuc_err = tot_enucnuc - np.einsum('i,ij,j', charges, 1./rr, charges) * .5
        logger.mlog(eda.stdout, "err_enucnuc", enucnuc_err)
    atm_enucnuc = atm_enucnuc, enuc1, enuc2
    return atm_enucnuc

def p2f(atm2bas_p):
    atm2bas_f = []
    for item in atm2bas_p:
        item_f = []
        for j in item:
            item_f.append(j+1)
        atm2bas_f.append(item_f)
    return atm2bas_f

def get_Ejk(eda, atm2bas_p, jk='jk', jktype='bas-eq'):
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
    
    if jktype=='bas-eq':
        vj,vk = scf.hf.get_jk(mol,dm)
        atom_ej = []
        atom_ek  =[]
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
    elif jktype=='atom-eq':
        atm2bas_f = p2f(atm2bas_p)
        num = []
        for i in range(len(atm2bas_f)):
            num.append(len(atm2bas_f[i]))
        num1 = (max(num))

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

    if eda.showinter==True:
        ejk1 = np.zeros(mol.natm)
        ejk2 = np.zeros((mol.natm, mol.natm))
        ejk3 = np.zeros((mol.natm, mol.natm, mol.natm))
        ejk4 = np.zeros((mol.natm, mol.natm, mol.natm, mol.natm))
        #ejk1,ejk2,ejk3,ejk4 = eda_inter.jk_inter(eda, atm2bas_p, 'jk') 
        atm_ejk = atm_ejk + (ejk1, ejk2, ejk3, ejk4)
    return atm_ejk

