'''

'''
# PySCF module
from pyscf import gto, scf, dft, symm, qmmm
from pyscf.gto import moleintor
from pyscf.gto.mole import inter_distance
from pyscf.scf import _vhf

# dftpart module
from new_eda import preri
from h1e_new import h1e, bgh1e, h1e_inter
import xceda, jkeda
import bg
import eda_inter
from kit import logger, dft_kit, gjf_kit, misc

# others
import numpy as np
import os,sys,subprocess 
import time
import copy


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
            os.system("echo '' > "+self.output+'-inter.log')
            self.stdout_inter = open(self.output+'-inter.log','a')
        #with open(self.output+'-eda.log','a') as f:
        logger.mlog(self.stdout,"method,basis: ", self.method)
        t1 = time.time()
        atm2bas_f, atm2bas_p = get_atm2bas(self.mol)
        if self.showinter:
            atm_e1, e1_1, e1_2, e1_3 = get_E1(self,atm2bas_f)
        else:
            atm_e1 = get_E1(self,atm2bas_f)
        atm_enuc, enuc2 = get_Enuc(self)
        t2 = time.time()
        #with open(self.output+'-eda.log','a') as f:
        logger.slog(self.stdout,"time for E1, E_nuc: %.5f\n", (t2-t1))
        if self.method[0] == 'hf':
            if self.showinter:
                atm_ej, atm_ek, ejk1, ejk2, ejk3, ejk4 = jkeda.get_Ejk(self, atm2bas_p,'jk')
                atm_ejk = atm_ej + atm_ek
                RR1 = e1_1 + ejk1
                RR2 = e1_2 + enuc2 + ejk2
                RR3 = e1_3 + ejk3
                RR4 = ejk4
                #RR_inter = eda_inter.get_RR_inter(e1_1+ejk1,e1_2+enuc2+ejk2)
            else:
                atm_ej, atm_ek = jkeda.get_Ejk(self, atm2bas_p,'jk')
                atm_ejk = atm_ej + atm_ek
            atm_E = atm_e1 + atm_enuc + atm_ejk
            t3 = time.time()
            #with open(self.output+'-eda.log','a') as f:
            logger.slog(self.stdout,"time for Ej, Ek: %.5f\n", (t3-t2))
        elif dft_kit.is_dft(self.method[0]):
            if self.showinter:
                atm_exc, atm_ej, ejxc1, ejxc2, ejxc3, ejxc4 = xceda.get_atmexc(self,atm2bas_p)
                RR1 = e1_1 + ejxc1
                RR2 = e1_2 + enuc2 + ejxc2
                RR3 = e1_3 + ejxc3
                RR4 = ejxc4
                #RR_inter = eda_inter.get_RR_inter(e1_1+ejxc1,e1_2+enuc2+ejxc2)
            else:
                atm_exc, atm_ej = xceda.get_atmexc(self,atm2bas_p)
            atm_E = atm_e1 + atm_enuc + atm_ej + atm_exc
        if 'charge' in self.method:
            if self.showinter:
                bg_corrxn, bg_corrxn_fake, bg2,bg3 = get_bg_corrxn(self, atm2bas_f, 'charge')
                #RC_inter = eda_inter.get_RC_inter(bg1,bg2)
                logger.log(self.stdout_inter, "bg2",bg2)
                #logger.mlog(self.stdout_inter, "bg3",bg3)
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
        inter_terms = []
        if self.showinter: 
            logger.log(self.stdout_inter, "RR1",RR1)
            logger.log(self.stdout_inter, "RR2",RR2)
            logger.mlog(self.stdout_inter, "RR3",RR3)
            logger.mlog(self.stdout_inter, "RR4",RR4)
            inter_terms = [RR1, RR2, RR3, RR4] 
            if 'charge' in self.method:
                inter_terms.append(bg2) 
                inter_terms.append(bg3) 
        return atm_E, totE, conv, inter_terms

def build(eda, gjf, method):

    starttime = time.time()
    #xyznam = sys.argv[1]
    mol = gto.Mole()
    ##with open(xyznam) as f:
    #   geom = f.read()
    mol.atom, coords, charges, mol.charge, mol.spin = gjf_kit.gjf_parser(gjf)
    if 'charge' in method:
        eda.bgchgs = (coords, charges)
        if eda.molchgs is None:
            logger.slog(eda.stdout, "Warning: center frag has no atomic charges")
    if 'qmmm' in method:
        eda.bgchgs = (coords, charges)
    mol.cart=('cart' in method)
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
    bg correction (nuc)  = Z_A* Q_B/R_AB
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
        #bg1 = np.zeros(mol.natm)
        bg2 = np.zeros((mol.natm, mol.natm))
        bg3 = 0.0 #np.zeros((mol.natm, mol.natm, mol.natm))
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
        return bg_corrxn, bg_corrxn_fake, bg2 , bg3
    else:
        return bg_corrxn, bg_corrxn_fake

def get_E1(eda, atm2bas):
    #atom_kinE = []
    mol = eda.mol
    #method = eda.method
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
                    kin1[a] += kin
                else:
                    kin2[a,b] += kin
                    kin2[b,a] += kin
    atom_kin = atom_kin/2
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_kinE=",atom_kin)
    #if eda.showinter:
    #    logger.log(eda.stdout_inter,"kin1",kin1)
    #    logger.log(eda.stdout_inter,"kin2",kin2)

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
        if mol.cart:
            int1en = moleintor.getints("int1e_nuc_cart",fakeatm[i],mol._bas, mol._env)
        else:
            int1en = moleintor.getints("int1e_nuc_sph",fakeatm[i],mol._bas, mol._env)
        int1enuc.append(int1en)
    #e1 = 0.0
    #atom_h1E = np.zeros(atom_number)
    if eda.showinter:
        atom_1enuc, e1_1, e1_2, e1_3 = h1e_inter(dm,bas2atm,int1enuc,mol.natm,nao)
        atom_1enuc = np.asarray(atom_1enuc)[0:mol.natm]
        e1_1 = np.asarray(e1_1) + kin1
        e1_2 = np.triu(np.asarray(e1_2) + kin2)
        #intersum = e1_1.sum() + e1_2.sum() + e1_3
        #e1_3 = np.asarray(e1_3)
    else:
        atom_1enuc = np.asarray(h1e(dm,bas2atm,int1enuc,mol.natm,nao))[0:mol.natm]
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_1enuc=",atom_1enuc)
    atm_e1 = atom_kin + atom_1enuc
    logger.log(eda.stdout,"atom_e1=",atm_e1)
    if eda.showinter:
        logger.log(eda.stdout_inter,"e1_1",e1_1)
        logger.log(eda.stdout_inter,"e1_2",e1_2)
        logger.mlog(eda.stdout_inter,"e1_3 ",e1_3)
        #e1sum = atm_e1.sum()
        #logger.mlog(eda.stdout_inter,"intersum ", intersum)
        #logger.mlog(eda.stdout_inter,"e1sum ", e1sum)
    #with open(eda.output+'-eda.log','a') as f:
    if eda.anal:
        tot_akin = atom_kin.sum()
        kin_err = tot_akin - np.einsum('ij,ji',dm,mol.intor_symmetric('int1e_kin'))
        tot_a1enuc = atom_1enuc.sum()
        a1enuc_err = tot_a1enuc - np.einsum('ij,ji',dm,moleintor.getints("int1e_nuc_sph",mol._atm,mol._bas, mol._env))
        logger.mlog(eda.stdout,"kin_err: ",kin_err)
        logger.mlog(eda.stdout,"1enuc_err: ",a1enuc_err)
    if eda.showinter: atm_e1 = atm_e1, e1_1, e1_2 , e1_3
    return atm_e1

def get_Enuc(eda):
    mol = eda.mol
    charges = mol.atom_charges()
    if True: #eda.showinter:
        #enuc1 = np.zeros(mol.natm)
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
    if eda.showinter:
        for i in range(mol.natm):
            for j in range(i+1, mol.natm):
                enuc2[i][j] = charges[i]*charges[j]/rr[i][j]
        logger.log(eda.stdout_inter,"enucnuc_2",enuc2)
    if eda.anal:
        tot_enucnuc = atm_enucnuc.sum()
        enucnuc_err = tot_enucnuc - np.einsum('i,ij,j', charges, 1./rr, charges) * .5
        logger.mlog(eda.stdout, "err_enucnuc", enucnuc_err)
    atm_enucnuc = atm_enucnuc, enuc2
    return atm_enucnuc
