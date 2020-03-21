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
        self.nao = None
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
        self.frag_list = []
        self.jktype = 'py'
        # mapping
        self.atm2bas_p = None
        self.atm2bas_f = None
        self.bas2atm = None
        self.bas2atm_f = None
        self.bas2frg = None

    def build(self):
        return build(self, self.gjf, self.method)

    def kernel(self):
        if self.built == False:
            self.build()
        t1 = time.time()
        #atm2bas_f, atm2bas_p = get_atm2bas(self.mol)
        if self.showinter:
            atm_e1, e1_1, e1_2 = get_E1(self)
        else:
            atm_e1 = get_E1(self)
        atm_enuc, enuc1, enuc2 = get_Enuc(self)
        t2 = time.time()
        #with open(self.output+'-eda.log','a') as f
        logger.slog(self.stdout,"time for E1, E_nuc: %.5f\n", (t2-t1))
        if self.method[0] == 'hf':
            if self.showinter:
                atm_ej, atm_ek, ejk1, ejk2, ejk3, ejk4 = jkeda.get_Ejk(self, 'jk', self.jktype)
                atm_ejk = atm_ej + atm_ek
                RR1 = e1_1 + enuc1 + ejk1
                RR2 = e1_2 + enuc2 + ejk2
                RR3 = ejk3
                RR4 = ejk4
                #RR_inter = eda_inter.get_RR_inter(e1_1+ejk1,e1_2+enuc2+ejk2)
            else:
                atm_ej, atm_ek = jkeda.get_Ejk(self, 'jk', self.jktype)
                atm_ejk = atm_ej + atm_ek
            atm_E = atm_e1 + atm_enuc + atm_ejk
            t3 = time.time()
            #with open(self.output+'-eda.log','a') as f:
            logger.slog(self.stdout,"time for Ej, Ek: %.5f\n", (t3-t2))
        elif dft_kit.is_dft(self.method[0]):
            if self.showinter:
                atm_exc, atm_ej, ejxc1, ejxc2, ejxc3, ejxc4 = xceda.get_atmexc(self)
                RR1 = e1_1 + ejxc1
                RR2 = e1_2 + enuc2 + ejxc2
                RR3 = e1_3 + ejxc3
                RR4 = ejxc4
                #RR_inter = eda_inter.get_RR_inter(e1_1+ejxc1,e1_2+enuc2+ejxc2)
            else:
                atm_exc, atm_ej = xceda.get_atmexc(self)
            atm_E = atm_e1 + atm_enuc + atm_ej + atm_exc
        if 'charge' in self.method:
            if self.showinter:
                bg_corrxn, bg_corrxn_fake, bg2,bg3 = get_bg_corrxn(self, 'charge')
                #RC_inter = eda_inter.get_RC_inter(bg1,bg2)
                logger.log(self.stdout_inter, "bg2",bg2)
                #logger.mlog(self.stdout_inter, "bg3",bg3)
            else:
                bg_corrxn, bg_corrxn_fake = get_bg_corrxn(self, 'charge')
            atm_E += bg_corrxn
        elif 'qmmm' in self.method:
            # show inter TBD
            bg_corrxn, bg_corrxn_fake = get_bg_corrxn(self, 'qmmm')
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
            '''logger.log(self.stdout_inter, "RR1",RR1)
            logger.log(self.stdout_inter, "RR2",RR2)
            logger.mlog(self.stdout_inter, "RR3",RR3)
            logger.mlog(self.stdout_inter, "RR4",RR4)
            inter_terms = [RR1, RR2, RR3, RR4] 
            if 'charge' in self.method:
                inter_terms.append(bg2) 
                inter_terms.append(bg3)''' 
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
    eda.nao = len(eda.dm)

    os.system("echo '' > "+eda.output+'-eda.log')
    eda.stdout = open(eda.output+'-eda.log','a')
    if eda.showinter:
        os.system("echo '' > "+eda.output+'-inter.log')
        eda.stdout_inter = open(eda.output+'-inter.log','a')
    #with open(self.output+'-eda.log','a') as f:
    logger.mlog(eda.stdout,"method,basis: ", eda.method)
    eda.atm2bas_f, eda.atm2bas_p = get_atm2bas(eda.mol)        
    # _p: bas starts from 0
    # _f: bas starts from 1 
    eda.bas2atm, eda.bas2atm_f = get_bas2atm(eda.atm2bas_f, eda.nao, eda.mol.natm) 
    #     atm starts from 0
    # _f: atm starts from 1
    eda.bas2frg = get_bas2frg(eda.bas2atm, eda.frag_list)           # frg starts from 1
    eda.atm2frg = get_atm2frg(mol.natm, eda.frag_list)
    eda.nfrag = len(eda.frag_list)
    eda.ncap = 0
    for f in eda.frag_list:
        if f.layer == 'cap': eda.ncap += 1

    if eda.verbose >= 6:
        logger.mlog(eda.stdout, "atm2bas_f", eda.atm2bas_f)
        logger.mlog(eda.stdout, "atm2bas_p", eda.atm2bas_p)
        #print(eda.bas2atm)
        logger.mlog(eda.stdout, "bas2atm", eda.bas2atm)
        logger.mlog(eda.stdout, "bas2frg", eda.bas2frg)
        logger.mlog(eda.stdout, "atm2frg", eda.atm2frg)
    
    eda.built = True
    return eda

def get_atm2frg(natm, frag_list):
    atm2frg = []
    for a in range(natm):
        for f in frag_list:
            if (a+1 in f.atm_insub):
                atm2frg.append(f.label)
    return atm2frg

def get_bas2frg(bas2atm, frag_list):
    bas2frg = []
    for atm in bas2atm:
        for f in frag_list:
            if (atm+1 in f.atm_insub):
                bas2frg.append(f.label)
    return bas2frg


def get_atm2bas(mol):
    atm2bas_p = []
    atm2bas_f = []
    #twoatom = []
    #threeatom = []
    #fouratom = []

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
    #print(atm2bas_f)
    #print(atm2bas_p)

    return atm2bas_f, atm2bas_p

def get_bas2atm(atm2bas,nao,natm):
    bas2atm = [] #np.zeros(nao, dtype=int)
    bas2atm_f = []
    for i in range(1,nao+1):
        for j in range(natm):
            if (i in atm2bas[j]):
                bas2atm.append( j ) # atm starts from 0
                bas2atm_f.append(j+1)
    return bas2atm, bas2atm_f


def get_bg_corrxn(eda, charge='charge'):
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
    atm2bas = eda.atm2bas_f
    bas2atm = eda.bas2atm
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
    #bas2atm = get_bas2atm(atm2bas,nao,mol.natm)
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

def get_E1(eda):
    #atom_kinE = []
    mol = eda.mol
    #method = eda.method
    #dm = eda.dm
    nao = eda.nao
    bas2atm = eda.bas2atm #get_bas2atm(atm2bas,nao,mol.natm)
    bas2frg = eda.bas2frg
    atm2frg = eda.atm2frg
    #print(bas2atm)
    kinmat = mol.intor_symmetric('int1e_kin')
    atom_kin = np.zeros(mol.natm)
    if eda.showinter:
        #kin1 = np.zeros(mol.natm)
        #kin2 = np.zeros((mol.natm, mol.natm))
        kin1 = np.zeros(eda.totnum_frag + 2)
        kin2 = np.zeros((eda.totnum_frag + 2, eda.totnum_frag +2))
    for i in range(nao):
        for j in range(nao):
            #print(dm,kin)
            kin = eda.dm[i,j] * kinmat[i,j]
            a = bas2atm[i]
            b = bas2atm[j]
            atom_kin[a] = atom_kin[a] + kin
            atom_kin[b] = atom_kin[b] + kin
            if eda.showinter:
                aa = bas2frg[i]
                bb = bas2frg[j]
                if aa==bb:
                    kin1[aa-1] += kin
                else:
                    kin2[aa-1,bb-1] += kin
                    kin2[bb-1,aa-1] += kin
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
        #print(atm2frg,mol.natm)
        atom_1enuc, e1_1, e1_2 = h1e_inter(eda.dm,bas2atm, bas2frg, atm2frg, int1enuc,mol.natm,nao,eda.totnum_frag+2)
        atom_1enuc = np.asarray(atom_1enuc)[0:mol.natm]
        e1_1 = np.asarray(e1_1) + kin1
        #print(e1_2)
        #print(kin1,kin2)
        e1_2 = np.triu(np.asarray(e1_2) + kin2)
        #intersum = e1_1.sum() + e1_2.sum() + e1_3
        #e1_3 = np.asarray(e1_3)
        #print(e1_3)
        #e1_3 = simp3(e1_3, eda.nfrag)
    else:
        atom_1enuc = np.asarray(h1e(eda.dm,bas2atm,int1enuc,mol.natm,nao))[0:mol.natm]
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_1enuc=",atom_1enuc)
    atm_e1 = atom_kin + atom_1enuc
    logger.log(eda.stdout,"atom_e1=",atm_e1)
    if eda.showinter:
        logger.log(eda.stdout_inter,"e1_1",e1_1)
        logger.log(eda.stdout_inter,"e1_2",e1_2)
        #logger.mlog(eda.stdout_inter,"e1_3 ",e1_3)
    #with open(eda.output+'-eda.log','a') as f:
    anal = True
    if anal:
        tot_akin = atom_kin.sum()
        #tot_kin = np.einsum('ij,ji',eda.dm,mol.intor_symmetric('int1e_kin'))
        tot_fkin = kin1.sum() + np.triu(kin2).sum()
        tot_a1enuc = atom_1enuc.sum()
        #tot_1enuc = np.einsum('ij,ji',eda.dm,moleintor.getints("int1e_nuc_sph",mol._atm,mol._bas, mol._env))
        tot_fe1 = e1_1.sum() + e1_2.sum()
        logger.mlog(eda.stdout,"a_e1 ", tot_akin + tot_a1enuc)
        logger.mlog(eda.stdout,"e1 ",tot_fe1)
    if eda.showinter: atm_e1 = atm_e1, e1_1, e1_2 #, e1_3
    return atm_e1

def simp3(e3, nfrag):
    e3simp = {}
    for f in range(nfrag):
        for g in range(f+1,nfrag):
            for h in range(g+1, nfrag):
                fgh = "%d,%d,%d"%(f,g,h)
                e3simp[fgh] = e3[f,g,h]
    return e3simp

def get_Enuc(eda):
    mol = eda.mol
    charges = mol.atom_charges()
    if True: #eda.showinter:
        enuc1 = np.zeros(eda.totnum_frag+2)
        enuc2 = np.zeros((eda.totnum_frag+2, eda.totnum_frag+2))
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
        for f in eda.frag_list:
            for g in eda.frag_list:
                if g.label < f.label:
                    continue
                else:
                    tem = 0.0
                    for i in f.atm_insub:
                        for j in g.atm_insub:
                            tem += charges[i-1]*charges[j-1]/rr[i-1][j-1]
                    if g.label == f.label:
                        enuc1[f.label] = tem
                    else:
                        enuc2[f.label, g.label] = tem
        logger.log(eda.stdout_inter,"enucnuc_1",enuc1)
        logger.log(eda.stdout_inter,"enucnuc_2",enuc2)
    if eda.anal:
        tot_enucnuc = atm_enucnuc.sum()
        enucnuc_err = tot_enucnuc - np.einsum('i,ij,j', charges, 1./rr, charges) * .5
        logger.mlog(eda.stdout, "err_enucnuc", enucnuc_err)
    atm_enucnuc = atm_enucnuc, enuc1, enuc2
    return atm_enucnuc
