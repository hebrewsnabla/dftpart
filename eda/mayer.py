'''

'''
# PySCF module
from pyscf import gto, scf, dft, symm, qmmm
from pyscf.gto import moleintor
from pyscf.gto.mole import inter_distance
from pyscf.scf import _vhf

# dftpart module
from .new_eda import preri
from .h1e_new import h1e, bgh1e, h1e_inter, bgh1e_inter
from . import xceda, jkeda, eda_inter, bg
from ..kit import logger, dft_kit, gjf_kit, misc
from ..gfea import gfea3

# others
import numpy as np
import os,sys,subprocess 
import time
import copy
from functools import partial

einsum = partial(np.einsum, optimize=True)
np.set_printoptions(precision=6, linewidth=160, suppress=True)
#EDATHRESH = 1e-8

class MayerEDA():
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
                 test-inter.log : eda inter energies (optional)

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
        self.showinter = True
        self.exclude_cap = True
        self.inter_thresh = 1e-7
        self.inter_print = 1e-4
        self.conv_thresh = 1e-7
        ### GFEA ####
        self.molchgs = None # mol charges
        #self.envchgs = None # env frag charges
        self.bgchgs = None # background charges
        self.cen = None
        self.env = None # That's NOT mol.env_ !
        self.frag_list = None
        self.lso = None
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
        #if self.showinter:
        atm_e1, e1_2 = get_E1(self)
        atm_enuc, enuc2 = get_Enuc(self)
        #else:
        #    atm_e1 = get_E1(self)
        #    atm_enuc = get_Enuc(self)
        t2 = time.time()
        #with open(self.output+'-eda.log','a') as f
        logger.slog(self.stdout,"time for E1, E_nuc: %.5f\n", (t2-t1))
        if self.method[0] == 'hf':
            #if self.showinter:
            ejk2 = get_Ejk(self, 'jk', self.jktype)
            #atm_ejk = atm_ej + atm_ek
            #RR1 = e1_1 + enuc1 + ejk1
            RR2 = e1_2 + enuc2 + ejk2
            #print(RR2)
            logger.slog(self.stdout,"e1_2: %.5f",e1_2.sum())
            logger.slog(self.stdout,"enuc2: %.5f",enuc2.sum())
            logger.slog(self.stdout,"ejk2: %.5f", ejk2.sum())
            logger.log(self.stdout,"RR2",RR2)
            #RR3 = ejk3.merge(e1_3).cut(self.inter_thresh)
            #RR4 = ejk4.cut(self.inter_thresh)
            #RR_inter = eda_inter.get_RR_inter(e1_1+ejk1,e1_2+enuc2+ejk2)
            #else:
            #    atm_ej, atm_ek = jkeda.get_Ejk(self, 'jk', self.jktype)
            #    atm_ejk = atm_ej + atm_ek
            #atm_E = atm_e1 + atm_enuc + atm_ejk
            t3 = time.time()
            #with open(self.output+'-eda.log','a') as f:
            logger.slog(self.stdout,"time for Ej, Ek: %.5f\n", (t3-t2))
        elif dft_kit.is_dft(self.method[0]):
            #if self.showinter:
            atm_exc, atm_ej, ejxc1, ejxc2 = xceda.get_atmexc(self)
            RR1 = e1_1 + ejxc1
            RR2 = e1_2 + enuc2 + ejxc2
            #    RR3 = e1_3.merge(ejxc3)
            #    RR4 = ejxc4
            #    #RR_inter = eda_inter.get_RR_inter(e1_1+ejxc1,e1_2+enuc2+ejxc2)
            #else:
            #    atm_exc, atm_ej = xceda.get_atmexc(self)
            atm_E = atm_e1 + atm_enuc + atm_ej + atm_exc
        if 'charge' in self.method:
            #if self.showinter:
            bg_corrxn, bg_corrxn_fake, bg2,bg3, bgbg2 = get_bg_corrxn(self, 'charge')
            #RC_inter = eda_inter.get_RC_inter(bg1,bg2)
            #logger.log(self.stdout_inter, "bg2",bg2)
            #logger.mlog(self.stdout_inter, "bg3",bg3)
            #else:
            #    bg_corrxn, bg_corrxn_fake = get_bg_corrxn(self, 'charge')
            atm_E += bg_corrxn
        elif 'qmmm' in self.method:
            # show inter TBD
            bg_corrxn, bg_corrxn_fake = get_bg_corrxn(self, 'qmmm')
            atm_E += bg_corrxn
        #atm_ehf = atm_e1 + atm_enuc + atm_ej
        #atm_E = anal(self, atm_ek, atm_exc)
        totE = np.sum(RR2)
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
            conv = abs(totE - scfE) < self.conv_thresh
        else:
            totE_fake = (atm_E - bg_corrxn + bg_corrxn_fake).sum()
            logger.slog(self.stdout,"fake tot Energy =%16.10f", (totE_fake))
            logger.slog(self.stdout,"Err of fake totE =%16.10f", (totE_fake - scfE))
            conv = abs(totE_fake - scfE) < self.conv_thresh
        ok = ['FAIL','OK']
        logger.slog(self.stdout,"EDA convergence -- "+ok[conv])

        exit()
        for i in range(atm_E.shape[0]):
            logger.slog(self.stdout,"%s %i %16.10f", self.mol.atom_symbol(i),i+1,atm_E[i])
        inter_terms = []
        if self.showinter: 
            RR1 = misc.mat2dict(RR1, self.inter_thresh, self.frag2layer)
            RR2 = misc.mat2dict(RR2, self.inter_thresh, self.frag2layer)
            #RR1p = RR1.cut(self.inter_print)
            RR2p = RR2.cut(self.inter_print)
            #RR3p = RR3.cut(self.inter_print)
            #RR4p = RR4.cut(self.inter_print)
            RR1E = sum(RR1.energies())
            RR2E = sum(RR2.energies())
            #RR3E = sum(RR3.energies())
            #RR4E = sum(RR4.energies())
            logger.mlog(self.stdout, 'RR1E, RR2E', [RR1E, RR2E])
            logger.ilog(self.stdout_inter, "RR1",RR1)
            logger.ilog(self.stdout_inter, "RR2",RR2p)
            #logger.ilog(self.stdout_inter, "RR3",RR3p)
            #logger.ilog(self.stdout_inter, "RR4",RR4p)
            inter_terms = [RR1, RR2] 
            if 'charge' in self.method:
                RC2 = misc.mat2dict(bg2, self.inter_thresh, self.frag2layer)
                RC3 = bg3.cut(self.inter_thresh)
                logger.ilog(self.stdout_inter, "RC2",RC2)
                logger.ilog(self.stdout_inter, "RC3",RC3)
                inter_terms.append(RC2) 
                inter_terms.append(RC3)
        return atm_E, totE, conv, inter_terms

    def get_frags(self):
        frags_list = []
        if self.lso is not None:
            subsys_lso, num_subsys, num_subsys_tot, frg_intot, num_atoms = gfea3.lso_parser(self.lso)
            for cenfrg in range(1,len(frg_intot)+1):
                f = gfea3.Frag()
                f.atm_intot = frg_intot[cenfrg-1]
                f.atm_insub = frg_intot[cenfrg-1]
                f.label = cenfrg
                f.layer = 'c'
                frags_list.append(f)
            self.totnum_frag = len(frg_intot)
        else:
            for i in range(self.mol.natm):
                f = gfea3.Frag()
                f.atm_intot = [i+1] 
                f.atm_insub = [i+1] 
                f.label = i+1
                f.layer = 'c'
                frags_list.append(f)
        self.frag_list = frags_list
        


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
    eda.cart = mol.cart

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
    if eda.showinter:
        if eda.frag_list is None:
            eda.get_frags()
        eda.bas2frg = get_bas2frg(eda.bas2atm, eda.frag_list)           # frg starts from 1
        eda.atm2frg = get_atm2frg(mol.natm, eda.frag_list)
        eda.nfrag = len(eda.frag_list)
        eda.ncap = 0
        eda.capatoms = []
        eda.frag2layer = {}
        for f in eda.frag_list:
            eda.frag2layer[f.label] = f.layer
            if f.layer == 'cap': 
                eda.ncap += 1
                eda.capatoms.append(f.atm_insub[0])
        eda.capbas_p = []
        eda.capbas_f = []
        for i in eda.capatoms:
            eda.capbas_p.append(eda.atm2bas_p[i-1])
            eda.capbas_f.append(eda.atm2bas_f[i-1])
    if eda.verbose >= 6:
        logger.mlog(eda.stdout, "atm2bas_f", eda.atm2bas_f)
        logger.mlog(eda.stdout, "atm2bas_p", eda.atm2bas_p)
        #print(eda.bas2atm)
        logger.mlog(eda.stdout, "bas2atm", eda.bas2atm)
        if eda.showinter:
            logger.mlog(eda.stdout, "bas2frg", eda.bas2frg)
            logger.mlog(eda.stdout, "atm2frg", eda.atm2frg)
            logger.mlog(eda.stdout, "cap atoms", eda.capatoms)
            logger.mlog(eda.stdout, "cap basis_p", eda.capbas_p)
      
    eda.built = True
    return eda

def get_atm2frg(natm, frag_list):
    atm2frg = []
    for a in range(natm):
        for f in frag_list:
            if f.layer == 'q':
                continue
            if (a+1 in f.atm_insub):
                atm2frg.append(f.label)
    return atm2frg

def get_bas2frg(bas2atm, frag_list):
    bas2frg = []
    for atm in bas2atm:
        for f in frag_list:
            if f.layer == 'q': 
                continue
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
    nao = eda.nao
    atm2bas = eda.atm2bas_f
    bas2atm = eda.bas2atm
    if eda.showinter:
        bas2frg = eda.bas2frg
    molchgs = eda.molchgs
    bgcoords, bgchgs = eda.bgchgs
    #cen = eda.cen
    #env = eda.env
    bg_corrxn = np.zeros(mol.natm)
    if eda.showinter:
        #bg1 = np.zeros(mol.natm)
        #bg2 = np.zeros((mol.natm, mol.natm))
        #bg3 = np.zeros((mol.natm, mol.natm, mol.natm))
        elecbg2 = np.zeros((eda.nfrag+2, eda.nfrag+2))
        nucbg2 = np.zeros((eda.nfrag+2, eda.nfrag+2))
        bg2 = np.zeros((eda.nfrag+2, eda.nfrag+2))
        bg3 = misc.EDict()
        bgbg2 = np.zeros((eda.nfrag+2, eda.nfrag+2))
    vchg = bg.inter_elecbg(mol, dm, bgcoords, bgchgs)
    #bas2atm = get_bas2atm(atm2bas,nao,mol.natm)
    if charge=='charge':
        elecbg = bgh1e(dm,bas2atm,vchg,mol.natm,nao)
        atom_elecbg = 2*np.asarray(elecbg)[0:mol.natm]
        atom_elecbg -= 0.5*bg.inter_bgbg(mol.atom_coords(), molchgs, bgcoords, bgchgs)
        atom_nucbg = bg.inter_nucbg(mol, bgcoords, bgchgs)
        if eda.showinter:
            for f in eda.frag_list:
                if f.layer is not 'q': 
                    continue
                chg_insub = misc.one2zero(f.atm_insub)
                vchg_frag = bg.inter_elecbg(mol, dm, bgcoords[chg_insub], bgchgs[chg_insub])
                felecbg_a, felecbg2, felecbg3 = bgh1e_inter(dm,bas2atm, bas2frg, vchg_frag, mol.natm,nao,eda.nfrag+2)
                #print(felecbg2)
                for g in eda.frag_list:
                    if g.layer is 'q': 
                        continue
                    elecbg2[g.label-1,f.label-1] = felecbg2[g.label-1]
                    nucbg2[g.label-1,f.label-1] = bg.inter_nucbg_f(mol, g.atm_insub, bgcoords, bgchgs, f.atm_insub)
                    for h in eda.frag_list:
                        if (h.layer is 'q') or (h.label <= g.label):
                            continue
                        bg3gh = felecbg3[g.label-1,h.label-1]
                        if abs(bg3gh) > 1e-12:
                            ghf = (g.label,h.label,f.label)
                            layers = (g.layer, h.layer, f.layer)
                            if 'cap' in layers:
                                continue
                            bg3[ghf] = [bg3gh, layers]
            bg2 = elecbg2 + nucbg2
            #bg2 = misc.mat2dict(bg2, eda.inter_thresh, eda.frag2layer)
            #bgbg, bgbg2 = bg.inter_bgbg(mol.atom_coords(), molchgs, bgcoords, bgchgs)
            #atom_elecbg_f = 2*np.asarray(elecbg)[0:mol.natm] - 0.5*bgbg
            #bg3 = eda_inter.simp3(bg3)
            #logger.mlog(eda.stdout_inter, "elecbg_a", atom_elecbg.sum())
            #logger.mlog(eda.stdout_inter, "elecbg_f", atom_elecbg_f.sum())
            if eda.verbose >= 6:
                logger.log(eda.stdout_inter, "elecbg2", elecbg2)
                logger.log(eda.stdout_inter, "nucbg2", nucbg2)
            logger.log(eda.stdout_inter, "bg2", bg2)
            logger.ilog(eda.stdout_inter, "bg3", bg3)
            #logger.mlog(eda.stdout_inter, "bgbg2", bgbg2)

    elif charge=='qmmm':
        elecbg = bgh1e(dm,bas2atm,vchg,mol.natm,nao)
        atom_elecbg = np.asarray(elecbg)[0:mol.natm]
        atom_nucbg = 0.5*bg.inter_nucbg(mol, bgcoords, bgchgs)
    
    logger.log(eda.stdout,"atom_elecbg=",atom_elecbg)
    logger.log(eda.stdout,"atom_nucbg=",atom_nucbg)

    bg_corrxn = atom_elecbg + atom_nucbg
    logger.log(eda.stdout,"atom_bg_correction=",bg_corrxn)
    bg_corrxn_fake = 2*np.asarray(elecbg)[0:mol.natm] + bg.inter_nucbg(mol, bgcoords, bgchgs)
    if eda.showinter:
        return bg_corrxn, bg_corrxn_fake, bg2 , bg3, bgbg2
    else:
        return bg_corrxn, bg_corrxn_fake

def get_E1(eda):
    #atom_kinE = []
    mol = eda.mol
    #method = eda.method
    dm = eda.dm
    nao = eda.nao
    #atm2bas = eda.atm2bas_p
    bas2atm = eda.bas2atm #get_bas2atm(atm2bas,nao,mol.natm)
    #if eda.showinter:
    #    bas2frg = eda.bas2frg
    #    atm2frg = eda.atm2frg
    #print(bas2atm)
    kinmat = mol.intor_symmetric('int1e_kin')
    #atom_kin = np.zeros(mol.natm)
    if eda.showinter:
        #kin1 = np.zeros(mol.natm)
        #kin2 = np.zeros((mol.natm, mol.natm))
        #kin1 = np.zeros(eda.nfrag)
        kin2 = np.zeros((mol.natm, mol.natm))
        e1_2 = np.zeros((mol.natm, mol.natm))
    #naorange = range(nao)
    #if eda.exclude_cap:
    #    naorange -= eda.capbas_p
    aoslice = gto.aoslice_by_atom(mol)
    for i in range(mol.natm):
        #kin1[i] = 
        ao0i, ao1i = aoslice[i][2:]
        for j in range(i, mol.natm):
            ao0j, ao1j = aoslice[j][2:]
            kin2[i,j] = einsum('ij,ji->', dm[ao0i:ao1i, ao0j:ao1j], kinmat[ao0j:ao1j, ao0i:ao1i])
    #with open(eda.output+'-eda.log','a') as f:
    #print(kin2)
    #logger.log(eda.stdout,"atom_kinE=",atom_kin)
    kin2 = kin2 + np.triu(kin2,1)
    if eda.showinter:
        #logger.log(eda.stdout_inter,"kin1",kin1)
        logger.log(eda.stdout,"kin2",kin2)

    fakeatm = []
    for n in range(mol.natm):
        #if eda.exclude_cap and (n+1 in eda.capatoms):
        #    continue
        back = copy.copy(mol._atm)
        for i in range(mol.natm):
            if (i!=n):
                back[i,0] = 0
        fakeatm.append(back)
    #fakeatm = np.array(fakeatm)
    #print(mol._atm)
    #print(fakeatm)
    int1enuc = []
    for i in range(mol.natm):
        #if eda.exclude_cap and (n+1 in eda.capatoms):
        #    #int1en = 
        #    continue
        if mol.cart:
            int1en = moleintor.getints("int1e_nuc_cart",fakeatm[i],mol._bas, mol._env)
        else:
            int1en = moleintor.getints("int1e_nuc_sph",fakeatm[i],mol._bas, mol._env)
        int1enuc.append(int1en)
    #e1 = 0.0
    for i in range(mol.natm):
        #kin1[i] = 
        ao0i, ao1i = aoslice[i][2:]
        for j in range(mol.natm):
            e1_2[i,j] = einsum('ij,ji->', dm[:,ao0i:ao1i], int1enuc[j][ao0i:ao1i, :])
    e1_2 = np.triu(e1_2) + np.tril(e1_2, -1).T 
    #print(e1_2)
    e1_2 += kin2
    #with open(eda.output+'-eda.log','a') as f:
    #logger.log(eda.stdout,"atom_1enuc=",atom_1enuc)
    atm_e1 = np.zeros(mol.natm)
    #logger.log(eda.stdout,"atom_e1=",atm_e1)
    if eda.showinter:
        #logger.log(eda.stdout_inter,"e1_1",e1_1)
        logger.log(eda.stdout,"e1_2",e1_2)
        #logger.ilog(eda.stdout_inter,"e1_3 ",e1_3)
    #with open(eda.output+'-eda.log','a') as f:

    return atm_e1, e1_2


def get_Enuc(eda):
    mol = eda.mol
    charges = mol.atom_charges()
    if eda.showinter:
        #enuc1 = np.zeros(eda.nfrag+2)
        enuc2 = np.zeros((mol.natm, mol.natm))
    #coords = mol.atom_coords()
    rr = inter_distance(mol)
    #print(rr)
    rr[np.diag_indices_from(rr)] = 1e200
    #if CHECK_GEOM and numpy.any(rr < 1e-5):
    #    for atm_idx in numpy.argwhere(rr<1e-5):
    #    logger.warn(mol, 'Atoms %s have the same coordinates', atm_idx)
    #    raise RuntimeError('Ill geometry')
    atm_enucnuc = einsum('i,ij,j->i', charges, 1./rr, charges) * .5
    logger.log(eda.stdout,"atom_enucnuc=",atm_enucnuc)
    enuc2 = einsum('i,ij,j->ij', charges, 1./rr, charges)
    enuc2 = np.triu(enuc2)
    #print(enuc2)
    logger.log(eda.stdout,"enuc2",enuc2)

    return atm_enucnuc, enuc2

def get_Ejk(eda, jk='jk', jktype='py'):
    mol = eda.mol
    dm = eda.dm
    aoslice = gto.aoslice_by_atom(mol)  
    atm2bas = eda.atm2bas_p
    ej_2 = np.zeros((mol.natm, mol.natm))
    ek_2 = np.zeros((mol.natm, mol.natm))
    for i in range(mol.natm):
        sh0i, sh1i, ao0i, ao1i = aoslice[i]
        for j in range(i, mol.natm):
            sh0j, sh1j, ao0j, ao1j = aoslice[j]
            vj = scf.jk.get_jk(mol, dm[:,ao0i:ao1i], scripts='ijkl,ji->kl', intor='int2e',
                    shls_slice=(sh0i,sh1i,0,mol.nbas,sh0j,sh1j,0,mol.nbas))
            vk = scf.jk.get_jk(mol, dm[:,ao0i:ao1i], scripts='ijkl,li->kj', intor='int2e',
                    shls_slice=(sh0i,sh1i,0,mol.nbas,sh0j,sh1j,0,mol.nbas))
            ej = einsum('ij,ji->', dm[:, ao0j:ao1j],vj ) * 0.5
            ek = einsum('ij,ji->', dm[:, ao0j:ao1j],vk ) * (-0.25)
            ej_2[i,j] = ej
            ek_2[i,j] = ek

    ejk2 = ej_2 + ek_2
    ejk2 = ejk2 + np.triu(ejk2, 1)
    #print(ejk2)
    logger.log(eda.stdout,"ejk2",ejk2)
    return ejk2

