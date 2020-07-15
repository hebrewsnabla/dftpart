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


#EDATHRESH = 1e-8

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
        self.showinter = False
        self.exclude_cap = True
        self.inter_thresh = 1e-7
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
        if self.showinter:
            atm_e1, e1_1, e1_2, e1_3 = get_E1(self)
            atm_enuc, enuc1, enuc2 = get_Enuc(self)
        else:
            atm_e1 = get_E1(self)
            atm_enuc = get_Enuc(self)
        t2 = time.time()
        #with open(self.output+'-eda.log','a') as f
        logger.slog(self.stdout,"time for E1, E_nuc: %.5f\n", (t2-t1))
        if self.method[0] == 'hf':
            if self.showinter:
                atm_ej, atm_ek, ejk1, ejk2, ejk3, ejk4 = jkeda.get_Ejk(self, 'jk', self.jktype)
                atm_ejk = atm_ej + atm_ek
                RR1 = e1_1 + enuc1 + ejk1
                RR2 = e1_2 + enuc2 + ejk2
                RR3 = ejk3.merge(e1_3).cut(self.inter_thresh)
                RR4 = ejk4.cut(self.inter_thresh)
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
                bg_corrxn, bg_corrxn_fake, bg2,bg3, bgbg2 = get_bg_corrxn(self, 'charge')
                #RC_inter = eda_inter.get_RC_inter(bg1,bg2)
                #logger.log(self.stdout_inter, "bg2",bg2)
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
            conv = (totE - scfE) < self.conv_thresh
        else:
            totE_fake = (atm_E - bg_corrxn + bg_corrxn_fake).sum()
            logger.slog(self.stdout,"fake tot Energy =%16.10f", (totE_fake))
            logger.slog(self.stdout,"Err of fake totE =%16.10f", (totE_fake - scfE))
            conv = (totE_fake - scfE) < self.conv_thresh
        ok = ['FAIL','OK']
        logger.slog(self.stdout,"EDA convergence -- "+ok[conv])

        for i in range(atm_E.shape[0]):
            logger.slog(self.stdout,"%s %i %16.10f", self.mol.atom_symbol(i),i+1,atm_E[i])
        inter_terms = []
        if self.showinter: 
            RR1 = misc.mat2dict(RR1, self.inter_thresh, self.frag2layer)
            RR2 = misc.mat2dict(RR2, self.inter_thresh, self.frag2layer)
            logger.ilog(self.stdout_inter, "RR1",RR1)
            logger.ilog(self.stdout_inter, "RR2",RR2)
            logger.ilog(self.stdout_inter, "RR3",RR3)
            logger.ilog(self.stdout_inter, "RR4",RR4)
            inter_terms = [RR1, RR2, RR3, RR4] 
            if 'charge' in self.method:
                RC2 = misc.mat2dict(bg2, self.inter_thresh, self.frag2layer)
                RC3 = bg3.cut(self.inter_thresh)
                logger.ilog(self.stdout_inter, "RC2",RC2)
                logger.ilog(self.stdout_inter, "RC3",RC3)
                inter_terms.append(RC2) 
                inter_terms.append(RC3)
        return atm_E, totE, conv, inter_terms

    def get_frags(self):
        subsys_lso, num_subsys, num_subsys_tot, frg_intot, num_atoms = gfea3.lso_parser(self.lso)
        frags_list = []
        for cenfrg in range(1,len(frg_intot)+1):
            f = gfea3.Frag()
            f.atm_intot = frg_intot[cenfrg-1]
            f.atm_insub = frg_intot[cenfrg-1]
            f.label = cenfrg
            f.layer = 'c'
            frags_list.append(f)
        self.totnum_frag = len(frg_intot)
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
    #dm = eda.dm
    nao = eda.nao
    bas2atm = eda.bas2atm #get_bas2atm(atm2bas,nao,mol.natm)
    if eda.showinter:
        bas2frg = eda.bas2frg
        atm2frg = eda.atm2frg
    #print(bas2atm)
    kinmat = mol.intor_symmetric('int1e_kin')
    atom_kin = np.zeros(mol.natm)
    if eda.showinter:
        #kin1 = np.zeros(mol.natm)
        #kin2 = np.zeros((mol.natm, mol.natm))
        kin1 = np.zeros(eda.nfrag + 2)
        kin2 = np.zeros((eda.nfrag + 2, eda.nfrag +2))
    naorange = range(nao)
    #if eda.exclude_cap:
    #    naorange -= eda.capbas_p
    for i in naorange:
        for j in naorange:
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
    if eda.showinter:
        logger.log(eda.stdout_inter,"kin1",kin1)
        logger.log(eda.stdout_inter,"kin2",kin2)

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
    #atom_h1E = np.zeros(atom_number)
    if eda.showinter:
        #print(atm2frg,mol.natm)
        #print(int1enuc[0][:3,:3])
        atom_1enuc, e1_1, e1_2, e1_3 = h1e_inter(eda.dm,bas2atm, bas2frg, atm2frg, int1enuc,mol.natm,nao,eda.nfrag+2)
        atom_1enuc = np.asarray(atom_1enuc)[0:mol.natm]
        logger.log(eda.stdout, "e1n_1", e1_1)
        e1_1 = np.asarray(e1_1) + kin1
        #print(e1_2)
        #print(kin1,kin2)
        e1_2 = np.triu(np.asarray(e1_2) + kin2)
        #intersum = e1_1.sum() + e1_2.sum() + e1_3
        #e1_3 = np.asarray(e1_3)
        #print(e1_3)
        e1_3 = eda_inter.simp3(e1_3, eda.nfrag, eda.frag2layer)
    else:
        atom_1enuc = np.asarray(h1e(eda.dm,bas2atm,int1enuc,mol.natm,nao))[0:mol.natm]
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_1enuc=",atom_1enuc)
    atm_e1 = atom_kin + atom_1enuc
    logger.log(eda.stdout,"atom_e1=",atm_e1)
    if eda.showinter:
        logger.log(eda.stdout_inter,"e1_1",e1_1)
        logger.log(eda.stdout_inter,"e1_2",e1_2)
        logger.ilog(eda.stdout_inter,"e1_3 ",e1_3)
    #with open(eda.output+'-eda.log','a') as f:
    #anal = True
    if eda.anal and eda.showinter:
        tot_akin = atom_kin.sum()
        tot_kin = np.einsum('ij,ji',eda.dm,mol.intor_symmetric('int1e_kin'))
        tot_fkin = kin1.sum() + np.triu(kin2).sum()
        tot_a1enuc = atom_1enuc.sum()
        tot_1enuc = np.einsum('ij,ji',eda.dm,moleintor.getints("int1e_nuc_sph",mol._atm,mol._bas, mol._env))
        tot_fe1 = e1_1.sum() + e1_2.sum() + sum(e1_3.energies())

        '''basis_range = range(53)
        dm1 = eda.dm[np.ix_(basis_range, basis_range)]
        inte1n1 = np.zeros((nao,nao))
        for i in range(7):
            inte1n1 += int1enuc[i]
        inte1n1 = inte1n1[np.ix_(basis_range, basis_range)]
        e1n1 = np.einsum('ij,ji', dm1, inte1n1)
        logger.mlog(eda.stdout, "old_e1_1", e1n1)'''

        logger.mlog(eda.stdout,"a_e1 ", tot_akin + tot_a1enuc)
        #logger.mlog(eda.stdout,"a_e1k ", tot_akin)
        #logger.mlog(eda.stdout,"a_e1n ", tot_a1enuc)
        logger.mlog(eda.stdout,"f_e1 ",tot_fe1)
        logger.mlog(eda.stdout,"t_e1", tot_kin + tot_1enuc)
    if eda.showinter: atm_e1 = atm_e1, e1_1, e1_2 , e1_3
    return atm_e1


def get_Enuc(eda):
    mol = eda.mol
    charges = mol.atom_charges()
    if eda.showinter:
        enuc1 = np.zeros(eda.nfrag+2)
        enuc2 = np.zeros((eda.nfrag+2, eda.nfrag+2))
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
            if f.layer == 'q': continue
            for g in eda.frag_list:
                if g.layer == 'q': continue
                if g.label < f.label:
                    continue
                else:
                    if g.label == f.label:
                        tem = 0.0
                        for i in f.atm_insub:
                            for j in f.atm_insub:
                                if j > i:
                                    tem += charges[i-1]*charges[j-1]/rr[i-1][j-1]
                        enuc1[f.label-1] = tem
                    else:
                        tem = 0.0
                        for i in f.atm_insub:
                            for j in g.atm_insub:
                                tem += charges[i-1]*charges[j-1]/rr[i-1][j-1]
                        enuc2[f.label-1, g.label-1] = tem
        logger.log(eda.stdout_inter,"enucnuc_1",enuc1)
        logger.log(eda.stdout_inter,"enucnuc_2",enuc2)
        if eda.anal:
            tenuc = np.einsum('i,ij,j', charges, 1./rr, charges) * .5 
            #aenuc = atm_enucnuc.sum()
            fenuc = enuc1.sum() + enuc2.sum()
            logger.slog(eda.stdout, "tenuc: %f", tenuc)
            #logger.slog(eda.stdout_inter, "aenuc: %f", aenuc)
            logger.slog(eda.stdout, "fenuc: %f", fenuc)

            aenuc = atm_enucnuc.sum()
            logger.slog(eda.stdout, "aenuc: %f", aenuc)
            enucnuc_err = aenuc - tenuc
            logger.mlog(eda.stdout, "Err_aenuc", enucnuc_err)
    if eda.showinter:
        atm_enucnuc = atm_enucnuc, enuc1, enuc2
    return atm_enucnuc


