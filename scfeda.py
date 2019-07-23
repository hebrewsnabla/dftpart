'''

'''

from pyscf import gto, scf, dft, symm, qmmm
#import edanew
from pyscf.scf import _vhf
#from frame_small2 import preri
#from frame_small6 import preri
from h1e_new import h1e
#from h1e import h1e
import xceda
#from frame_small5 import preri
from pyscf.gto import moleintor
from pyscf.gto.mole import inter_distance
import numpy as np
import os,sys,subprocess
import time
import copy
#from gfea import gfea2_2
from QCKit import logger, dft_kit, gjf_kit
class EDA():
    def __init__(self):
        self.mol = None
        self.mf = None 
	#eda.atm_ = None     
        #eda.atml = None
        #eda.bas_ = None 
        #eda.basl = None 
        #eda.env_ = None
        #eda.envl = None
        self.dm = None
        #eda.nao = None
        #eda.nbas = None
        #self.atomlist = None
        #self.realatomlabel = None
        # must specify
        self.gjf = None
        self.method = None
        #self.spinlist = None
        self.output = None
        self.verbose = 4
        self.built = False
    def build(self):
        return build(self, self.gjf, self.method)

    def kernel(self):
        if self.built == False:
            self.build()
        os.system("echo '' > "+self.output+'-eda.log')
        self.stdout = open(self.output+'-eda.log','a')
        #with open(self.output+'-eda.log','a') as f:
        logger.mlog(self.stdout,"method,basis: ", self.method)
        t1 = time.time()
        atm2bas_f, atm2bas_p = get_atm2bas(self.mol)
        atm_e1 = get_E1(self,atm2bas_f)
        atm_enuc = get_Enuc(self)
        t2 = time.time()
        #with open(self.output+'-eda.log','a') as f:
        logger.slog(self.stdout,"time for E1, E_nuc: %.5f\n", (t2-t1)) 
        if self.method[0] == 'hf':
            atm_ej, atm_ek = get_Ejk(self, atm2bas_p,'jk')
            atm_E = atm_e1 + atm_enuc + atm_ej + atm_ek
            t3 = time.time()
            #with open(self.output+'-eda.log','a') as f:
            logger.slog(self.stdout,"time for Ej, Ek: %.5f\n", (t3-t2))
        elif dft_kit.is_dft(self.method[0]):
            atm_exc, atm_ej = xceda.get_atmexc(self,atm2bas_p) 
            atm_E = atm_e1 + atm_enuc + atm_ej + atm_exc
        #atm_ehf = atm_e1 + atm_enuc + atm_ej
        #atm_E = anal(self, atm_ek, atm_exc)
        totE = np.sum(atm_E)
        #outfile = xyznam[0:-3] + "atom"
        if '(' in self.output:
            op = self.output.replace('(','\(')
            self.output = op.replace(')','\)')
        scfE = subprocess.getstatusoutput("grep 'converged SCF energy' "+self.output+'-pyscf.log')[1]
        scfE = float(scfE.strip().split()[-1])
        time.sleep(20)
        #with open(self.output+'-eda.log','a') as f:
        logger.slog(self.stdout,"tot Energy =%16.10f", (totE))
        logger.slog(self.stdout,"SCF Energy =%16.10f", (scfE))
        logger.slog(self.stdout,"Err of totE =%16.10f", (totE - scfE))
        if (totE - scfE) < 1e-8:
            conv = True
            logger.slog(self.stdout,'    OK')
        else:
            conv = False
            logger.slog(self.stdout,'    FAIL')
        for i in range(atm_E.shape[0]):
            logger.slog(self.stdout,"%s %i %16.10f", self.mol.atom_symbol(i),i+1,atm_E[i])
        #CH3E = atm_E[:4].sum()
        ##with open("CH3_anal.txt",'a') as f:
        #    logger.slog(self.stdout,"%s %.10f\n" % (self.output,CH3E))
        #endtime = time.time()
        #print('timeTot=',(endtime-starttime))
        return atm_E, totE, conv

def build(eda, gjf, method):
    
    starttime = time.time()
    #xyznam = sys.argv[1]
    mol = gto.Mole()
    ##with open(xyznam) as f:
    #   geom = f.read()
    mol.atom, coords, charges, charge, spin = gjf_kit.gjf_parser(gjf)
    
    #mol.cart=True
    mol.basis = method[1]
    #mol.symmetry = 1
    mol.output = eda.output + '-pyscf.log'
    mol.verbose = eda.verbose
    mol.build()
    eda.mol = mol
    #atom_number = len(mol._atom)
    #atom_coord = np.zeros((atom_number,3))
    #for i in range(atom_number):
    #    atom_coord[i] = mol.atom_coord(i)
    #eda.atomlist = np.ones(mol.natm,dtype=int)
    
    #eda.realatomlabel = list(range(mol.natm))
    
    if method[0] == 'hf':
        mf = scf.RHF(mol)
    elif dft_kit.is_dft(method[0]):
        mf = dft.RKS(mol)
        mf.xc = method[0]
    mf.kernel()
    eda.mf = mf
    
    #pyscf_time = time.time()
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


#b = time.time()
#print("combination-time=",b-pyscf_time)
def get_E1(eda, atm2bas):
    #atom_kinE = []
    mol = eda.mol
    dm = eda.dm
    nao = len(dm)
    bas2atm = get_bas2atm(atm2bas,nao,mol.natm)
    kinmat = mol.intor_symmetric('int1e_kin')
    atom_kin = np.zeros(mol.natm)
    for i in range(nao):
        for j in range(nao):
            #print(dm,kin)
            kin = dm[i,j] * kinmat[i,j]
            a = int(bas2atm[i])
            b = int(bas2atm[j])
            atom_kin[a] = atom_kin[a] + kin
            atom_kin[b] = atom_kin[b] + kin
    atom_kin = atom_kin/2
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_kinE=",atom_kin)
    
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
    atom_1enuc = np.asarray(h1e(dm,bas2atm,int1enuc,mol.natm,nao))[0:mol.natm]
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_1enuc=",atom_1enuc)
    atm_e1 = atom_kin + atom_1enuc
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_e1=",atm_e1)
    if eda.verbose > 5:
        tot_akin = atom_kin.sum()
        kin_err = tot_akin - np.einsum('ij,ji',dm,mol.intor_symmetric('int1e_kin'))
        tot_a1enuc = atom_1enuc.sum()
        a1enuc_err = tot_a1enuc - np.einsum('ij,ji',dm,moleintor.getints("int1e_nuc_sph",mol._atm,mol._bas, mol._env))
        logger.mlog(eda.stdout,"kin_err: ",kin_err)
        logger.mlog(eda.stdout,"1enuc_err: ",a1enuc_err)
        
    return atm_e1

def get_Enuc(eda):
    mol = eda.mol
    charges = mol.atom_charges()
    #coords = mol.atom_coords()
    rr = inter_distance(mol)
    #print(rr)
    rr[np.diag_indices_from(rr)] = 1e200
    #if CHECK_GEOM and numpy.any(rr < 1e-5):
    #    for atm_idx in numpy.argwhere(rr<1e-5):
    #    logger.warn(mol, 'Atoms %s have the same coordinates', atm_idx)
    #    raise RuntimeError('Ill geometry')
    atm_enucnuc = np.einsum('i,ij,j->i', charges, 1./rr, charges) * .5 
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout,"atom_enucnuc=",atm_enucnuc)
    if eda.verbose > 5:
        tot_enucnuc = atm_enucnuc.sum()
        enucnuc_err = tot_enucnuc - np.einsum('i,ij,j', charges, 1./rr, charges) * .5
        logger.mlog(eda.stdout, "err_enucnuc", enucnuc_err)
    return atm_enucnuc

def p2f(atm2bas_p):
    # TBD
    return atm2bas_p

def get_Ejk(eda, atm2bas_p, jk='jk', jktype='bas-eq'):
    '''
    jk: j -> only calc E_j
        jk -> calc E_j, E_k
    jktype: def partition coeff for (ij|kl)
            atom-eq : find atom label a,b,c,d of i,j,k,l, and suppose atoms are equal
                      e.g. (ij|kl) -> from atom 1,2,3 -> coeff {1/3, 1/3, 1/3}
            bas-eq: each basis in (ij|kl) takes 1/4
    '''
    e_coul = []
    #vhf = []
    
    #eri_s8 = mol.intor("int2e_sph",aosym='s8')
    #eri = mol.intor("int2e_sph")   
    #print(max(atm2bas))
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
            return np.array(atom_ej), np.array(atom_ek)
        elif jk=='j':
            return np.array(atom_ej)
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
        # ---------------------------------------------------
        atm_ = mol._atm.T
        atml = np.shape(mol._atm)[0]
        bas_ = mol._bas.T
        basl = np.shape(mol._bas)[0]
        env_ = mol._env
        envl = np.shape(mol._env)[0] 
        nao = len(dm)
        nbas = mol._bas.shape[0]
        #with open(eda.output+'-eda.log','a') as f:
        #if eda.verbose > 5:
        #logger.mlog(eda.stdout,"atm_",atm_)
        #logger.mlog(eda.stdout,"bas_",bas_)
        #logger.mlog(eda.stdout,"env_",env_)
        logger.slog(eda.stdout,"nao=%d",nao)
        logger.slog(eda.stdout,"nbas=%d",nbas)

        #t1 = time.time()
        atom_ej, atom_ek = preri(atm_,atml,bas_,basl,env_,envl,nao,nbas,dm,atm2bas_f,singleitem,num1)
        #t2 = time.time()
        #print("T2-T1=",t2-t1)
        #atom_energy = np.array(atom_energy)
        ##print("Atom_energy=",atom_energy)
        #with open(eda.output+'-eda.log','a') as f:
        if jk=='jk':
            atom_ej = np.array(atom_ej)[0:mol.natm]
            atom_ek = np.array(atom_ek)[0:mol.natm]
            logger.log(eda.stdout,"Atom_ej=",atom_ej)
            logger.log(eda.stdout,"Atom_ek=",atom_ek)
            if eda.verbose > 5:
                tot_aej = atom_ej.sum()
                tot_aek = atom_ek.sum()
                ej_err = tot_aej - 0.5*np.einsum('ij,ji',dm,eda.mf.get_j(mol,dm))
                ek_err = tot_aek - 0.5*(-0.5)*np.einsum('ij,ji',dm,eda.mf.get_k(mol,dm))
                logger.mlog(eda.stdout, "err_ej,ek: ",(ej_err,ek_err))
            return atom_ej, atom_ek
        elif jk=='j':
            atom_ej = np.array(atom_ej)[0:mol.natm]
            logger.log(eda.stdout,"Atom_ej=",atom_ej)
            if eda.verbose > 5:
                tot_aej = atom_ej.sum()
                ej_err = tot_aej - np.einsum('ij,ji',dm,eda.mf.get_j(mol,dm))
                logger.mlog(eda.stdout, "err_ej: ",(ej_err))
                return atom_ej



