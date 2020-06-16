"""
GFEA, Generalized Fragment Energy Assembler

v3  Jul 2019

Shirong Wang <wangshirongnju@gmail.com>

References:
  FEA -- J. Chem. Phys. 124, 154102 (2006)  doi.org/10.1063/1.2186997

"""

from pyscf import gto, scf, dft, symm, qmmm, lib
#from eda import eda2_2, getqq
#from gau2py import gau2py
import numpy as np
import sys, os
import time
import subprocess

from kit import logger, gjf_kit, misc
from dftpart import scfeda
import labc

class GFEA():
    def __init__(self):
        # general settings
        self.verbose = 4
        self.max_memory = 4000
        self.output = "test"
        self.stdout = None
        # GFEA input (frg style)
        self.inputstyle = 'frg'
        self.method = []
        self.gjfname = None
        self.gjf = None # auto detect
        self.frg = None # auto detect
        self.lso = None # auto detect
        self.cha = None # auto detect
        self.labc = None # auto detect
        self.axyz = None
        self.dm0 = 'pyscf' 
        ## atomlist style (deprecated) ##
        #self.subsys = None
        self.atomlist_tot = None
        self.spinlist_tot = None
        self.chglist_tot = None

        # GFEA variables
        self.E_GFEA = 0
        self.atom_E = None
        self.totmol = None
        self.method = None
        self.built = False
        self.subsys = None
        self.gjflist = None  
        #self.fchklist = None  
        #self.gmslist = None  #
        # GFEA settings
        self.showinter = False
        self.do_deriv = False
        self.exclude_cap = True
        #self.Gau_version = 'G16'
        #self.fchklist = None
    def build(self):
        os.system("echo '' > "+self.output+'-gfea.log')
        self.stdout = open(self.output+'-gfea.log','a')
        if self.inputstyle == 'atomlist':
            inp = self.subsys, self.atomlist_tot, self.spinlist_tot, self.chglist_tot, self.gjfname
        elif self.inputstyle == 'frg':
            if self.frg is None:
                self.frg = self.gjfname + '/' + self.gjfname + '.frg'
            if self.lso is None:
                self.lso = self.gjfname + '/' + self.gjfname + '.lso'
            if self.axyz is None:
                self.axyz = self.gjfname + '/' + self.gjfname + '.axyz'
            self.cha = self.gjfname + '/' + self.gjfname + '.cha'
            self.gjf = self.gjfname + '.gjf'
            self.labc = self.gjfname + '/' + self.gjfname + '.labc'
        else:
            logger.slog(f,"Invalid Input Style!")
        #E_GFEA, atom_E =  kernel(self)
        #f.close()
        totmol = gto.Mole()
        totmol.atom, coo, cha, totmol.charge, totmol.spin = gjf_kit.gjf_parser(self.gjf)
        totmol.build()
        self.totmol = totmol
        self.atom_E = np.zeros(totmol.natm)

        if self.inputstyle == 'atomlist':
            subsys, atomlist_tot, spinlist_tot, chglist_tot, gjfname = inp
            num_subsys = len(subsys)
            num_frag = len(atomlist_tot)
            gjflist, fchklist, gmslist = get_subgjf(gjfname, num_subsys)
        elif self.inputstyle == 'frg':
            #frg, lso, axyz, gjfname = self.frg, self.lso, self.axyz, self.gjfname
            #num_subsys = get_num_subsys(lso)
            self.subsys_lso, self.num_subsys_prm, self.num_subsys_tot, self.frg_intot, self.num_atoms = lso_parser(self.lso)
            #print(self.subsys_lso)
            if self.do_deriv:
                self.num_subsys = self.num_subsys_tot
            else:
                self.num_subsys = self.num_subsys_prm
            self.gjflist, fchklist, gmslist = get_subgjf(self.gjfname, self.num_subsys)
            logger.mlog(self.stdout, "atoms in frg (tot order)\n", self.frg_intot)
            logger.mlog(self.stdout, "subsys", self.subsys_lso)
            self.num_frag = len(self.frg_intot)
            logger.slog(self.stdout, "number of frags: %d", self.num_frag)
            logger.slog(self.stdout, "number of subsystems: %d", self.num_subsys)
        if self.dm0 == 'Gaussian':
            run_Gau(gjflist)
            # to be done
            #
        elif self.dm0 == 'fchk':
            readpunch = ('punchbasis' in method)
            dmlist, submflist, coords, charges = fchk2dm(
                gjflist, fchklist, gmslist, verbose, max_memory, method, readpunch)

        built = True

    def kernel(self):
        if self.built == False: self.build()

        if self.inputstyle == 'atomlist':
            for i in len(submflist):
                submf = submflist[i]
                dm = dmlist[i]
                atomlist = []
                spinlist = []
                chglist = []
                backlist = []
                subsys_i = subsys_lso[i][0] + subsys_lso[i][1]
                subsys_i.sort()
                # combine cen and env
                for label in subsys_i:
                    atomlist.append(atomlist_tot[label - 1])  # list of num_atom
                    # in real atom frags
                    spinlist.append(spinlist_tot[label - 1])
                    chglist.append(chglist_tot[label - 1])
                backlabel = [
                    item for item in range(1,
                                           len(atomlist_tot) + 1)
                    if item not in subsys_i
                ]
                for label in backlabel:
                    backlist.append(atomlist_tot[label - 1])  # list of num_chg
                    # in background charge frags

                EDA_data_tot += subsys_EDA('atomlist', submf, dm, subsys[i],
                                           atomlist, spinlist, chglist, backlabel,
                                           backlist, method, coords, charges,
                                           'subsys', verbose)

        elif self.inputstyle == 'frg':
            #atomlist_tot_frg, spinlist_tot_frg, chglist_tot_frg = frg_parser(
            #    self.frg)
            #logger.mlog(self.stdout, "atomlist_tot", atomlist_tot_frg)
            #num_frag = len(atomlist_tot_frg)
            #logger.mlog(self.stdout, "num_frag", num_frag)
            lab = labc.labc_parser(self.labc, "lab")
            lac = labc.labc_parser(self.labc, "lac")
            if ('charge' in self.method) or ('qmmm' in self.method):
                self.chglist = cha_parser(self.cha)
            all_intert = []
            for i in range(self.num_subsys):
                logger.slog(self.stdout, "## Do EDA on subsys %d ############", i+1)

                #atomlist = []
                _cen = self.subsys_lso[i][0]
                _env = self.subsys_lso[i][1]
                _bg = get_bglabel(self.num_frag, _cen, _env)
                subsys_i = _cen + _env
                subsys_i.sort()
                logger.mlog(self.stdout, "subsys_i", subsys_i)
                logger.mlog(self.stdout, "_cen", _cen)
                logger.mlog(self.stdout, "_env", _env)
                logger.mlog(self.stdout, "_bg", _bg)
                #frag_list = get_frags(self.frg_intot, _cen, _env)
                #for label in subsys_i:
                try:
                    lab_i = lab[i]
                except:
                    lab_i = [0]*len(lac[i])
                cen_intot, cen_insub, frag_list, molchgs = get_frags(self, lab_i, lac[i], _cen, _env, _bg)
                try:
                    logger.mlog(self.stdout, "atomlist_lab", lab[i])
                    logger.mlog(self.stdout, "cen_intot", cen_intot)
                    logger.mlog(self.stdout, "cen_insub", cen_insub)
                except:
                    pass
                if 'charge' in self.method:
                    logger.log(self.stdout, "molchgs", molchgs)
                logger.slog(self.stdout, "----------------------------------------------------------")
                logger.slog(self.stdout, "frag  layer        atm_intot               atm_insub                    selfchg")
                for f in frag_list:
                    logger.slog(self.stdout, "%d      %s     "%(f.label, f.layer), endl=False)
                    logger.mlog(self.stdout, "", f.atm_intot, endl=False)
                    logger.mlog(self.stdout, "         ", f.atm_insub, endl=False)
                    logger.log(self.stdout, "        ", f.selfchg, digits=4, newl=False)
                logger.slog(self.stdout, "----------------------------------------------------------")

                subeda = scfeda.EDA()
                subeda.method = self.method
                subeda.gjf = self.gjflist[i]
                if 'charge' in self.method:
                    subeda.molchgs = molchgs
                subeda.output = subeda.gjf[:-4]
                subeda.verbose = self.verbose
                subeda.showinter = self.showinter
                subeda.exclude_cap = self.exclude_cap
                subeda.frag_list = frag_list
                subeda.totnum_frag = self.num_frag
                subatm_E, subE, conv, intert = subeda.kernel()
                all_intert.append(intert)
                #logger.log(self.stdout, "subatm_E", subatm_E)
                if conv==False:
                    logger.slog(self.stdout, "Error: EDA not converged")
                logger.slog(self.stdout, "center atom energies:")
                if i > (self.num_subsys_prm-1):
                    continue # skip deriv subsys
                for i in range(len(cen_intot)):
                    self.atom_E[cen_intot[i]-1] = subatm_E[cen_insub[i]-1]
                    logger.slog(self.stdout, "%d %.10f", cen_intot[i], self.atom_E[cen_intot[i]-1])
                logger.slog(self.stdout, "## END ###########################")
                
        self.E_GFEA = self.atom_E.sum()
        logger.log(self.stdout, "atom_E", self.atom_E)
        logger.slog(self.stdout, "E_GFEA (atom) = %.10f", self.E_GFEA)
        if self.showinter:
            #pass
            tot_inter1, tot_inter2, tot_inter3, tot_inter4 = alloc_inter(all_intert)   
            logger.ilog(self.stdout, 'inter1', tot_inter1) 
            logger.ilog(self.stdout, 'inter2', tot_inter2) 
            logger.ilog(self.stdout, 'inter3', tot_inter3) 
            logger.ilog(self.stdout, 'inter4', tot_inter4) 
            self.tot_inter = [tot_inter1, tot_inter2, tot_inter3, tot_inter4]

        return self.E_GFEA, self.atom_E, self.tot_inter

def alloc_inter(intert_list):
    #import priority
    inter1 = misc.EDict()
    inter2 = misc.EDict()
    inter3 = misc.EDict()
    inter4 = misc.EDict()

    for sub in intert_list:
        RR1,RR2,RR3,RR4,RC2,RC3 = sub[0],sub[1],sub[2],sub[3],sub[4],sub[5]
        inter1.update(RR1)
        inter2.update(RR2, RC2)
        inter3.update(RR3, RC3)
        inter4.update(RR4)
    return inter1, inter2, inter3, inter4




class Frag():
    def __init__(self):
        self.label = None
        self.atm_intot = []
        self.atm_insub = []
        self.layer = None
        self.selfchg = []


def get_frags(gfea, atomlist_lab, atomlist_lac, _cen, _env, _bg):
    #atomlist_lab = lab[i]
    #atomlist_lac = lac[i]
    cen_intot = []
    cen_insub = []
    env_intot = []
    env_insub = []
    frags_list = []
    for cenfrg in _cen:
        f = Frag()
        f.atm_intot = gfea.frg_intot[cenfrg-1]
        f.label = cenfrg
        f.layer = 'c'
        frags_list.append(f)
    for envfrg in _env:
        f = Frag()
        f.atm_intot = gfea.frg_intot[envfrg-1]
        f.label = envfrg
        f.layer = 'e'
        frags_list.append(f)
    for bgfrg in _bg:
        f = Frag()
        f.atm_intot = gfea.frg_intot[bgfrg-1]
        f.label = bgfrg
        f.layer = 'q'
        frags_list.append(f)
    
    #if len(_cen) > 0:
    #    maxfrag = max([max(_cen),max(_env)])
    #else:
    #    maxfrag = max(_env)
    maxfrag = gfea.num_frag

    caplabel = maxfrag
    for l in range(len(atomlist_lac)):
        label = atomlist_lac[l]
        if label==0:
            capf = Frag()
            capf.layer = 'cap'
            capf.label = caplabel + 1
            caplabel += 1
            capf.atm_intot = [-1]
            capf.atm_insub = [l+1]
            frags_list.append(capf)

    molchgs = np.zeros(len(atomlist_lab))
    for label in atomlist_lab:
        if label != 0:
            cen_intot.append(label)
            cen_insub.append(atomlist_lab.index(label)+1)
            if 'charge' in gfea.method:
                molchgs[atomlist_lab.index(label)] = gfea.chglist[label-1]
        #else:
        #    if 'charge' in gfea.method:
        #        molchgs[atomlist_lab.index(label)] = 0.0
    list_chg = get_listchg(gfea.num_atoms, atomlist_lac)
    #print(gfea.num_atoms)
    #print(atomlist_lac)
    #print(list_chg)
    for f in frags_list:
        if (f.layer == 'c') or (f.layer == 'e'):
            #f.selfchg = np.zeros(len(f.atm_intot))
            for label in f.atm_intot:
                f.atm_insub.append(atomlist_lac.index(label)+1)
                if 'charge' in gfea.method:
                    f.selfchg.append(gfea.chglist[label-1])
        elif f.layer == 'q':
            for label in f.atm_intot:
                f.atm_insub.append(list_chg.index(label)+1)
                if 'charge' in gfea.method:
                    f.selfchg.append(gfea.chglist[label-1])


    return cen_intot, cen_insub, frags_list, molchgs

def get_bglabel(num_frag, _cen, _env):
    _bg = []
    for i in range(1,num_frag+1):
        if i in _cen:
            continue
        if i in _env:
            continue
        _bg.append(i)
    _bg.sort()
    return _bg

def get_listchg(num_atoms, lac):
    list_chg = []
    for i in range(1,num_atoms+1):
        if i not in lac:
            list_chg.append(i)
    return list_chg


def cha_parser(cha):
    with open(cha,'r') as f:
        data = f.readlines()
    chglist = []
    for line in data:
        line = line.strip().split()
        chglist.append(float(line[2]))
    return np.array(chglist)

def get_subgjf(gjfname, num_subsys):
    gjflist = []
    fchklist = []
    gmslist = []
    for i in range(num_subsys):
        gjflist.append(gjfname + '_subsys/' + gjfname + '_' + str(i + 1) + '.gjf')
        fchklist.append(gjfname + '_' + str(i + 1) + '.fchk')
        gmslist.append(gjfname + '_' + str(i + 1) + '.7')
    return gjflist, fchklist, gmslist


def run_Gau_cmd(gjf):
    os.system(Gau_version + ' ' + gjf)


def run_Gau(gjflist):
    dmlist = []
    for gjf in gjflist:
        run_Gau_cmd(gjf)


def fchk2dm(gjflist,
            fchklist,
            gmslist,
            verbose,
            max_memory,
            method,
            readpunch=False):
    dmlist = []
    submflist = []
    if verbose > 4:
        logger.log(f,"read fchk: ", fchklist)
    for fchk in fchklist:
        gms = gmslist[fchklist.index(fchk)]
        gjf = gjflist[fchklist.index(fchk)]
        if readpunch == True:
            gau2py.get_geombasis(gms, 'submol')
            nwcfile = gms[:-1] + 'nwc'
            with open(nwcfile, 'r') as f:
                #    #logger.log(f,nwcfile)
                ff = f.read()
            submol = gto.Mole()
            exec(ff, {'submol': submol, 'gto': gto})

            #logger.log(f,globals())
        # will generate submol =  gto.Mole() and get submol.atom, submol.basis automatically
        atom, coords, charges, molcharge, spin = gjf_parser(gjf)
        if readpunch == False:
            submol = gto.Mole()
            #submol.atom = atom
            submol.basis = method[1]
        submol.atom = atom
        submol.charge = molcharge
        submol.spin = spin
        submol.verbose = verbose

        submol.max_memory = max_memory
        submol.cart = ('cart' in method)
        submol.output = gjf[:-4] + '-pyscf.log'
        submol.build()

        if method[0] == 'hf':
            submf = scf.RHF(submol)
        elif is_dft(method[0]):
            submf = dft.RKS(submol)
            submf.xc = method[0]
            submf.grids.atom_grid = (99, 590)
        submf = scf.addons.remove_linear_dep_(submf, threshold=1e-6, lindep=1e-7)
        #if 'charge' in method:
        #    #h1e_hf = submf.get_hcore()
        #    #charges, coords = gjf_parser(gjf)
        #    submf = qmmm.mm_charge(submf, coords, charges, unit='a')

        dm = gau2py.get_dm(submf, fchk)
        dmlist.append(dm)
        submflist.append(submf)
        #logger.log(f,dm.shape, submf.get_ovlp().shape)
    return dmlist, submflist, coords, charges




def frg_parser(frg):
    with open(frg, 'r') as f:
        lines = f.readlines()
    atomlist_tot_frg = []
    spinlist_tot_frg = []
    chglist_tot_frg = []
    for line in lines:
        line = line.split()
        if len(line) != 4:
            continue
        #logger.log(f,line)
        atom = line[2].strip().strip('(').strip(')').split(',')
        atomlist_frg = [int(item) for item in atom]
        atomlist_tot_frg.append(atomlist_frg)
        spinlist_tot_frg.append(
            int(line[1].strip()) - 1)  # use pyscf spin = 2S, not 2S+1
        chg = line[3].strip()
        if chg == '+': chgvalue = 1
        elif '+' in chg: chgvalue = int(chg.strip('+'))
        elif '-' in chg: chgvalue = -int(chg.strip('-'))
        else: chgvalue = int(chg)
        chglist_tot_frg.append(chgvalue)
    return atomlist_tot_frg, spinlist_tot_frg, chglist_tot_frg


def lso_parser(lso, style='by_sub'):
    r''' subsys info in .lso file like
         -------------------
         1-3    4,5
         4      2,3,5,6
         ...
         -------------------

         would be parsed into

         [[[1],[2,3,4,5]],
          [[2],[1,3,4,5]],
          [[3],[1,2,4,5]],
          [[4],[2,3,5,6]],
          ...]
         if style='by_cen' (i.e. every center frag gets a line)

         or
         [[[1,2,3],[4,5]],
          [[4],[2,3,5,6]],
          ...]
         if style='by_sub' (i.e. every subsystem calculated gets a line)

    '''
    with open(lso, 'r') as f:
        lines = f.readlines()
    subsys_table = []
    subsys_table_deriv = []
    frg_table = []

    read = False
    for line in lines:
        if '#final capped fragments' in line:
            read = True
            continue
        if read:
            subsys_table.append(line)
        if 'Num of capped fragments:' in line:
            break
    
    read = False
    for line in lines:
        if 'Deriv' in line:
            read = True
            continue
        if read:
            subsys_table_deriv.append(line)
            if '===' in line:
                break
    
    read = False
    for line in lines:
        if '#fragment atoms' in line:
            read = True
            continue
        if read:
            frg_table.append(line)
            if '#end' in line:
                break
            
    subsys_data = subsys_table[3:-2]
    subsys_lso = []
    for line in subsys_data:
        #logger.log(f,line)
        line = line.split()
        cen = line[0].strip()
        env = line[1].strip()
        cen_list = comma_parser(cen)
        env_list = comma_parser(env)
        if style == 'by_sub':
            subsys_lso.append([cen_list, env_list])
        elif style == 'by_cen':
            for i in cen_list:
                #logger.log(f,cen_list)
                cen_rest = [item for item in cen_list if item is not i]
                real_env = env_list + cen_rest
                real_env.sort()
                subsys_lso.append([[i], real_env])
    for line in subsys_table_deriv[:-1]:
        line = line.split()
        #cen = line[0].strip()
        env = line[4].strip()
        cen_list = []
        env_list = comma_parser(env)
        #if style == 'by_sub':
        subsys_lso.append([cen_list, env_list])
        #elif style == 'by_cen':
        #    for i in cen_list:
        #        cen_rest = [item for item in cen_list if item is not i]
        #        real_env = env_list + cen_rest
        #        real_env.sort()
        #        subsys_lso.append([[i], real_env])

    #frg_data = frg_table[3:-2]
    #print(frg_table)
    frg_intot = []
    for line in frg_table[3:-2]:
        line = line.strip().split()
        frg_i = line[5]
        if len(line)==7: frg_i = [line[5],line[6]]
        #print(frg_i)
        frg_i = comma_parser(frg_i)
        frg_intot.append(frg_i)

    line = subprocess.getstatusoutput("grep 'Num of capped fragments' " + lso)
    #logger.log(f,line[1].split())
    num_subsys_prm = int(line[1].split()[-1].strip())
    line2 = subprocess.getstatusoutput("grep 'Number of subsystems:' " + lso)
    num_subsys_tot = int(line2[1].split()[-1].strip())
    line3 = subprocess.getstatusoutput("grep 'Number of atoms:' " + lso)
    num_atoms = int(line3[1].split()[-1].strip())
    

    return subsys_lso, num_subsys_prm, num_subsys_tot, frg_intot, num_atoms

def comma_parser(exprs):
    if isinstance(exprs, list):
        out = []
        for item in exprs:
            out += comma_parser(item)
        return out
    else:
        exprs = exprs.strip('(').strip(')')
        exprs = exprs.split(',')
        expr_list = []
        for item in exprs:
            if '-' in item:
                item = item.split('-')
                item_head = int(item[0])
                item_tail = int(item[1])
                item_list = list(range(item_head, item_tail + 1))
            else:
                item_list = [int(item)]
            expr_list += item_list
    return expr_list


def axyz_parser(axyz, subsyslabel, _atomlist):
    with open(axyz, 'r') as f:
        axyzdata = f.readlines()
    i = 0
    j = -1
    for line in axyzdata:
        # determine the range of subsys we want
        j += 1
        if line.strip() == "":
            i += 1
        if i == subsyslabel + 1:
            startline = j + 3
            #logger.log(f,startline)
            break
    i = 0
    j = -1
    for line in axyzdata:
        j += 1
        if line.strip() == "":
            i += 1
        if i == subsyslabel + 2:
            endline = j - 1
            break

    subsyslines = axyzdata[startline:endline + 1]
    #for i in subsyslines:
    #    logger.log(f,i)

    sequence = []
    cup = 0
    for line in subsyslines:
        line = line.split()
        if '*' not in line[-1]:
            sequence.append(int(line[-1].strip()))
        else:
            cup += 1
    atomlist_frg = []
    for i in _atomlist:
        for j in sequence:
            if i == j:
                atomlist_frg.append(sequence.index(j))
                matched = True
        if not matched:
            logger.slog(f,"Cannot find %dth atom in axyz\n", i)
    #logger.log(f,subsyslines[0])
    #logger.log(f,_atomlist)
    #logger.log(f,sequence)
    #logger.log(f,atomlist_frg)
    return atomlist_frg, cup



#def build(gfea):


def ctr(label):
    return set(label[0] + label[1] + label[2])




def gfea_anal(EDA_data_tot, num_frag, qq_data, gfea_type='smart', verbose=5):
    if verbose > 4:
        logger.slog(f,"## c-c interaction (Q-Q) ##")
        logger.slog(f,"----------------------------------")
        for item in qq_data:
            logger.mlog(f, "  ", item, qq_data[item])
        logger.slog(f,"----------------------------------")
    logger.slog(f,"## collect GFEA data ##")
    logger.slog(f,"-----------------------------------")
    GFEA_data = []
    E_GFEA = 0.0
    cmd = get_cmd_nat(num_frag)
    for item in cmd:
        bestprior = 20
        candidate = []
        candilabel = []
        for jtem in EDA_data_tot:
            E = jtem[0]
            label = jtem[1]
            if ctr(label) == set(item):
                if priority(label) < bestprior:
                    candidate = [E]
                    candilabel = [label]
                    bestprior = priority(label)
                elif priority(label) == bestprior:
                    candidate.append(E)
                    candilabel.append(label)
        if len(candidate) == 0:
            E_cand = 0.0
        else:
            E_cand = np.mean(candidate)

        if is_cb(candilabel):
            E_cand *= 2
            E_cand -= qq_data[(candilabel[0][0][0],candilabel[0][2][0])]
        GFEA_data.append([E_cand, item])
        if verbose > 4:
            logger.mlog(f,"  ", E_cand, item, candidate, candilabel)
        else:
            logger.mlog(f,"  ", E_cand, item)
        E_GFEA += E_cand

    logger.slog(f,"-----------------------------------")

    logger.slog(f,"\n## GFEA Energy = %18.10f ##", E_GFEA)

    logger.slog(f,"## decomposit energy on frags ##")
    logger.slog(f,"-----------------------------------")
    frag_E = []
    for i in range(1, num_frag + 1):
        eps = 0.0
        for item in GFEA_data:
            if i in item[1]:
                eps += item[0] / len(item[1])
        frag_E.append([eps, i])
        logger.mlog(f,"  ", eps, '   ', i)
    logger.slog(f,"-----------------------------------")

    return E_GFEA, GFEA_data, frag_E



def has_ovlp(a, b):
    if len(list(set(a) & set(b))) != 0:
        return True
    else:
        return False
