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

from kit import logger, gjf_kit
from dftpart import scfeda
import labc

class GFEA():
    def __init__(self):

        self.inputstyle = 'frg'
        self.verbose = 4
        self.max_memory = 4000
        self.output = "test"
        self.stdout = None
        #self.stdout = sys.stdout
        #self.stdout = mol.stdout
        self.E_GFEA = 0
        self.atom_E = None

        self.gjfname = None
        self.gjf = None
        ## atomlist style ##
        self.subsys = None
        self.atomlist_tot = None
        self.spinlist_tot = None
        self.chglist_tot = None
        self.gjflist = None  # optional
        self.fchklist = None  # optional
        self.gmslist = None  #
        ################
        ## frg style ##
        self.frg = None
        self.lso = None
        self.axyz = None
        ###############

        self.totmol = None
        self.method = None
        self.dm0 = 'pyscf'
        self.built = False
        self.showinter = False
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
            self.subsys_lso, self.num_subsys, self.frg_intot = lso_parser(self.lso)
            self.gjflist, fchklist, gmslist = get_subgjf(self.gjfname, self.num_subsys)
            logger.mlog(self.stdout, "atoms in frg (tot order)\n", self.frg_intot)
            logger.mlog(self.stdout, "subsys", self.subsys_lso)
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
            for i in range(self.num_subsys):
                logger.slog(self.stdout, "## Do EDA on subsys %d ############", i+1)

                #atomlist = []
                _cen = self.subsys_lso[i][0]
                _env = self.subsys_lso[i][1]
                subsys_i = _cen + _env
                subsys_i.sort()
                logger.mlog(self.stdout, "subsys_i", subsys_i)
                logger.mlog(self.stdout, "_cen", _cen)
                logger.mlog(self.stdout, "_env", _env)
                #frag_list = get_frags(self.frg_intot, _cen, _env)
                #for label in subsys_i:
                cen_intot, cen_insub, frag_list, molchgs = get_frags(self, lab[i], lac[i], _cen, _env)
                logger.mlog(self.stdout, "atomlist_lab", lab[i])
                logger.mlog(self.stdout, "cen_intot", cen_intot)
                logger.mlog(self.stdout, "cen_insub", cen_insub)
                if 'charge' in self.method:
                    logger.log(self.stdout, "molchgs", molchgs)
                logger.slog(self.stdout, "----------------------------------------------------------")
                logger.slog(self.stdout, "frag  layer        atm_intot               atm_insub                    selfchg")
                for f in frag_list:
                    logger.slog(self.stdout, "%d      %s     "%(f.label, f.layer), endl=False)
                    logger.mlog(self.stdout, "", f.atm_intot, endl=False)
                    logger.mlog(self.stdout, "         ", f.atm_insub, endl=False)
                    logger.log(self.stdout, "        ", f.selfchg)
                logger.slog(self.stdout, "----------------------------------------------------------")

                subeda = scfeda.EDA()
                subeda.method = self.method
                subeda.gjf = self.gjflist[i]
                if 'charge' in self.method:
                    subeda.molchgs = molchgs
                subeda.output = subeda.gjf[:-4]
                subeda.verbose = self.verbose
                subeda.showinter = self.showinter
                subeda.frag_list = frag_list
                subatm_E, subE, conv, intert = subeda.kernel()
                #logger.log(self.stdout, "subatm_E", subatm_E)
                if conv==False:
                    logger.slog(self.stdout, "Error: EDA not converged")
                logger.slog(self.stdout, "center atom energies:")
                for i in range(len(cen_intot)):
                    self.atom_E[cen_intot[i]-1] = subatm_E[cen_insub[i]-1]
                    logger.slog(self.stdout, "%d %.10f", cen_intot[i], self.atom_E[cen_intot[i]-1])
                logger.slog(self.stdout, "## END ###########################")
                
        self.E_GFEA = self.atom_E.sum()
        logger.log(self.stdout, "atom_E", self.atom_E)
        logger.slog(self.stdout, "E_GFEA = %.10f", self.E_GFEA)
        #if self.showinter:
            

        return self.E_GFEA, self.atom_E

class Frag():
    def __init__(self):
        self.label = None
        self.atm_intot = []
        self.atm_insub = []
        self.layer = None
        self.selfchg = []


def get_frags(gfea, atomlist_lab, atomlist_lac, _cen, _env):
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
    nfrag = len(frags_list)
    for label in atomlist_lac:
        caplabel = nfrag
        if label==0:
            capf = Frag()
            capf.layer = 'cap'
            capf.label = caplabel + 1
            caplabel += 1
            capf.atm_intot = [-1]
            capf.atm_insub = [atomlist_lac.index(label)+1]
            frags_list.append(capf)

    molchgs = np.zeros(len(atomlist_lab))
    for label in atomlist_lab:
        if label != 0:
            cen_intot.append(label)
            cen_insub.append(atomlist_lab.index(label)+1)
            if 'charge' in gfea.method:
                molchgs[atomlist_lab.index(label)] = gfea.chglist[label-1]
        else:
            if 'charge' in gfea.method:
                molchgs[atomlist_lab.index(label)] = 0.0
    for f in frags_list:
        if f.layer is not 'cap':
            #f.selfchg = np.zeros(len(f.atm_intot))
            for label in f.atm_intot:
                f.atm_insub.append(atomlist_lac.index(label)+1)
                if 'charge' in gfea.method:
                    f.selfchg.append(gfea.chglist[label-1])
    return cen_intot, cen_insub, frags_list, molchgs

def one2zero(alist):
    if isinstance(alist, list):
        blist = []
        for item in alist:
            blist.append(one2zero(item))
    elif isinstance(alist, int):
        blist = alist - 1
    return blist

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
            submf.xc = method[0]
            submf.grids.atom_grid = (99, 590)
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
    frg_table = []
    read = False
    for line in lines:
        if '#final capped fragments' in line:
            read = True
            continue
        if read is True:
            subsys_table.append(line)
        if 'Num of capped fragments:' in line:
            break
    read = False
    for line in lines:
        if '#fragment atoms' in line:
            read = True
            continue
        if read is True:
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
    num_subsys = int(line[1].split()[-1].strip())

    return subsys_lso, num_subsys, frg_intot

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

def subsys_EDA(inputstyle,
               submf,
               dm,
               subsys,
               atomlist=[],
               spinlist=[],
               chglist=[],
               backlabel=[],
               backlist=[],
               method=['m062x', 'cc-pvtz', 'cart', 'charge'],
               coords=[],
               charges=[],
               edatype=[],
               verbose=4):
    if verbose > 4:
        logger.slog(f,"verbose = %d", verbose)
    return EDA(inputstyle, submf, dm, subsys, atomlist, spinlist, chglist,
               backlabel, backlist, method, coords, charges, edatype, verbose)


def EDA(inputstyle,
        submf,
        dm,
        subsys,
        atomlist=[],
        spinlist=[],
        chglist=[],
        backlabel=[],
        backlist=[],
        method=['m062x', 'cc-pvtz', 'cart', 'charge'],
        coords=[],
        charges=[],
        edatype='debug',
        #rctype='gebf',
        verbose=4):
    r"""
    Kwargs:
        edatype: 'full'  -> include all terms
                 'debug' -> include 'c'-related terms,
                            and 'e/e/e', 'e/e/e/e' if available
                 'smart' -> include 'debug' terms,
                            drop zero and less accurate terms
        verbose: 4 -> normal
                 5 -> print complex-frag energy
    """
    if inputstyle == 'frg':
        _cen = subsys[0]
        _env = subsys[1]
        _subsys = subsys[0] + subsys[1]
        _subsys.sort()

        s = ""
        for i in subsys[0]:
            s += (str(i) + ' ')
        cen_str = "/".join(s.split())
        #_subsys = count_from_zero(np.array(_subsys))

        ss = ""
        for i in _subsys:
            ss += (str(i) + ' ')
        sub_str = ",".join(ss.split())
    elif inputstyle == 'atomlist':
        _subsys = subsys
        cen_str = ''
        ss = ""
        for i in subsys:
            ss += (str(i) + ' ')
        sub_str = ",".join(ss.split())

    logger.slog(f,"## EDA in subsystem constructed by frag " + sub_str + " ##\n")

    EDA_data = []
    E = 0.0
    if verbose > 4:
        logger.mlog(f,"atomlist", atomlist)
    #logger.log(f,verbose)
    if verbose > 4:
        logger.slog(f,"## frag energy in real atom part ##\n")
    method_real = [item for item in method if item is not 'charge']
    # method for real atom part, removing 'charge'
    for i in _subsys:
        E = eda2_2.get_energy(inputstyle, submf, dm, [_subsys.index(i)],
                              atomlist, spinlist, chglist, method_real)
        if i in _cen:
            EDA_data.append([E[0], [[i], [], []], _cen])
        elif i in _env:
            EDA_data.append([E[0], [[], [i], []], _cen])
        else:
            logger.slog(f,"Error: label of frag not found in _cen or _env")
        if verbose > 4:
            logger.mlog(f, "  ", E[0], i, cen_stri, sep='\t\t')

    for i in _subsys:
        for j in _subsys:
            if j > i:
                E = eda2_2.get_energy(
                    inputstyle, submf, dm,
                    [_subsys.index(i), _subsys.index(j)], atomlist, spinlist,
                    chglist, method_real)
                if i in _cen:
                    if j in _cen:
                        EDA_data.append([E[0], [[i, j], [], []], _cen])
                    else:
                        EDA_data.append([E[0], [[i], [j], []], _cen])
                else:
                    if j in _cen:
                        EDA_data.append([E[0], [[j], [i], []], _cen])
                    else:
                        EDA_data.append([E[0], [[], [i, j], []], _cen])
                if verbose > 4:
                    logger.mlog(f,
                        "  ", E[0],
                        (i, j),
                        cen_str,
                        sep='\t\t')

    def ijkl_label(ijkl, _cen, _env):
        c = []  # center
        e = []  # environment
        b = []  # background
        for item in ijkl:
            if item in _cen:
                c.append(item)
            elif item in _env:
                e.append(item)
            else:
                b.append(item)
        return [c, e, b]

    for i in _subsys:
        for j in _subsys:
            for k in _subsys:
                if j > i and k > j:
                    E = eda2_2.get_energy(
                        inputstyle, submf, dm,
                        [_subsys.index(i),
                         _subsys.index(j),
                         _subsys.index(k)], atomlist, spinlist, chglist,
                        method_real)
                    EDA_data.append(
                        [E[0], ijkl_label([i, j, k], _cen, _env), _cen])
                    if verbose > 4:
                        logger.mlog(f,
                            "  ", E[0],
                            (i, j, k),
                            cen_str,
                            sep='\t\t')

    for i in _subsys:
        for j in _subsys:
            for k in _subsys:
                for l in _subsys:
                    if j > i and k > j and l > k:
                        E = eda2_2.get_energy(inputstyle, submf, dm, [
                            _subsys.index(i),
                            _subsys.index(j),
                            _subsys.index(k),
                            _subsys.index(l)
                        ], atomlist, spinlist, chglist, method_real)
                        EDA_data.append(
                            [E[0],
                             ijkl_label([i, j, k, l], _cen, _env), _cen])
                        if verbose > 4:
                            logger.log(f,
                                "  ", E[0],
                                (i, j, k, l),
                                cen_str,
                                sep='\t\t')

    logger.slog(f,"## inter energy in real atom part ##")
    EDA_RR_data = eda2_2.data2inter(EDA_data)  # real-real interaction
    logger.slog(f,"---------------------------------------------------------")
    logger.slog(f,'        Energy       interaction     center')
    for item in EDA_RR_data:
        logger.mlog(f,"  ", item[0], item[1], item[2], sep='\t\t')
    logger.slog(f,"---------------------------------------------------------\n\n")

    if 'charge' in method:
        logger.log(f,"coords", coords)
        logger.slog(f,"## iter energy between real atom frags and charges ##")
        logger.slog(f,"---------------------------------------------------------")
        logger.slog(f,'        Energy       interaction     center')
        EDA_RC_data = []  # real atom-charge interaction
        logger.mlog(f, "backlabel", backlabel)
        for i in _subsys:
            if edatype == 'subsys':
                if i not in _cen:
                    continue
            for j in backlabel:
                backrange = backlist[backlabel.index(j)]
                logger.mlog(f, "atomlist/backrange of subsys", atomlist[_subsys.index(i)], backrange)
                E = eda2_2.RC_inter(inputstyle, submf, dm, [_subsys.index(i)],
                                    atomlist, spinlist, chglist, 'qmmm',
                                    coords[np.ix_(backrange)],
                                    charges[np.ix_(backrange)])
                EDA_RC_data.append([E, ijkl_label([i, j], _cen, _env), _cen])
                logger.mlog(f,
                    "  ", E, (i, j), "*", cen_str, sep='\t\t')

        for i in _subsys:
            for j in _subsys:
                if j > i:
                    if edatype == 'subsys':
                        if i not in _cen and j not in _cen:
                            continue
                    for k in backlabel:
                        backrange = backlist[backlabel.index(k)]
                        E = eda2_2.RC_inter(
                            inputstyle, submf, dm,
                            [_subsys.index(i),
                             _subsys.index(j)], atomlist, spinlist, chglist,
                            'qmmm', coords[np.ix_(backrange)],
                            charges[np.ix_(backrange)])
                        EDA_RC_data.append(
                            [E, ijkl_label([i, j, k], _cen, _env), _cen])
                        logger.mlog(f,
                            "  ", E,
                            (i, j, k), "*",
                            cen_str,
                            sep='\t\t')
        for i in _subsys:
            if edatype == 'subsys':
                if i not in _cen:
                    continue
            for j in backlabel:
                for k in backlabel:
                    if k > j:
                        #backrange = backlist[backlabel.index(j)] \
                        #            + backlist[backlabel.index(k)]
                        #E = eda.RC_inter(inputstyle, submf, dm,
                        #                 [_subsys.index(i)], atomlist,
                        #                 spinlist, chglist, 'qmmm',
                        #                 coords[np.ix_(backrange)],
                        #                 charges[np.ix_(backrange)])
                        E = 0.0
                        EDA_RC_data.append(
                            [E, ijkl_label([i, j, k], _cen, _env), _cen])
                        logger.mlog(f,
                            "  ", E,
                            (i, j, k), "*",
                            cen_str,
                            sep='\t\t')

        for i in _subsys:
            for j in _subsys:
                for k in _subsys:
                    if edatype == 'subsys':
                        if i not in _cen and j not in _cen and k not in _cen:
                            continue
                    if j > i and k > j:
                        for l in backlabel:
                            #backrange = backlist[backlabel.index(l)]
                            #E = eda.RC_inter(inputstyle, submf, dm, [
                            #    _subsys.index(i),
                            #    _subsys.index(j),
                            #    _subsys.index(k)], atomlist, spinlist, chglist,
                            #                 'qmmm', coords[np.ix_(backrange)],
                            #                 charges[np.ix_(backrange)])
                            E = 0.0
                            EDA_RC_data.append([
                                E,
                                ijkl_label([i, j, k, l], _cen, _env), _cen
                            ])
                            logger.mlog(f,
                                "  ", E,
                                (i, j, k, l), "*",
                                cen_str,
                                sep='\t\t')
        for i in _subsys:
            for j in _subsys:
                if j > i:
                    if edatype == 'subsys':
                        if i not in _cen and j not in _cen:
                            continue
                    for k in backlabel:
                        for l in backlabel:
                            if l > k:
                                #backrange = backlist[backlabel.index(k)] \
                                #            + backlist[backlabel.index(l)]
                                #E = eda.RC_inter(
                                #    inputstyle, submf, dm,
                                #    [_subsys.index(i),
                                #     _subsys.index(j)], atomlist, spinlist,
                                #    chglist, 'qmmm', coords[np.ix_(backrange)],
                                #    charges[np.ix_(backrange)])
                                E = 0.0
                                EDA_RC_data.append([
                                    E,
                                    ijkl_label([i, j, k, l], _cen, _env), _cen
                                ])
                                logger.mlog(f,
                                    "  ", E,
                                    (i, j, k, l), "*",
                                    cen_str,
                                    sep='\t\t')
        for i in _subsys:
            if edatype == 'subsys':
                if i not in _cen:
                    continue
            for j in backlabel:
                for k in backlabel:
                    for l in backlabel:
                        if k > j and l > k:
                            #backrange = backlist[backlabel.index(j)] \
                            #            + backlist[backlabel.index(k)] \
                            #            + backlist[backlabel.index(l)]
                            #E = eda.RC_inter(inputstyle, submf, dm,
                            #                 [_subsys.index(i)], atomlist,
                            #                 spinlist, chglist, 'qmmm',
                            #                 coords[np.ix_(backrange)],
                            #                 charges[np.ix_(backrange)])
                            E = 0.0
                            EDA_RC_data.append([
                                E,
                                ijkl_label([i, j, k, l], _cen, _env), _cen
                            ])
                            logger.mlog(f,
                                "  ", E,
                                (i, j, k, l), "*",
                                cen_str,
                                sep='\t\t')
        logger.slog(f,"---------------------------------------------------------\n\n")
    if 'charge' in method:
        return EDA_RR_data + EDA_RC_data
    else:
        return EDA_RR_data


def NN_input(inputstyle,
             submf,
             dm,
             subsys,
             atomlist=[],
             spinlist=[],
             chglist=[],
             method=[
                 'm062x',
                 'cc-pvtz',
                 'cart',
             ]):
    EDA_RR_data = EDA(
        inputstyle,
        submf,
        dm,
        subsys,
        atomlist,
        spinlist,
        chglist,
        method=method)
    sumE = 0.0
    for item in EDA_RR_data:
        sumE += item[0]
    logger.mlog(f,"## tot_EDA_sum for subsys ", subsys, " = %18.10f ##" % sumE)
    EDA_NN = []
    for item in subsys:
        E = 0.0
        for jtem in EDA_RR_data:
                if len(jlabel) == 2:
                    E += jtem[0] * 0.5
                if len(jlabel) == 3:
                    E += jtem[0] * (1 / 3)
                if len(jlabel) == 4:
                    E += jtem[0] * 0.25
        EDA_NN.append([E, item])
    return EDA_NN
