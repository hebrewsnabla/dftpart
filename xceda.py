import numint_sep
import gen_grid_sep
from pyscf import dft
import numpy as np
import time
import scfeda
from QCKit import logger

ok = ['FAIL','OK']
TOT_THRESH = 1e-9

def analysis(ks, kwd=['atom','dexc']):
    '''
    kwd: atom   output every atomic xc E
         dexc   output exc_tot - sum(atom_exc)
    '''
    dm = ks.make_rdm1()
    ni = numint_sep.NumInt()
    mol = ks.mol
    grids = gen_grid_sep.Grids(ks.mol)
    tot_exc = dft.numint.nr_rks(ni, mol, ks.grids, ks.xc, dm, max_memory=mol.max_memory)[1]
    atom_exc = np.hstack(numint_sep.nr_rks_sep(ni, mol, grids, ks.xc, dm)[3])
    atom_exc_tot = atom_exc.sum()
    dexc = atom_exc_tot - tot_exc
    
    out = open(ks.mol.output[:-10]+'-xceda.log','a')
    if 'atom' in kwd:
        for i in range(ks.mol.natm):
            out.write("%s  &  %.12f  \\\\ \n" % (ks.mol.atom_symbol(i), atom_exc[i])) 
    if 'dexc' in kwd:
        exc_conv = int(abs(dexc) < TOT_THRESH)
        out.write("E_xc err: %g  %s\n" % (dexc, ok[exc_conv]))
    out.close()
    return atom_exc, dexc

def get_atmexc(eda,atm2bas):
    t1 = time.time()
    dm = eda.dm
    ks = eda.mf
    ni = ks._numint
    mol = ks.mol
    grids = gen_grid_sep.Grids(ks.mol)
    atom_exc = numint_sep.nr_rks_sep(ni, mol, grids, ks.xc, dm)[3]
    #with open(eda.output+'-eda.log','a') as f:
    logger.log(eda.stdout, "Atom_exc(pure):", atom_exc)
    if eda.verbose > 5:
        tot_aexc = atom_exc.sum()
        err_exc = tot_aexc - dft.numint.nr_rks(ni,mol,ks.grids,ks.xc,dm)[1]
        logger.mlog(eda.stdout,"err_exc",err_exc)
    t2 = time.time()
    #with open(eda.output+'-eda.log','a') as f:
    logger.slog(eda.stdout, "time for Exc: %.5f\n", (t2-t1))
    omega, alpha, hyb = ni.rsh_and_hybrid_coeff(ks.xc, mol.spin)
    if ks.omega is not None: omega = ks.omega
    if abs(hyb) > 1e-10:
        atom_ej, atom_ek = scfeda.get_Ejk(eda,atm2bas,'jk')
        atom_exc += hyb * atom_ek
        #if abs(omega) > 1e-10:
        #    eklr = get_eklr(mol, dm, omega, hermi=1)*(alpha-hyb)
        #    atom_exc += eklr
        #with open(eda.output+'-eda.log','a') as f:
        logger.log(eda.stdout, "Atom_exc(hydrid):", atom_exc)
        t3 = time.time()
        #with open(eda.output+'-eda.log','a') as f:
        logger.slog(eda.stdout, "time for Ej, Ek: %.5f\n", (t3-t2))
    else:
        atom_ej = scfeda_get_Ejk(eda,atm2bas,'j')
        t3 = time.time()
        #with open(eda.output+'-eda.log','a') as f:
        logger.slog(eda.stdout, "time for Ej: %.5f\n", (t3-t2))
    return atom_exc, atom_ej
"""
def get_eklr(mol, dm, omega, hermi):
    nao = dm.shape[-1]
    with mol.with_range_coulomb(omega):
        intor = mol._add_suffix('int2e')
        eklr = 
    return eklr
"""
