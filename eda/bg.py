from pyscf import gto, lib, df
import numpy as np

def RC_inter(inputstyle, submf, dm_tot, fraglabels=[], atomlist=[], spinlist=[],
             chglist=[], method='qmmm', coords=[], charges=[]):
    t1 = datetime.datetime.now()
    mol = submf.mol
    frag, E_nuc, basis_range = build_frag(inputstyle, mol, fraglabels,
                                          atomlist, spinlist, chglist)
    frag.cart = mol.cart
    dm = dm_tot[np.ix_(basis_range, basis_range)]
    if len(fraglabels) == 1:
        inter_RC = inter_elecbg(frag, dm, coords, charges) + inter_nucbg(
            frag, coords, charges)
    elif len(fraglabels) == 2:
        frag1, E_nuc1, basis_range1 = build_frag(
            inputstyle, mol, [fraglabels[0]], atomlist, spinlist, chglist)
        frag2, E_nuc2, basis_range2 = build_frag(
            inputstyle, mol, [fraglabels[1]], atomlist, spinlist, chglist)
        frag1.cart = mol.cart
        frag2.cart = mol.cart
        dm1 = dm_tot[np.ix_(basis_range1, basis_range1)]
        dm2 = dm_tot[np.ix_(basis_range2, basis_range2)]
        inter_RC = inter_elecbg(frag, dm, coords, charges) - inter_elecbg(
            frag1, dm1, coords, charges) - inter_elecbg(
                frag2, dm2, coords, charges)
    t3 = datetime.datetime.now()
    #print("get inter E between ", fraglabels, "and () in %f s"%(t3-t2).seconds)
    return inter_RC


def inter_elecbg(mol, dm, coords, charges):
    r"""
    inter energy between electrons in real frags and background charges
    """
    if mol.cart:
        intor = 'int3c2e_cart'
    else:
        intor = 'int3c2e_sph'
    nao = mol.nao
    #max_memory = mol.max_memory - lib.current_memory()[0]
    #blksize = int(min(max_memory * 1e6 / 8 / nao**2, 200))
    #print(coords,'\n',len(coords))
    vc = np.zeros((nao, nao))
    for i in range(0, charges.size):
        fakemol = gto.fakemol_for_charges(coords[i:i+1])
        j3c = df.incore.aux_e2(mol, fakemol, intor=intor, aosym='s2ij')
        v = np.einsum('xk,k->x', j3c, -charges[i:i+1])
        v = lib.unpack_tril(v)
        #print(vc.shape, v.shape)
        vc += v
    #E = np.einsum('ij,ji', v, dm)
    #print("inter_elecbg: %f" % E)
    return vc

def inter_nucbg(mol, coords, charges):
    nuc = np.zeros(mol.natm)
    #print(coords)
    for j in range(mol.natm):
        q2, r2 = mol.atom_charge(j), mol.atom_coord(j)
        r = lib.norm(r2 - coords, axis=1)
        #print(r2)
        nuc[j] = q2 * (charges / r).sum()
    #print("inter_nucbg: %f" % nuc)
    return nuc

def inter_nucbg_f(mol, atm_insub, coords, charges, chg_insub):
    #nuc = np.zeros(mol.natm)
    nuc = 0.0
    #print(coords)
    for j in atm_insub:
        q2, r2 = mol.atom_charge(j-1), mol.atom_coord(j-1)
        for k in chg_insub:
            #r = lib.norm(r2 - coords, axis=1)
            #print(r2)
            r = lib.norm(r2 - coords[k-1])
            q1 = charges[k-1]
            nuc += q2 * (q1 / r)
    #print("inter_nucbg: %f" % nuc)
    return nuc

def inter_bgbg(coords1, charges1, coords2, charges2):
    E = np.zeros(len(charges1))
    chg1 = len(coords1)
    chg2 = len(coords2)
    for i in range(chg1):
        for j in range(chg2):
            r = np.linalg.norm(coords1[i] - coords2[j])
            E[i] += charges1[i] * charges2[j] / r
    return E

