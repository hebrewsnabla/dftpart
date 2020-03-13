import numpy as np
from pyscf import lib

def gjf_parser(gjf):
    with open(gjf, 'r') as f:
        gjflines = f.readlines()
    geom = ""
    coords = []
    charges = []
    molcharge = 0
    spin = 0
    for line in gjflines:
        l = line.strip().split()
        if len(l) == 4 and l[0].isalpha() and ('lsqc' not in l):
            #logger.log(f,'--',line,'--')
            geom += line
        if len(l) == 4 and ((l[0][0].isdigit()) or (l[0][0] == '-')):
            coords.append(float(l[0]))
            coords.append(float(l[1]))
            coords.append(float(l[2]))
            charges.append(float(l[3]))
        if len(l) == 2:
            try:
                a = int(l[0])
            except:
                continue
            else:
                molcharge = int(l[0])
                spin = int(l[1]) - 1
    #oldcharges = charges
    charges = np.array(charges)
    #oldcoords = coords
    coords = np.array(coords).reshape(-1, 3) / lib.parameters.BOHR
    return geom, coords, charges, molcharge, spin
