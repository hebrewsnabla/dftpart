import numpy as np
import sys
import labc as la

def read_gfea(gfealog):
    with open(gfealog,'r') as f:
        data = f.readlines()
    for line in data:
        if 'atom_E' in line:
            start = data.index(line)+1
        if 'E_GFEA' in line:
            end = data.index(line)
    atom_E = []
    for line in data[start:end]:
        line = line.strip().split()
        for item in line:
            atom_E.append(float(item))
    atmE_gfea = np.array(atom_E)
    return atmE_gfea

def read_eda(edalog):
    with open(edalog,'r') as f:
        data = f.readlines()
    for line in data:
        if 'OK' in line:
            start = data.index(line)+1
    atom_E = []
    atomsymb = []
    for line in data[start:]:
        line = line.strip().split()
        if len(line)==3:
            atom_E.append(float(line[2]))
            atomsymb.append(line[0])
    atmE_eda = np.array(atom_E)
    return atmE_eda, atomsymb

def get_cenintot(labc):
    lab = la.labc_parser(labc)
    cen_intot = []
    for i in range(len(lab)):
        atomlist_lab = lab[i]
        cen_intot_i = []
        for label in atomlist_lab:
            if label != 0:
                cen_intot_i.append(label)
        cen_intot.append(cen_intot_i)
    return cen_intot


def compare(name, style='by_atom'):
    compare2(name+'-gfea.log', name+'-eda.log', name+'.labc', style)

def compare2(gfealog, edalog, labc, style='by_atom'):
    atmE_gfea = read_gfea(gfealog)
    atmE_eda, atomsymb = read_eda(edalog)
    print("atom       GFEA/a.u.         EDA/a.u.      diff/mH ")
    for i in range(len(atomsymb)):
        ss = ("%d"%(i+1)).ljust(3) + atomsymb[i]
        ss += ("%.10f"%atmE_gfea[i]).rjust(18)
        ss += ("%.10f"%atmE_eda[i]).rjust(18)
        diff = (atmE_gfea[i]-atmE_eda[i])*1000.0
        ss += ("%.3f"%diff).rjust(10)
        print(ss+'\n')
    if style=='by_frag':
        print("\nfrag       GFEA/a.u.         EDA/a.u.      diff/mH ")
        cen_intot = get_cenintot(labc)
        for cen in cen_intot:
            fragE_gfea = 0.0
            fragE_eda = 0.0
            for i in cen:
                fragE_gfea += atmE_gfea[i-1]
                fragE_eda += atmE_eda[i-1]
            diff = (fragE_gfea-fragE_eda)*1000.0
            ss = ("%d"%(cen_intot.index(cen)+1)).ljust(3) 
            ss += ("%.10f"%fragE_gfea).rjust(18)
            ss += ("%.10f"%fragE_eda).rjust(18)
            #diff = (atmE_gfea[i]-atmE_eda[i])*1000.0
            ss += ("%.3f"%diff).rjust(10)
            print(ss+'\n')


if __name__ == '__main__':
    if sys.argv[1]=='1':
        compare(sys.argv[2],sys.argv[2])
    elif sys.argv[1]=='3':
        compare2(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
