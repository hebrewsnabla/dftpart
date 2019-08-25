import numpy as np
import sys

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

def compare(name):
    compare2(name+'-gfea.log', name+'-eda.log')

def compare2(gfealog, edalog):
    atmE_gfea = read_gfea(gfealog)
    atmE_eda, atomsymb = read_eda(edalog)
    print("atom    GFEA/a.u.         EDA/a.u.      diff/mH ")
    for i in range(len(atomsymb)):
        ss = atomsymb[i]
        ss += ("%.10f"%atmE_gfea[i]).rjust(18)
        ss += ("%.10f"%atmE_eda[i]).rjust(18)
        diff = (atmE_gfea[i]-atmE_eda[i])*1000.0
        ss += ("%.3f"%diff).rjust(10)
        print(ss+'\n')

if __name__ == '__main__'
    if sys.argv[1]=='1':
        compare(sys.argv[2])
    elif sys.argv[1]=='2':
        compare2(sys.argv[2], sys.argv[3])
