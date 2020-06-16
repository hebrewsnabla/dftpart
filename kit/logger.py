'''
logger for GFEA

Jun/09/2019
'''

import numpy as np
from pyscf.tools.dump_mat import dump_rec
import simplejson as json
json.encoder.FLOAT_REPR = lambda o: "{:.8f}".format(o) #repr(round(o,8))

def slog(out, txt, *args, endl=True):
    out.write(txt % args)
    if endl:
        out.write('\n')
    out.flush()

def mlog(out, txt, *args, sep=" ", endl=True):
    s = txt + ": "
    for item in args:
        s += sep
        s += json.dumps(item)
    if endl: s += '\n'
    out.write(s)
    out.flush()

def ilog(out, txt, idict):
    s = txt + ": \n"
    s += 'inter term      layer        energy\n'
    for k,v in idict.items():
        s += json.dumps(k)
        s += '    '
        lenv = len(v)//2
        if len(v)==1:
            s += '            %.6f\n' % v[0]
            continue
        for i in range(lenv):
            if i>0:
                s+= '             '
            s += json.dumps(v[2*i+1])
            s += '    %.6f\n' % v[2*i]
    out.write(s)
    out.flush()
    

def log(out, txt, obj, label=None, ncol=10, digits=8, newl=True):
    if isinstance(obj, list):
        mat = np.array(obj)
        dump_mat(out, txt, mat, label, ncol, digits)
    elif isinstance(obj, np.ndarray):
        dump_mat(out, txt, obj, label, ncol, digits, newl)
    else:
        out.write("obj type must be list or np.ndarray.\n")
    out.flush()

def dump_1darray(out, txt, mat, label=None, ncol=10, digits=10, newl=True):
    if newl:
        out.write(txt + '\n')
    if label is 'n':
        label = ['#%d'%i for i in range(0,ncol)]
        out.write('%s\n' % ' '.join(label))
    templ = " %" + ".%df " % digits
    for i in range(len(mat)):
        out.write(templ % mat[i])
        if i%ncol == (ncol-1):
            out.write('\n')
    out.write('\n')    
    #out.flush()
# for 2D np.array, 
# use pyscf.tools.dump_mat.dump_rec(out, mat, label, label2, ncol, digits, start)


def dump_mat(out, txt, mat, label=None, ncol=10, digits=10, newl=True):
    if len(mat.shape) == 1:
        dump_1darray(out, txt, mat, label, ncol, digits, newl)
    elif len(mat.shape) == 2:
        if label is None:
            out.write(txt +'\n')
            dump_rec(out, mat, None, None, ncol, digits)
    else:
        log(out, "Cannot dump \'%s\', list dim must be 1 or 2\n", txt)
