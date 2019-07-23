'''
logger for GFEA

Jun/09/2019
'''

import numpy as np
from pyscf.tools.dump_mat import dump_rec
import simplejson as json
json.encoder.FLOAT_REPR = lambda o: "{:.10f}".format(o) #repr(round(o,10))

def slog(out, txt, *args):
    out.write(txt % args)
    out.write('\n')

def mlog(out, txt, *args, sep=" "):
    s = txt + ": "
    for item in args:
        s += sep
        s += json.dumps(item)
    out.write(s+"\n")

def log(out, txt, obj, label=None, ncol=10):
    if isinstance(obj, list):
        mat = np.array(obj)
        dump_mat(out, txt, mat, label, ncol)
    elif isinstance(obj, np.ndarray):
        dump_mat(out, txt, obj, label, ncol)
    else:
        out.write("obj type must be list or np.ndarray.\n")

def dump_1darray(out, txt, mat, label=None, ncol=10):
    out.write(txt + '\n')
    if label is 'n':
        label = ['#%d'%i for i in range(0,ncol)]
        out.write('%s\n' % ' '.join(label))
    for i in range(len(mat)):
        out.write(' %s' % str(mat[i]))
        if i%ncol == (ncol-1):
            out.write('\n')
    out.write('\n')    
# for 2D np.array, 
# use pyscf.tools.dump_mat.dump_rec(out, mat, label, label2, ncol, digits, start)


def dump_mat(out, txt, mat, label=None, ncol=10):
    if len(mat.shape) == 1:
        dump_1darray(out, txt, mat, label, ncol)
    elif len(mat.shape) == 2:
        if label is None:
            out.write(txt +'\n')
            dump_rec(out, mat, None, None, ncol)
    else:
        log(out, "Cannot dump \'%s\', list dim must be 1 or 2\n", txt)
