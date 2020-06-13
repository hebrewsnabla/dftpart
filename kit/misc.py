
def p2f(atm2bas_p):
    atm2bas_f = []
    for item in atm2bas_p:
        item_f = []
        for j in item:
            item_f.append(j+1)
        atm2bas_f.append(item_f)
    return atm2bas_f

def one2zero(alist):
    if isinstance(alist, list):
        blist = []
        for item in alist:
            blist.append(one2zero(item))
    elif isinstance(alist, int):
        blist = alist - 1
    return blist

def mat2dict(mat, thresh, f2layer):
    mdict = {}
    if mat.ndim==1:
        n = mat.shape[0]
        for i in range(n):
            if abs(mat[i]) > thresh:
                layer = f2layer[i]
                mdict[(i)] = [mat[i], (layer)]
    elif mat.ndim==2:
        n,m = mat.shape[0], mat.shape[1]
        for i in range(n):
            for j in range(i,n):
                ilayer = f2layer[i]
                jlayer = f2layer[j]
                if abs(mat[i,j]) > thresh: 
                    mdict[(i,j)] = [mat[i,j], (ilayer, jlayer)]
    return mdict

def dict_merge(*dicts):
    sum_dict = {}
    for d in dicts:
        for k,v in d.items():
            if k in sum_dict.keys():
                sum_dict[k] += v
            else:
                sum_dict[k] = v
    return sum_dict

def dict_cut(olddict, thresh):
    newdict = {}
    for k,v in olddict.items():
        if abs(v) > thresh:
            newdict[k] = v
    return newdict