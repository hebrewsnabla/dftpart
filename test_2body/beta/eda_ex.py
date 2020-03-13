from dftpart import scfeda
import numpy as np

frg = [
    [33,34,35,36,37,38,59,60,61,62,63,64,65,66,67,68,69],
    [1,2,3,4,5,6,7,8,9,10,28,29,51,52,53,54,55,56,57,58],
    [11,12,13,14,15,16,17,18,19,20,30,31,43,44,45,46,47,48,49,50],
    [21,22,23,24,25,26,27,32,39,40,41,42,70,71]
]


test1 = scfeda.EDA()
test1.gjf = 'beta.gjf'
test1.method = ['hf','6-31gss']
test1.output = 'beta'
#test1.cart = False
#test1.build()
test1.showinter=True
atm_E, totE, conv, inter_terms = test1.kernel()
        
RR1, RR2, RR3, RR4 = inter_terms

E1 = np.zeros((4))
E2 = np.zeros((4,4))
#E3 = np.zeros((4,4,4))
#E4 = np.zeros((4,4,4,4))
E3 = RR3
E4 = RR4

for f in range(4):
    for i in frg[f]:
        E1[f] += RR1[i-1]

for f1 in range(4):
    for f2 in range(f1,4):
        if f2==f1:
            for i in frg[f1]:
                for j in frg[f2]:
                    if i < j:
                        E1[f1] += RR2[i-1,j-1]
                    else:
                        E1[f1] += RR2[j-1,i-1]
        else:
            for i in frg[f1]:
                for j in frg[f2]:
                    if i < j:
                        E2[f1,f2] += RR2[i-1,j-1]
                    else:
                        E2[f1,f2] += RR2[j-1,i-1]
print(E1)
print(E2)
print(E3)
print(E4)
