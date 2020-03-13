from dftpart import scfeda
import numpy as np

frg = [[1,2,3,4,5,6,7,8,9,10,28,29,51,52,53,54,55,56,57,58],
    [11,12,13,14,15,16,17,18,19,20,30,31,43,44,45,46,47,48,49,50],
    [21,22,23,24,25,26,27,32,39,40,41,42,70,71],
    [33,34,35,36,37,38,59,60,61,62,63,64,65,66,67,68,69]]
cen_intot1 = [1, 2, 3, 5, 6, 7, 8, 10, 51, 55, 33, 35, 36, 37, 59, 63, 66, 4, 9, 28, 29, 34, 38, 52, 53, 54, 56, 57, 58, 60, 61, 62, 64, 65, 67, 68, 69]
cen_insub1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 21, 22, 23, 24, 25, 26, 27, 28, 29, 33, 34, 37, 38, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58]
cen_intot2 = [11, 12, 13, 15, 16, 17, 18, 20, 43, 47, 21, 22, 23, 25, 26, 27, 39, 70, 14, 19, 24, 30, 31, 32, 40, 41, 42, 44, 45, 46, 48, 49, 50, 71]
cen_insub2 = [11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 36, 37, 38, 40, 41, 42, 43, 44, 45, 46, 47, 48, 55]
molchgs1 = np.zeros((58))

tot2sub = {}
for i in range(len(cen_intot1)):
    tot2sub[cen_intot1[i]] = cen_insub1[i]
for i in range(len(cen_intot2)):
    tot2sub[cen_intot2[i]] = cen_insub2[i]


test1 = scfeda.EDA()
test1.gjf = 'beta_subsys/beta_1.gjf'
test1.method = ['hf','6-31gss', 'charge']
test1.output = 'beta_1'
#test1.cart = False
#test1.build()
test1.molchgs = molchgs1
test1.showinter=True
atm_E1, totE1, conv1, inter_terms1 = test1.kernel()        
RR1a, RR2a, RR3a, RR4a, bg2, bg3 = inter_terms1

#test2 = scfeda.EDA()
#test2.gjf = 'beta_subsys/beta_2.gjf'
#test2.method = ['hf','6-31gss']
#test2.output = 'beta_2'
##test1.cart = False
##test1.build()
#test2.showinter=True
#atm_E2, totE2, conv2, inter_terms2 = test1.kernel()
#RR1b, RR2b, RR3b, RR4b = inter_terms2

# sub1: [1,4] [2]
# sub2: [2,3] [1]

E1a = np.zeros((4))
E1b = np.zeros((4))
E2a = np.zeros((4,4))
E2b = np.zeros((4,4))
#E3 = np.zeros((4,4,4))
#E4 = np.zeros((4,4,4,4))
E3 = RR3a#,RR3b]
E4 = RR4a#,RR4b]

t = {
3: 0, 0: 1, 1: 2, 2: 3
}



for f in [0,3]:
    for i in frg[f]:
        ff = t[f]
        ii = tot2sub[i] 
        E1a[ff] += RR1a[ii-1]
#for f in [1,2]:
#    for i in frg[f]:
#        ff = t[f]
#        ii = tot2sub[i] 
#        E1b[ff] += RR1b[ii-1]


for f1 in [3,0,1]:
    for f2 in [3, 0,1]:
        if f2 == f1: 
            for i in frg[f1]:
                for j in frg[f2]:
                    ii = tot2sub[i]
                    jj = tot2sub[j]
                    ff1 = t[f1]
                    if ii < jj:
                        E1a[ff1] += RR2a[ii-1,jj-1]
                    else:
                        E1a[ff1] += RR2a[jj-1,ii-1]
        else:
            for i in frg[f1]:
                for j in frg[f2]:
                    ii = tot2sub[i]
                    jj = tot2sub[j]
                    ff1 = t[f1]
                    ff2 = t[f2]
                    if ii < jj:
                        E2a[ff1,ff2] += RR2a[ii-1,jj-1]
                    else:
                        E2a[ff1,ff2] += RR2a[jj-1,ii-1]
#for f1 in [0,1,2]:
#    for f2 in [0,1,2]:
#        if f2 == f1: 
#            for i in frg[f1]:
#                for j in frg[f2]:
#                    ii = tot2sub[i]
#                    jj = tot2sub[j]
#                    ff1 = t[f1]
#                    if ii < jj:
#                        E1b[ff1] += RR2b[ii-1,jj-1]
#                    else:
#                        E1b[ff1] += RR2b[jj-1,ii-1]
#        else:
#            for i in frg[f1]:
#                for j in frg[f2]:
#                    ii = tot2sub[i]
#                    jj = tot2sub[j]
#                    ff1 = t[f1]
#                    ff2 = t[f2]
#                    if ii < jj:
#                        E2b[ff1,ff2] += RR2b[ii-1,jj-1]
#                    else:
#                        E2b[ff1,ff2] += RR2b[jj-1,ii-1]

print(E1a)
#print(E1b)
print(E2a)
#print(E2b)
print(E3)
print(E4)
