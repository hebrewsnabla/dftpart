from dftpart import scfeda

mollist = [
#'c12h26'
#'c2h6',
'c4h10-i','c6h14-ia','c6h14-ib','c6h14-ic','c6h14-t'
#'c8h18','c12h26'
]
basislist = [
#'6-31g',
'6-31gss','6-311gss'#,'6-311+gss',
#'cc-pvdz'
,'cc-pvtz',
#'def2-svp','def2-tzvp'
]
for mol in mollist:
    for basis in basislist:
        test1 = scfeda.EDA()
        test1.gjf = './'+mol+'.gjf'
        test1.method = ['m062x',basis]
        test1.output = mol.upper()+'_'+basis
        #test1.build()
        test1.kernel()
        

