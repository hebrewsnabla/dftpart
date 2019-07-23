from dftpart import scfeda

mollist = [
#'c12h26',
'c2h6',
#'c4h10','c6h14',
#'c8h18','c12h26',
#'c10h22'
]
basislist = [
#'6-31g',
'6-31gss',
#'6-311gss','6-311+gss',
#"6311g(2df,2pd)"
#'cc-pvdz','cc-pvtz',
#'def2-svp','def2-tzvp'
]
for mol in mollist:
    for basis in basislist:
        test1 = scfeda.EDA()
        test1.gjf = './'+mol+'.xyz'
        test1.method = ['m062x',basis]
        test1.output = mol.upper()+'_'+basis+'_0724'
        #test1.build()
        test1.kernel()
        

