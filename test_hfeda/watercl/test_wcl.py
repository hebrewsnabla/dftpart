from dftpart import scfeda

mollist = [
#'5','6','7'
'8','9'
#,'c4h10','c6h14'
]
basislist = [
#'6-31g',
'6-31gss','6-311gss',
#'6-311+gss',
'cc-pvdz','cc-pvtz',
'def2-svp','def2-tzvp'
]
for mol in mollist:
    for basis in basislist:
        test1 = scfeda.EDA()
        test1.gjf = './TIP4P-'+mol+'.xyz'
        test1.method = ['m062x',basis]
        test1.output = 'TIP4P-'+mol.upper()+'_'+basis
        #test1.build()
        test1.kernel()

