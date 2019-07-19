from dftpart import scfeda

test1 = scfeda.EDA()
test1.gjf = '../1.xyz'
test1.method = ['m062x','6-31gss']
test1.output = 'C2H6_6-31gss'
#test1.build()
test1.kernel()

