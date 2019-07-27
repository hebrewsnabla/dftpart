from dftpart import scfeda

test1 = scfeda.EDA()
test1.gjf = './c2h6.xyz'
test1.method = ['m062x','cc-pvtz']
test1.output = 'C2H6_cc-pvtz_b'
test1.verbose = 6
#test1.build()
test1.kernel()

