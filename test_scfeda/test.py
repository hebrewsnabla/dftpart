from dftpart import scfeda

test1 = scfeda.EDA()
test1.gjf = './c4h10.xyz'
test1.method = ['m062x','6-31g']
test1.output = 'C4H10_6-31g_b'
test1.verbose = 6
#test1.build()
test1.kernel()

