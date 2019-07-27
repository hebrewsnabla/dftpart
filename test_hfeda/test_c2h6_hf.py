from dftpart import scfeda

test1 = scfeda.EDA()
test1.gjf = '../1.xyz'
test1.method = ['hf','sto-3g']
test1.output = 'C2H6_hf'
#test1.build()
test1.kernel()

