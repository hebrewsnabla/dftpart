from dftpart import scfeda

test1 = scfeda.EDA()
test1.gjf = 'c8b.gjf'
test1.method = ['hf','6-31gss']
test1.output = 'c8b'
test1.verbose = 9
#test1.build()
test1.showinter=True
test1.lso = 'c8b/c8b.lso'
test1.kernel()
        

