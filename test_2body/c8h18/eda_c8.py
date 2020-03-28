from dftpart import scfeda

test1 = scfeda.EDA()
test1.gjf = 'c8h18.gjf'
test1.method = ['hf','6-31gss']
test1.output = 'c8h18'
test1.verbose = 9
#test1.build()
test1.showinter=True
test1.lso = 'c8h18/c8h18.lso'
test1.kernel()
        

