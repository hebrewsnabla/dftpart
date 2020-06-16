from dftpart import scfeda

test1 = scfeda.EDA()
test1.gjf = 'w19.gjf'
test1.method = ['hf','6-31gs', 'cart']
test1.output = 'w19'
test1.verbose = 9
#test1.build()
test1.showinter=True
test1.lso = 'w19/w19.lso'
test1.kernel()
        

