from dftpart import scfeda

test1 = scfeda.EDA()
test1.gjf = 'c8.xyz'
test1.method = ['hf','3-21g']
test1.output = 'c8'
#test1.build()
test1.showinter=True
test1.kernel()
        

