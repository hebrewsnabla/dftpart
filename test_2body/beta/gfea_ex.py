from dftpart import gfea3

fr0 = gfea3.GFEA()
fr0.inputstyle = 'frg'
fr0.method = ['hf', '6-31gss','charge']
fr0.gjfname = 'beta'
fr0.cart = False
fr0.showinter = True
fr0.kernel()

