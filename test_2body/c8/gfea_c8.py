from dftpart import gfea3

fr0 = gfea3.GFEA()
fr0.inputstyle = 'frg'
fr0.method = ['hf', '3-21g','charge']
fr0.gjfname = 'c8'
fr0.showinter = True
fr0.kernel()

