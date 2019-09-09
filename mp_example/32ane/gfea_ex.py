from dftpart import gfea3

fr0 = gfea3.GFEA()
fr0.inputstyle = 'frg'
fr0.method = ['hf', '6-31gs']
fr0.gjfname = '32ane'
fr0.output = '32ane'
fr0.kernel()

