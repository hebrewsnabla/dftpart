from dftpart import gfea3_mp

fr0 = gfea3_mp.GFEA()
fr0.inputstyle = 'frg'
fr0.method = ['hf', '6-311gss']
fr0.gjfname = '32ane'
fr0.output = '32ane'
fr0.rack = 'rack1'
fr0.nodes = 2
fr0.kernel()

