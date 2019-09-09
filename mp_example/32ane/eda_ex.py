from dftpart import scfeda

fr0tot = scfeda.EDA()
fr0tot.method = ['hf','6-31gs']
fr0tot.output = '32ane_tot'
fr0tot.gjf = '32ane.gjf'
fr0tot.kernel()
