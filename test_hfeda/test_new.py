from dftpart import scfeda

ok = ['FAIL','OK']
mollist = [
#'c12h26',
'c2h6',
'c4h10','c6h14',
'c8h18','c12h26',
'c10h22'
]
basislist = [
#'6-31g',
'6-31gss','6-311gss','6-311+gss',
#"6311g(2df,2pd)"
'cc-pvdz','cc-pvtz',
'def2-svp','def2-tzvp'
]
with open("CH3_new.txt",'a') as f:
    for basis in basislist:
        for mol in mollist:
            test1 = scfeda.EDA()
            test1.gjf = './'+mol+'.xyz'
            test1.method = ['m062x',basis]
            test1.output = mol.upper()+'_'+basis+'_new'
            #test1.build()
            atme, tote, conv = test1.kernel()
            CH3E = atme[:4].sum()
            f.write("%s %.10f  %s\n" % (test1.output,CH3E,ok[conv]))
