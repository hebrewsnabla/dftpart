import os,sys

src = sys.argv[1] 
obj = src[:-4]
env = sys.argv[2]
if env=='wsr':
    lib="/home/wsr/pyscf/lib/"

omp=""
if sys.argv[3]==1:
    omp="--f90flags=\'-fopenmp\' -lgomp"
os.system("f2py -m "+obj+" -c "+src+" -L"+lib+" -lcint --fcompiler=gfortran "+omp)
