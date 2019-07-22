import os,sys

src = sys.argv[1] 
obj = src[:-4]
env = sys.argv[2]
if env=='wsr':
    lib="/home/wsr/pyscf/lib/"
elif env=='ubuntu':
    lib="/home/wsr/anaconda3/lib/python3.7/site-packages/pyscf/lib/"
elif env=='huawei':
    lib="/home/liwei01/liaokang/anaconda3/lib/python3.6/site-packages/pyscf/lib"
omp=""
if sys.argv[3]==1:
    omp="--f90flags=\'-fopenmp\' -lgomp"
os.system("f2py -m "+obj+" -c "+src+" -L"+lib+" -lcint --fcompiler=gfortran "+omp)
