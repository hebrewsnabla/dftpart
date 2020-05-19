# HF/DFT Energy Decomposition and GFEA

* Prerequisites
    - Python 3.6 (and some common modules, which are provided in Anaconda3)
    - PySCF 
    - simplejson
    - GEBF (see [here](itcc.nju.edu.cn/lsqc))
* Installation
Specify the path to PySCF in Makefile, and
```
make all
```
then add
```
export PYTHONPATH=$PYTHONPATH:/path/above/dftpart/
```
in `.bashrc`
