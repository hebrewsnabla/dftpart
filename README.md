# HF/DFT Energy Decomposition and GFEA

* Prerequisites
    - Python 3 (and some common modules, which are provided in Anaconda3)
    - PySCF 
    - simplejson
    - GEBF (see [here](https://itcc.nju.edu.cn/lsqc))
* Installation

Specify the path to libcint library in Makefile, and
```
make all
```
then add
```
export PYTHONPATH=$PYTHONPATH:/path/above/dftpart/
```
in `.bashrc`
