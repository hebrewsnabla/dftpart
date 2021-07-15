# HF/DFT Energy Decomposition and GFEA

## Prerequisites
* Python3, numpy
* PySCF 
* simplejson

EDA/GFEA also requires
* libcint
* GEBF (see [here](https://itcc.nju.edu.cn/lsqc))

## Features
* Energy decomposition up to 4-body
* Mayer decomposition (1-body and 2-body)
* Generalized Fragment Energy Assembler (GFEA)
* Availability:
  - HF, pure DFT, hybrid DFT

## Installation
If you only need MayerEDA, skip the compilation.

For EDA/GFEA, specify the path to libcint library in Makefile, and
```
make all
```

Finally, add
```
export PYTHONPATH=$PYTHONPATH:/path/above/dftpart/
```
in `.bashrc`
