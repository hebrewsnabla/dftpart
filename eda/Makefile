FC = f2py
FCFLAGS_G = --fcompiler=gfortran --f90flags='-fopenmp' -lgomp

## huawei cluster
PYSCF = /share/home/srwang/pyscf-1.7.5/pyscf
CINT_FLAGS = -L$(PYSCF)/lib/deps/lib -I$(PYSCF)/lib/deps/include -lcint

## develop machine
#CINT_FLAGS = -L/home/jean/pyscf_deps/lib -lcint 

#FCFLAGS_I = --fcompiler=intelem --compiler=intelem -liomp5

h1e_new: h1e_new.f90
	$(FC) $(FCFLAGS_G) -m h1e_new -c h1e_new.f90 
new_eda: new_eda.f90
	$(FC) $(FCFLAGS_G) $(CINT_FLAGS) -m new_eda -c new_eda.f90  

all: h1e_new new_eda

clean:
	rm *.so
