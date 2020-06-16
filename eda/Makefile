FC = f2py
FCFLAGS_G = --fcompiler=gfortran --f90flags='-fopenmp' -lgomp
PYSCF = /home/liwei01/zcheng/pyscf/pyscf
FCFLAGS_I = -L$(PYSCF)/lib/deps/lib -lcint --fcompiler=intelem --compiler=intelem -liomp5

h1e_new: h1e_new.f90
	$(FC) $(FCFLAGS_G) -m h1e_new -c h1e_new.f90 
new_eda: new_eda.f90
	$(FC) $(FCFLAGS_I) -m new_eda -c new_eda.f90  

all: h1e_new new_eda

clean:
	rm *.so