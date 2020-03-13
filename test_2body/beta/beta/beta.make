# make -s -j[N] -f beta.make
fname = beta
keyname = ${fname}.keys
files = beta_1.gebf \
	beta_2.gebf \
	beta_3.gebf

all: $(files)
$(files):%.gebf:%.sxyz
	gauinput.py ${keyname} $(@:%.gebf=%.sxyz) g16 second
	g16 $(@:%.gebf=%.gjf) $(@:%.gebf=%.log)
	gebf.py $(@:%.gebf=%.chk) >& ${fname}.pr
