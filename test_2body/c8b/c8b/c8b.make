# make -s -j[N] -f c8b.make
fname = c8b
keyname = ${fname}.keys
files = c8b_1.gebf \
	c8b_2.gebf \
	c8b_3.gebf

all: $(files)
$(files):%.gebf:%.sxyz
	gauinput.py ${keyname} $(@:%.gebf=%.sxyz) g16 second
	g16 $(@:%.gebf=%.gjf) $(@:%.gebf=%.log)
	gebf.py $(@:%.gebf=%.chk) >& ${fname}.pr
