# make -s -j[N] -f c8.make
fname = c8
keyname = ${fname}.keys
files = c8_1.gebf \
	c8_2.gebf \
	c8_3.gebf

all: $(files)
$(files):%.gebf:%.sxyz
	gauinput.py ${keyname} $(@:%.gebf=%.sxyz) g16 second
	g16 $(@:%.gebf=%.gjf) $(@:%.gebf=%.log)
	gebf.py $(@:%.gebf=%.chk) >& ${fname}.pr
