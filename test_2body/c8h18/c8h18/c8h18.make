# make -s -j[N] -f c8h18.make
fname = c8h18
keyname = ${fname}.keys
files = c8h18_1.gebf \
	c8h18_2.gebf \
	c8h18_3.gebf

all: $(files)
$(files):%.gebf:%.sxyz
	gauinput.py ${keyname} $(@:%.gebf=%.sxyz) g16 second
	g16 $(@:%.gebf=%.gjf) $(@:%.gebf=%.log)
	gebf.py $(@:%.gebf=%.chk) >& ${fname}.pr
