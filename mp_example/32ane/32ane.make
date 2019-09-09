# make -s -j[N] -f 32ane.make
fname = 32ane
keyname = ${fname}.keys
files = 32ane_1.gebf \
	32ane_2.gebf \
	32ane_3.gebf \
	32ane_4.gebf \
	32ane_5.gebf \
	32ane_6.gebf \
	32ane_7.gebf \
	32ane_8.gebf \
	32ane_9.gebf \
	32ane_10.gebf \
	32ane_11.gebf

all: $(files)
$(files):%.gebf:%.sxyz
	gauinput.py ${keyname} $(@:%.gebf=%.sxyz) g16 second
	g16 $(@:%.gebf=%.gjf) $(@:%.gebf=%.log)
	gebf.py $(@:%.gebf=%.chk) >& ${fname}.pr
