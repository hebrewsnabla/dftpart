# make -s -j[N] -f w19.make
fname = w19
keyname = ${fname}.keys
files = w19_1.gebf \
	w19_2.gebf \
	w19_3.gebf \
	w19_4.gebf \
	w19_5.gebf \
	w19_6.gebf \
	w19_7.gebf

all: $(files)
$(files):%.gebf:%.sxyz
	gauinput.py ${keyname} $(@:%.gebf=%.sxyz) g16 second 4
	g16 $(@:%.gebf=%.gjf) $(@:%.gebf=%.log)
	gebf.py $(@:%.gebf=%.chk) >& ${fname}.pr
