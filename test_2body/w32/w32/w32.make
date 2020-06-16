# make -s -j[N] -f w32.make
fname = w32
keyname = ${fname}.keys
files = w32_1.gebf \
	w32_2.gebf \
	w32_3.gebf \
	w32_4.gebf \
	w32_5.gebf \
	w32_6.gebf \
	w32_7.gebf \
	w32_8.gebf \
	w32_9.gebf \
	w32_10.gebf \
	w32_11.gebf \
	w32_12.gebf \
	w32_13.gebf \
	w32_14.gebf \
	w32_15.gebf \
	w32_16.gebf \
	w32_17.gebf \
	w32_18.gebf \
	w32_19.gebf \
	w32_20.gebf \
	w32_21.gebf \
	w32_22.gebf \
	w32_23.gebf \
	w32_24.gebf \
	w32_25.gebf \
	w32_26.gebf \
	w32_27.gebf \
	w32_28.gebf \
	w32_29.gebf \
	w32_30.gebf \
	w32_31.gebf \
	w32_32.gebf \
	w32_33.gebf \
	w32_34.gebf \
	w32_35.gebf \
	w32_36.gebf \
	w32_37.gebf \
	w32_38.gebf \
	w32_39.gebf \
	w32_40.gebf \
	w32_41.gebf \
	w32_42.gebf \
	w32_43.gebf \
	w32_44.gebf \
	w32_45.gebf \
	w32_46.gebf \
	w32_47.gebf \
	w32_48.gebf \
	w32_49.gebf \
	w32_50.gebf \
	w32_51.gebf \
	w32_52.gebf \
	w32_53.gebf \
	w32_54.gebf \
	w32_55.gebf \
	w32_56.gebf \
	w32_57.gebf \
	w32_58.gebf \
	w32_59.gebf \
	w32_60.gebf \
	w32_61.gebf \
	w32_62.gebf \
	w32_63.gebf \
	w32_64.gebf \
	w32_65.gebf \
	w32_66.gebf \
	w32_67.gebf \
	w32_68.gebf \
	w32_69.gebf \
	w32_70.gebf \
	w32_71.gebf \
	w32_72.gebf \
	w32_73.gebf \
	w32_74.gebf \
	w32_75.gebf \
	w32_76.gebf \
	w32_77.gebf \
	w32_78.gebf \
	w32_79.gebf \
	w32_80.gebf \
	w32_81.gebf \
	w32_82.gebf \
	w32_83.gebf \
	w32_84.gebf \
	w32_85.gebf \
	w32_86.gebf \
	w32_87.gebf \
	w32_88.gebf \
	w32_89.gebf \
	w32_90.gebf \
	w32_91.gebf \
	w32_92.gebf \
	w32_93.gebf \
	w32_94.gebf \
	w32_95.gebf \
	w32_96.gebf \
	w32_97.gebf \
	w32_98.gebf \
	w32_99.gebf \
	w32_100.gebf

all: $(files)
$(files):%.gebf:%.sxyz
	gauinput.py ${keyname} $(@:%.gebf=%.sxyz) g16 second
	g16 $(@:%.gebf=%.gjf) $(@:%.gebf=%.log)
	gebf.py $(@:%.gebf=%.chk) >& ${fname}.pr
