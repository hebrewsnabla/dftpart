
#export lsroot=/home/liwei01/liaokang/software/lsqc_bin
export lsroot=/home/liwei01/wsr/lsqc-2.4.2020Mar-build
export PATH=$lsroot/bin:$PATH
export PYTHONPATH=$lsroot/bin:$PYTHONPATH
export LSQC_SCRDIR=/tmp/liwei01/gaussian/liaokang
export CIM_BASDIR=$lsroot/basis
export LD_LIBRARY_PATH=$lsroot/lib:$LD_LIBRARY_PATH

lsqc -v
bsub -n 24 -q $2 lsqc $1 &

