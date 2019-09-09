#! /usr/bin/env bash
export GAUSS_SCRDIR=/tmp/liwei01/gaussian/liaokang/
mkdir -p /tmp/liwei01/gaussian/liaokang/69924092/
make -s -j4 -f 32ane.make
rm -rf /scratch/liwei01/gaussian/69924092
