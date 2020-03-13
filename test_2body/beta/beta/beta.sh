#! /usr/bin/env bash
export GAUSS_SCRDIR=/tmp/liwei01/gaussian/liaokang/
mkdir -p /tmp/liwei01/gaussian/liaokang/44443996/
make -s -j4 -f beta.make
rm -rf /tmp/liwei01/gaussian/liaokang/44443996
