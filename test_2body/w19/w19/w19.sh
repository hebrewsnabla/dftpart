#! /usr/bin/env bash
export GAUSS_SCRDIR=/tmp/liwei01/gaussian/liaokang/
mkdir -p /tmp/liwei01/gaussian/liaokang/64063490/
make -s -j4 -f w19.make
rm -rf /tmp/liwei01/gaussian/liaokang/64063490
