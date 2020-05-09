#! /usr/bin/env bash
export GAUSS_SCRDIR=/tmp/liwei01/gaussian/liaokang/
mkdir -p /tmp/liwei01/gaussian/liaokang/19214293/
make -s -j4 -f c8.make
rm -rf /tmp/liwei01/gaussian/liaokang/19214293
