#! /usr/bin/env bash
export GAUSS_SCRDIR=/tmp/liwei01/gaussian/liaokang/
mkdir -p /tmp/liwei01/gaussian/liaokang/62621485/
make -s -j2 -f c8b.make
rm -rf /tmp/liwei01/gaussian/liaokang/62621485
