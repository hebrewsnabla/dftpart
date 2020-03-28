#! /usr/bin/env bash
export GAUSS_SCRDIR=/tmp/liwei01/gaussian/liaokang/
mkdir -p /tmp/liwei01/gaussian/liaokang/57988244/
make -s -j4 -f c8h18.make
rm -rf /tmp/liwei01/gaussian/liaokang/57988244
