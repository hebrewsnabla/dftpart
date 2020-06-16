#! /usr/bin/env bash
export GAUSS_SCRDIR=/tmp/liwei01/gaussian/wsr/
mkdir -p /tmp/liwei01/gaussian/wsr/89961990/
make -s -j2 -f w32.make
rm -rf /tmp/liwei01/gaussian/wsr/89961990
