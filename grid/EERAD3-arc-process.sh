#!/bin/bash
for a in $(ls output/*tgz); do
tar -xvzf $a
done
mkdir -p logcombine
cd calc
set -x
for a in $(ls -1 ../cards/E*combine); do
b=$(basename $a)
../../eerad3_combine -i $a  &> ../logcombine/$b'.log'
done

../../eerad3_dist -i ../cards/eerad3_dist.input
