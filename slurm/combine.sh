#!/bin/bash
set -x
mkdir -p calc
mkdir -p log
mkdir -p logcombine
cd calc
for a in $(ls -1 ../cards/E0*combine); do
b=$(basename $a)
../eerad3_combine -i $a  &> ../logcombine/$b'.log'
done

#exit

for a in $(ls -1 ../cards/Et*combine); do
b=$(basename $a)
../eerad3_combine -i $a  &> ../logcombine/$b'.log'
done

../eerad3_dist -i ../cards/eerad3_dist.input