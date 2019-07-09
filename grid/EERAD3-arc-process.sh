#!/bin/bash

sh EERAD3-arc-get.sh
for a in $(ls output/*tgz); do
tar -xzf $a
done
mkdir -p logcombine
cd calc
set -x
for a in $(ls -1 ../cards/E*combine); do
b=$(basename $a)
../../eerad3_combine -i $a  &> ../logcombine/$b'.log' &
done
wait
../../eerad3_dist -i ../cards/eerad3_dist.input
cd ..
cp calc/outxxxxxxxxx.xNNLO.* ../calc
cd ../doc/
sh bin/convert_predictions.sh
cd doc/
export PATH=/afs/cern.ch/sw/XML/texlive/latest/bin/x86_64-linux/:$PATH
pdflatex comparison.tex
cd ..
svn commit -m "OK" ./