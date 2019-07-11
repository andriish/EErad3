#!/bin/bash

#sh EERAD3-arc-get.sh
#for a in $(ls output/*tgz); do
#tar --skip-old-files -xzf  $a
#done
mkdir -p logcombine
cd calc
set -x
for a in $(ls -1 ../cards/E*combine | grep '.L0'); do
b=$(basename $a)
../../eerad3_combine -i $a  -R 7.0 &> ../logcombine/$b'.log' &
done

for a in $(ls -1 ../cards/E*combine | grep  '.NL0' ); do
b=$(basename $a)
../../eerad3_combine -i $a  -R 10.0 &> ../logcombine/$b'.log' &
done

for a in $(ls -1 ../cards/E*combine | grep  '.NNL0' ); do
b=$(basename $a)
../../eerad3_combine -i $a  -R 10.0 &> ../logcombine/$b'.log' &
done


wait
#exit
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