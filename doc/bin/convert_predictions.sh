#!/bin/bash
mkdir -p processed/colorful/
BINING=OPAL3
epem1=predictions/ee3jet--LO.yoda
epem2=predictions/ee3jet--NLO.yoda
epem3=predictions/ee3jet--NLOcont.yoda
epem4=predictions/ee3jet--NNLO.yoda
epem5=predictions/ee3jet--NNLOcont.yoda
ls -lah  $epem1 $epem2 $epem3 $epem4 $epem5
cat $epem1 | grep -m1  -A 500 "Path=/alphasEEC/EEC_"$BINING  | grep -m1 -B 300 "END YODA_SCATTER2D"  |grep '0e' | sed 's@\t@        @g' | tr -s ' ' > processed/colorful/LO.dat
cat $epem2 | grep -m1  -A 500 "Path=/alphasEEC/EEC_"$BINING  | grep -m1 -B 300 "END YODA_SCATTER2D"  |grep '0e' | sed 's@\t@        @g' | tr -s ' ' > processed/colorful/NLO.dat
cat $epem3 | grep -m1  -A 500 "Path=/alphasEEC/EEC_"$BINING  | grep -m1 -B 300 "END YODA_SCATTER2D"  |grep '0e' | sed 's@\t@        @g' | tr -s ' ' > processed/colorful/NLOcont.dat
cat $epem4 | grep -m1  -A 500 "Path=/alphasEEC/EEC_"$BINING  | grep -m1 -B 300 "END YODA_SCATTER2D"  |grep '0e' | sed 's@\t@        @g' | tr -s ' ' > processed/colorful/NNLO.dat
cat $epem5 | grep -m1  -A 500 "Path=/alphasEEC/EEC_"$BINING  | grep -m1 -B 300 "END YODA_SCATTER2D"  |grep '0e' | sed 's@\t@        @g' | tr -s ' ' > processed/colorful/NNLOcont.dat


cat processed/colorful/LO.dat | cut -f 1,4 -d' ' > 1.temp
cat processed/colorful/NLO.dat | cut -f 4 -d' ' > 2.temp
cat processed/colorful/NNLO.dat | cut -f 4,5 -d' ' > 3.temp

paste -d ' '  1.temp 2.temp 3.temp > processed/colorful/ALL.dat

mkdir -p processed/eerad3/


#cat ../log/E00.y1d8.iL0.LO.input_0.log | grep -m 1 -A 200 'EEC distribution' | grep -B 200 'sum' | tr -s ' ' | grep -v '[a-z]' | grep '.' > processed/eerad3/LOl.dat

cat ../calc/outxxxxxxxxx.xNNLO.ELa| tr -s ' ' | cut -f 2,3,6 -d' ' > processed/eerad3/LO.dat
cat ../calc/outxxxxxxxxx.xNNLO.ELa| tr -s ' ' | cut -f 2,4,6 -d' ' > processed/eerad3/NLO.dat
cat ../calc/outxxxxxxxxx.xNNLO.ELa | tr -s ' '| cut -f 2,5,6 -d' ' > processed/eerad3/NNLO.dat

cat ../calc/outxxxxxxxxx.xNNLO.ELa | tr -s ' '| cut -f 2,3,4,5,6 -d' ' >processed/eerad3/ALL.dat


paste -d ' '  processed/colorful/ALL.dat processed/eerad3/ALL.dat > processed/BOTH.dat

export PATH=/afs/cern.ch/sw/XML/texlive/latest/bin/x86_64-linux/:$PATH
cd doc
#pdflatex   comparison.tex