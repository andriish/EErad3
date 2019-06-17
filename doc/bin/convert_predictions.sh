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
