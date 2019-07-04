#!/bin/bash
source  /cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/etc/profile.d/setup-c7-ui-example.sh
#export LCG_GFAL_INFOSYS=grid-bdii.desy.de:2170
## lcg-info --vo zeus --list-ce
#lcg-info --list-ce  --vo zeus --attrs  OS,OSRelease,OSVersion
#
#voms-proxy-init -voms zeus --vomses ~/vomses/zeus-grid-voms.desy.de -valid=500:00
#declate -a SITES = ( creamce.gina.sara.nl:8443/cream-pbs-medium  atlas-cream02.na.infn.it:8443/cream-pbs-atlas 
echo "Submits the jobs. Zero arguments."
set -x

mkdir -p GP
cp ../eerad3* ./GP
cp -R cards ./GP
tar cfz EERAD3000.tgz GP

TOP=$(pwd)
TIME=$(date +%s)
TMPSUBMIT=$TOP/ARC"_"$TIME
mkdir -p $TMPSUBMIT
cd $TMPSUBMIT
echo $TIME > submit.txt
touch sites.txt

uberftp  grid-gftp2.rzg.mpg.de "cd /pnfs/rzg.mpg.de/data/zeus/group/andriish; mkdir EERAD3000"
uberftp  grid-gftp2.rzg.mpg.de "put ../EERAD3000.tgz /pnfs/rzg.mpg.de/data/zeus/group/andriish/EERAD3000/EERAD3000.tgz"
cp $TOP/EERAD3-arc-main.sh ./
chmod +x EERAD3-arc-main.sh
echo '
&
(executable = "EERAD3-arc-main.sh")
(stdout = "std.out")
(stderr = "std.err")
(inputFiles = ( "dir.txt" "") 
              ("EERAD3000.tgz" "gsiftp://grid-gftp2.rzg.mpg.de:2811/pnfs/rzg.mpg.de/data/zeus/group/andriish/EERAD3000.tgz")  
              )              
'> ./grid.jdl.template

CARDS=( E00.y1d8.iL0.LO.input   
        Etx.y1d8.iN1.NLOicol1.input  Etx.y1d8.iN2.NLOicol2.input Etx.y1d8.iN3.NLOicol3.input  
        Etx.y1d6.iZ1.NNLOicol1.input  
        Etx.y1d6.iZ2.NNLOicol2.input  
        Etx.y1d6.iZ3.NNLOicol3.input  
        Etx.y1d6.iZ4.NNLOicol4.input  
        Etx.y1d6.iZ5.NNLOicol5.input 
        Etx.y1d6.iZ6.NNLOicol6.input  
        )

for CARD in "${CARDS[@]}"
do

for i in `seq 200 205`;
do
rm -rf dir.txt
touch dir.txt
echo   $CARD     >> dir.txt
echo   $i        >> dir.txt
cat grid.jdl.template >grid.jdl
echo '
(outputFiles = ( "'$CARD'_'$i'.tar.gz" "gsiftp://grid-gftp2.rzg.mpg.de:2811/pnfs/rzg.mpg.de/data/zeus/group/andriish/EERAD3000/'$CARD'_'$i'.tgz")
              )
' >>grid.jdl
arcsub -j jobs.xml  -f  grid.jdl -z $TOP/client.conf
done

done