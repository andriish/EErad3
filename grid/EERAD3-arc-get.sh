#!/bin/bash
#for SL7   
source  /cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v4/etc/profile.d/setup-c7-ui-example.sh
mkdir -p output
cd output
uberftp  grid-gftp2.rzg.mpg.de "cd /pnfs/rzg.mpg.de/data/zeus/group/andriish/EERAD3000; ls" | tr -s ' ' | tr -d '\r'| cut -f 9 -d' ' | grep input | sort > computed.txt
ls -1 | sort > downloaded.txt

comm -1 -3 downloaded.txt computed.txt > todownload.txt
set -x
for a in $(cat todownload.txt); do
echo "->"$a"<-"
uberftp  grid-gftp2.rzg.mpg.de "cd /pnfs/rzg.mpg.de/data/zeus/group/andriish/EERAD3000; get "$a";"

done