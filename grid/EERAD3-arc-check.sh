#!/bin/bash
#for SL7   
source  /cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v4/etc/profile.d/setup-c7-ui-example.sh
echo "Check jobs from DB. One argument, the DB name, is needed."
a=$1
if [ "$#" -ne 1 ]; then
    echo "One argument, the DB name, is needed."
    a=$(ls -tr | grep ARC_ | tail -n 1)
fi

cd $a
for b in $( cat sites.txt | sort | uniq ); do
arcsync -j jobs.xml $b
done
arcstat -j jobs.xml