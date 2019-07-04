#!/bin/bash
#for SL7   
source  /cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v4/etc/profile.d/setup-c7-ui-example.sh
echo "Cancel jobs from DB. One argument, the DB name, is needed."
cd $1
for b in $( arcstat -j jobs.xml -p  | grep gsiftp ); do
arckill -j jobs.xml $b
done





