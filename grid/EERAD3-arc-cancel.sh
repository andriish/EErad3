#!/bin/bash
#for SL7   
source  /cvmfs/grid.cern.ch/centos7-ui-4.0.3-1_umd4v3/etc/profile.d/setup-c7-ui-example.sh
echo "Cancel jobs from DB. One argument, the DB name, is needed."
cd $1
for b in $( arcstat -j jobs.xml -p  | grep gsiftp ); do
arckill -j jobs.xml $b
done


exit
   56 files in directory ./JRTMC000/herwig71openloopsmergeddipolec/130

  148 files in directory ./JRTMC000/sherpaopenloopsc/130
  304 files in directory ./JRTMC000/sherpaopenloopsc/35


   16 files in directory ./JRTMC000/sherpaopenloopsc/133
   10 files in directory ./JRTMC000/sherpaopenloopsc/161
    6 files in directory ./JRTMC000/sherpaopenloopsc/183
    8 files in directory ./JRTMC000/sherpaopenloopsc/172
    1 files in directory ./JRTMC000/sherpaopenloopsc/136


  225 files in directory ./JRTMC000/sherpaopenloopsl/172


    1 files in directory ./JRTMC000/sherpaopenloopsl/194
    1 files in directory ./JRTMC000/sherpaopenloopsl/189
    1 files in directory ./JRTMC000/sherpaopenloopsl/183
    1 files in directory ./JRTMC000/sherpaopenloopsl/200
    1 files in directory ./JRTMC000/sherpaopenloopsl/196
    1 files in directory ./JRTMC000/sherpaopenloopsl/192
    1 files in directory ./JRTMC000/sherpaopenloopsl/205
    1 files in directory ./JRTMC000/sherpaopenloopsl/206
    1 files in directory ./JRTMC000/sherpaopenloopsl/202
    1 files in directory ./JRTMC000/sherpaopenloopsl/207


