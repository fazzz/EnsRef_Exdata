#!/bin/sh

opt=( dummy sim exp )
nopt=${#opt[*]}
if [ $# -le `expr ${nopt} - 2` ]; then
    echo "USAGE: $0" ${opt[*]:1:${nopt}}
    echo "Your Command"$*
    exit
fi

num=1
while [ $num -le `expr ${nopt} - 1` ]; do
    eval ${opt[$num]}=$1
    shift 1
    num=`expr $num + 1`
done

awk -F"_" '{print $2}' ${sim} | awk -F"-" '{print $1}' > temp0

paste temp0 ${sim} > temp1

cat ${exp} | tr \[a-z\] \[A-Z\] | awk '{print $1$2 " " $3}' > temp2

echo "exp sim(ave) dihedral"

join -1 1 -2 1 temp2 temp1 | awk '{printf("%5.2f %5.2f %5.4f %40s\n",$2,$3,($2-$3)*($2-$3),$4)}'
