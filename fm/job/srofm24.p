#!/bin/bash

#$ -q single.q
#$ -j yes
#$ -cwd

t0=$(date +%s.%N)
t0_string=$(date)

#########################
for var in `seq 5 25`
do
    ./srofm 2.4 $var 0.25
done
#########################

t1=$(date +%s.%N)
t1_string=$(date)

t=$(echo "$t1 - $t0"|bc)
h=$(echo "($t/3600)"|bc)
m=$(echo "($t%3600)/60"|bc)
s=$(echo "($t%3600)%60"|bc)

echo "#--------------------------------------------------"
echo "# Job Start : $t0_string"
echo "# Job End   : $t1_string"
echo "# Elapsed time : ${h}h ${m}m ${s}s"
echo "#--------------------------------------------------"
