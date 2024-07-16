#!/bin/bash

prot=$1
target=$2

for (( i=1 ; i<=$3 ; i++ ));
do
sbatch submitp.sh ${prot} ${target} $i
done
