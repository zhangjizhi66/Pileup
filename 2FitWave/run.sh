#!/bin/bash

# ./run 314 270000 54   CSNSwave202304
# ./run [218-223] 200000 40   CSNSwave202310

runnum=$1

totalentry=$2

kernel=1
while(( $kernel<=$3 ))
do
   ./ana $runnum `expr $kernel \* $2 / $3 - $2 / $3` `expr $kernel \* $2 / $3 - 1` &
   #echo "./ana $runnum `expr $kernel \* $2 / $3 - $2 / $3` `expr $kernel \* $2 / $3 - 1`"
    
   let "kernel++"
done