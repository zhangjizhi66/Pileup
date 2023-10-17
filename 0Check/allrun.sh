#!/bin/bash

for file in ../../CSNSwave202310/event_*_*.root
do
    temp=${file#*_}
    runnum=${temp%%_*}
    
    ./ana $runnum &
done