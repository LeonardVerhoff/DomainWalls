#!/bin/bash

MODIFY="t"

#sed -i "s/^$MODIFY = .*$/$MODIFY = 0.4/" input.dat

for ((VAL = 1 ; VAL < 10 ; VAL++)); 
do
    echo $VAL
    sed -i "s/^$MODIFY = .*$/$MODIFY = 0.$VAL/" input.dat
    ./monte_carlo

done
    
#./monte_carlo