#!/bin/bash
make
#The structuring set for a 3x3 square
SQUARE="-1,-1/ -1,0/ -1,1/ 0,-1/ 0,0/ 0,1/ 1,-1/ 1,0/ 1,1/"
COORDS=$SQUARE
THRESH="120"
./filtering $1 $2 $3 $THRESH $COORDS
