#!/bin/bash
make #Compile the cpp files
SRC="seq/seq1"
FILE1="f222.bmp"
FILE2="f223.bmp"
B_EXT="25"
S_EXT="25"
DELTA="100"

./motion "$SRC/$FILE1" "$SRC/$FILE2" $B_EXT $S_EXT $DELTA
