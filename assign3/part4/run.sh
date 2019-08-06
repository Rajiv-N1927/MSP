#!/bin/bash
make #Compile the cpp files
SRC="seq/seq1"
FILE1="f222"
OUT="$FILE1-out"
EXT="11"
SIG="2"
NKP="100"

./motion "$SRC/$FILE1.bmp" "$SRC/$OUT.bmp" $EXT $SIG $NKP
