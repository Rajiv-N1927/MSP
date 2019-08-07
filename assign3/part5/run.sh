#!/bin/bash
make #Compile the cpp files
SRC="seq/seq1"
FILE1="f222"
FILE2="f223"
OUT="$FILE1-$FILE2"
EXT="11"
SA_EXT="30"
SIG="2"
NKP="5000"

lldb ./motion "$SRC/$FILE1.bmp" "$SRC/$FILE2.bmp" "$SRC/$OUT.bmp" $EXT $SA_EXT $SIG $NKP
