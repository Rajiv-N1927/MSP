#!/bin/bash
make #Compile the cpp files
SRC="seq/seq1"
FILE1="f222"
FILE2="f223"
OUT="$FILE1-$FILE2"
B_EXT="25"
S_EXT="25"
DELTA="100"

lldb ./motion "$SRC/$FILE1.bmp" "$SRC/$FILE2.bmp" "$SRC/$OUT.bmp" $B_EXT $S_EXT $DELTA
