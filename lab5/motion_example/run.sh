#!/bin/bash
make #Compile the cpp files
SRC="calendar"
FILE1="mb14.bmp"
FILE2="mb15.bmp"
OUT="out.bmp"

lldb ./motion "$SRC/$FILE1" "$SRC/$FILE2" "$SRC/$OUT"
