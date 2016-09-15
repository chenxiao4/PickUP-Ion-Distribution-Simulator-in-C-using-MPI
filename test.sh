#!/bin/bash

ofile="output"
infile="in"
fmt=".log"

for i in 98 99 00 01 02 03 04 05 06 07 08 09 10
do
    fname="$ofile$i$fmt"
    iname="$infile$i"
done