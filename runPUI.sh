#!/bin/bash

make clean
make

args=("$@")
fname=$args

tail parameter.h

echo $fname
`./pui -s <$fname &>output.log &`
