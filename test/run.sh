#!/bin/bash

for f in `ls *.jl`
do
    echo ">>> Running test $f..."
    julia $f
done
