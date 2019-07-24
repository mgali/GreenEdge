#!/bin/bash
f=$1
for i in {1..203}
do
   mkdir $f/$(printf 'cast%03i' "${i}")
done
