#!/bin/bash

for energy in 200 500 1000; do
  for xy in "0" "5mm"; do
    for particle in "electron" "positron"; do
      for targetIndex in 0 1 ; do
        bsub -q lx -n2 python -u train.py $targetIndex $energy $xy $particle
      done
    done
  done
done

