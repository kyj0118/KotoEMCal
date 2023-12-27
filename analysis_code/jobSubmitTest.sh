#!/bin/bash

for energy in 200 500 1000; do
  for xy in "0" "5mm"; do
    for particle in "electron" "positron"; do
      bsub -q lx -n2 python -u test.py $energy $xy $particle
    done
  done
done
