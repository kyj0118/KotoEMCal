#!/bin/bash

OUTPUTPATH=../root

mkdir -p $OUTPUTPATH
for pbwid  in 150 170
do
    for (( seed=1; seed<=10; seed++ ))
    do
	for energy in 100 200 300 400 500 600 700 800 900 1000 1500 2000
	do
	    seed0=`printf "%04d\n" $seed`
	    echo bsub -q s ./exampleB5_Pb${pbwid}um RandomGeneration.mac $OUTPUTPATH/gamma${energy}MeV_Pb${pbwid}um_random_${seed0} $seed $energy 0
	done
    done
    
    
    for (( seed=1; seed<=9; seed++ ))
    do
	for energy in 100 200 300 400 500 600 700 800 900 1000 1500 2000
	do
	    seed0=`printf "%04d\n" $seed`
	    echo bsub -q s ./exampleB5_Pb${pbwid}um RandomGeneration100k.mac $OUTPUTPATH/gamma${energy}MeV_Pb${pbwid}um_stepdeg_${seed0} $seed $energy 1
	done
    done


    # for energy training
    for (( seed=11; seed<=50; seed++ ))
    do
	seed0=`printf "%04d\n" $seed`
	bsub -q s ./exampleB5_Pb${pbwid}um RandomGeneration.mac $OUTPUTPATH/gammaUniformEnergy_Pb${pbwid}um_random_${seed0} $seed 0 2
    done

done

