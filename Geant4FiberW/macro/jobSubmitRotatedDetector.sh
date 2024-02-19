#!/bin/bash

# option usage
# ./KotoEMCal -$option1=$arg1 -$option2=$arg2  ...  or ./KotoEMCal "-$option1 $arg1" "-$option2 $arg2" ...

# option list
# 1. useGPS: [true or false]  Whether use General Particle Gun or not
# 2. mac: [Name of macro file] (string)
# 3. seed: [Random seed for generation] (long int)
# 4. nEvents: [Number of events] (int)
# 5. outFileName: [Name of output root file] (string)

# option list for PrimaryGenerator
# 6. particle: [name of particle] (string)
# 7. energyOpt: [fixed or uniform or step]
# 7.a(fixed): [energy or momentum]: [Kinetic energy or momentum in unit of MeV] (double)
#                e.g. "-energyOpt=fixed -momentum=1000"
# 7.b(uniform): [energyMin and energyMax] or [momentumMin and MomentumMax] (double)
#                e.g. "-energyOpt=uniform -energyMin=0 -energyMax=1000" or "-energyOpt=uniform -momentumMin=0 -momentumMax=1000"
# 7.c(step): [energyBegin(double) and energyStepSize(double) and energyNstep(int)] or [momentumBegin(double) and momentumStepSize(double) and momentumNstep(int)]
#                e.g. "-energyOpt=step -momentumBegin=100 -momentumStepSize=100 -momentumNstep=10"
#
# 8. angleOpt: [fixed or uniformSolid or uniformTheta or stepTheta]
# 8.a(fixed): [theta and phi (deg)](double) or [directionX and directionY and directionZ](not need to be normalized)
#                e.g. "-angleOpt=fixed -theta=10 -phi=0" or "-angleOpt=fixed -directionX=0.2 -directionY=-0.3 -directionZ=1"
# 8.b(uniformSolid): [thetaMin and thetaMax (deg)](double)  uniform cos(theta) generation
#                e.g. "-angleOpt=uniformSolid -thetaMin=0 -thetaMax=50"
# 8.c(uniformTheta): [thetaMin and thetaMax (deg)](double)  uniform theta generation
#                e.g. "-angleOpt=uniformTheta -thetaMin=0 -thetaMax=50"
# 8.d(stepTheta): [thetaBegin(double) and thetaStepSize(double) and thetaNstep(int)] unit in deg
#                e.g. "-angleOpt=stepTheta -thetaBegin=0 -thetaStepSize=5 -thetaNstep=9"
# 9.(optional) fixedPhi: [Azimuthal angle phi in unit of deg] (double)
#
# 10. position [minX and maxX] or posX]: [Beam Generation position in unit of mm] (double) (same for Y and Z)
#                e.g. "-minX=-5 -maxX=5 -minY=-5 -maxY=5 -posZ=0"
# 11. (optional) detRotationAngle: Detector rotation angle (for ELPH beam test setup)


outputPATH="../root_for_ELPH"
mkdir -p $outputPATH

optionParticle="-particle=e+"
optionNEvents="-nEvents=10000"
optionOthers="-useGPS=false"

#for rotAngle in 0 15 30 45; do
for rotAngle in 30 20 10 0; do
  optionDetRotate="-detRotationAngle=${rotAngle}"

  for energy in 800 600 400 200; do
    optionEnergy="-energyOpt=fixed -momentum=${energy}"
    for ((seed = 1; seed <= 10; seed++)); do
      seed_prefix0=$(printf "%04d\n" $seed)
      optionSeed="-seed=$seed"

      optionAngle="-angleOpt=uniformTheta -thetaMin=0 -thetaMax=50"

      optionPosition="-minX=-5 -maxX=5 -minY=-5 -maxY=5 -posZ=0"
      optionOutput="-outFileName=$outputPATH/positron${energy}MeV_rotated${rotAngle}_random_position5mm_${seed_prefix0}"
      bsub -q s ./KotoEMCal $optionEnergy $optionAngle $optionPosition $optionParticle $optionNEvents $optionSeed $optionOutput $optionDetRotate $optionOthers
      optionPosition="-posX=0 -posY=0 -posZ=0"
      optionOutput="-outFileName=$outputPATH/positron${energy}MeV_rotated${rotAngle}_random_position0_${seed_prefix0}"
      bsub -q s ./KotoEMCal $optionEnergy $optionAngle $optionPosition $optionParticle $optionNEvents $optionSeed $optionOutput $optionDetRotate $optionOthers

      theta=0
      phi=0
      optionAngle="-angleOpt=fixed -theta=${theta} -phi=${phi}"
      optionEnergy="-energyOpt=fixed -momentum=${energy}"
      optionPosition="-minX=-5 -maxX=5 -minY=-5 -maxY=5 -posZ=0"
      optionOutput="-outFileName=$outputPATH/positron${energy}MeV_rotated${rotAngle}_theta${theta}_phi${phi}_position5mm_${seed_prefix0}"
      bsub -q s ./KotoEMCal $optionEnergy $optionAngle $optionPosition $optionParticle $optionNEvents $optionSeed $optionOutput $optionDetRotate $optionOthers

      optionPosition="-posX=0 -posY=0 -posZ=0"
      optionOutput="-outFileName=$outputPATH/positron${energy}MeV_rotated${rotAngle}_theta${theta}_phi${phi}_position0_${seed_prefix0}"
      bsub -q s ./KotoEMCal $optionEnergy $optionAngle $optionPosition $optionParticle $optionNEvents $optionSeed $optionOutput $optionDetRotate $optionOthers
    done
  done
done
