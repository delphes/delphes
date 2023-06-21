#!/bin/bash

#unset LD_LIBRARY_PATH
#unset PYTHONHOME
#unset PYTHONPATH
#source /cvmfs/sft.cern.ch/lcg/releases/LCG_88/Python/2.7.13/x86_64-slc6-gcc49-opt/Python-env.sh

PDG=${1}
NEVTS=${2}
PMIN=${3}
PMAX=${4}
ETAMIN=${5}
ETAMAX=${6}
SEED=${7}
OUTPUTDIR=${8}
HOMEDIR=${9}
DELPHES_PATH=${10}
DELPHES_CARD=${11}
CFG=${12}
LOG=${13}

JOBDIR="job_${PDG}_${SEED}"
mkdir $JOBDIR
cd $JOBDIR
cp ${CFG} .
cp -r ${HOMEDIR}/src .



################################################################################
## 1. run Delphes
################################################################################

OUTLHE="events_${PDG}_${SEED}.lhe"
python3 src/flatGunLHEventProducer.py \
  --pdg $PDG \
  --nevts $NEVTS \
  --pmin $PMIN \
  --pmax $PMAX \
  --etamin $ETAMIN \
  --etamax $ETAMAX \
  --seed $SEED \
  --output $OUTLHE \
    $LOG \


#cp ${OUTLHE} ${OUTPUTDIR}

################################################################################
## 2. run Delphes
################################################################################




cp -r ${DELPHES_PATH}/cards/* .
cp -r ${DELPHES_PATH}/examples/Pythia8/configLHE.cmnd .

OUTDELPHES="events_${PDG}_${SEED}.root"

BEAMSTR="Beams:LHEF = ${OUTLHE}"
echo $BEAMSTR >> configLHE.cmnd
echo 'Main:numberOfEvents = '${NEVTS} >> configLHE.cmnd
echo 'Main:timesAllowErrors = '${NEVTS} >> configLHE.cmnd

## run DelphesLHE if single that require no hadronisation
if [[ "$PDG" == "11" || "$PDG" == "13" || "$PDG" == "22" || "$PDG" == "130" || "$PDG" == "211" ]]; then
  ${DELPHES_PATH}/DelphesLHEF ${DELPHES_CARD} ${OUTDELPHES} ${OUTLHE}
else
  ${DELPHES_PATH}/DelphesPythia8 ${DELPHES_CARD} configLHE.cmnd ${OUTDELPHES}
fi

#cp ${OUTDELPHES} ${OUTPUTDIR}

################################################################################
## 3. run validation macro
################################################################################

cp ${DELPHES_PATH}/libDelphes.so .
OUTVAL="val_${PDG}_${SEED}.root"
mv src/run_loop.py .
python3 run_loop.py ${OUTDELPHES} ${OUTVAL} ${CFG}
cp ${OUTVAL} ${OUTPUTDIR}
