##############################################################################################
#
# This code produces at set of validation plots for a given detector card. 
#
# In order to run this you need to compile Delphes with Pythia8 first, see: 
# 
# https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/Pythia8
#
# After you (re-)compiled Delphes with Pythia8 you are ready to go, execute from Delphes main dir:
#
# ./examples/validation.sh [detector_card] [number_of_events] 
# 
#  e.g.
# 
# ./examples/validation.sh delphes_card_CMS.tcl 100000
#
# Note that the more events you specify, the more accurate the controls plots will be ...
#
############################################################################################

#! /bin/sh

EXPECTED_ARGS=2
E_BADARGS=65

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./examples/validation.sh [detector_card] [number_of_events]"
  echo "for instance: ./examples/validation.sh delphes_card_CMS.tcl 10000"
  exit $E_BADARGS
fi

card=$1
nEvents=$2
validationCard=cards/validation_$card
output=validation_${card%.*}.root

sed 's/delphes_card_CMS.tcl/'$card'/g' cards/validation_card.tcl  > $validationCard

function runParticleGun {
  name=$1
  pid=$2
  cmnd="examples/Pythia8/configParticleGun_$name.cmnd"
  sed '/Main:spareMode1/s/=[[:space:]]*[0-9]*/= '$pid'/' examples/Pythia8/configParticleGun.cmnd > $cmnd
  sed '/Main:numberOfEvents/s/=[[:space:]]*[0-9]*/= '$nEvents'/' examples/Pythia8/configParticleGun.cmnd  > $cmnd 
  ./DelphesPythia8 $validationCard $cmnd delphes_ParticleGun_$name.root
}

runParticleGun electron 11
runParticleGun muon 13
runParticleGun photon 22
runParticleGun jet 1
runParticleGun bjet 5
runParticleGun taujet 15

./Validation delphes_ParticleGun_electron.root delphes_ParticleGun_muon.root delphes_ParticleGun_photon.root delphes_ParticleGun_jet.root delphes_ParticleGun_bjet.root delphes_ParticleGun_taujet.root $output
