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
# ./examples/validation.sh cards/delphes_card_CMS.tcl 100000
#
# Note that the more events you specify, the more accurate the controls plots will be ...
# This said, 500k events should be ok for most cases. 
#
############################################################################################

#! /bin/sh

EXPECTED_ARGS=2
E_BADARGS=65

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./examples/validation.sh [detector_card] [number_of_events]"
  echo "for instance: ./examples/validation.sh cards/delphes_card_CMS.tcl 10000"
  exit $E_BADARGS
fi

cardbase=$(basename $1)
carddir=$(dirname $1)
nEvents=$2
validationCard=$carddir/validation_$cardbase
output=validation_${cardbase%.*}.root
mainoutputdir=report_${cardbase%.*}
outputroot=report_${cardbase%.*}/root
cardlabel=${cardbase%.*}
version=$(cat VERSION)
outpdf=$mainoutputdir/${output%.*}.pdf

mkdir -p $outputroot
mkdir -p $mainoutputdir/www/fig

sed 's/delphes_card_CMS.tcl/'$cardbase'/g' cards/validation_card.tcl  > $validationCard

function runParticleGun {
  name=$1
  pid=$2
  cmnd="examples/Pythia8/configParticleGun_$name.cmnd"
  rootfile="particleGun_${name}_${cardlabel}.root"
  sed '/Main:spareMode1/s/=[[:space:]]*[0-9]*/= '$pid'/' examples/Pythia8/configParticleGun.cmnd > examples/Pythia8/tmp.cmnd
  sed '/Main:numberOfEvents/s/=[[:space:]]*[0-9]*/= '$nEvents'/' examples/Pythia8/tmp.cmnd  > $cmnd 
  ./DelphesPythia8 $validationCard $cmnd $outputroot/$rootfile
  
}

runParticleGun pion 211
runParticleGun electron 11
runParticleGun muon 13
runParticleGun photon 22
runParticleGun neutron 2112
runParticleGun jet 1
runParticleGun bjet 5
runParticleGun cjet 4
runParticleGun taujet 15

./Validation $outputroot/particleGun_pion_$cardlabel.root $outputroot/particleGun_electron_$cardlabel.root $outputroot/particleGun_muon_$cardlabel.root $outputroot/particleGun_photon_$cardlabel.root $outputroot/particleGun_neutron_$cardlabel.root $outputroot/particleGun_jet_$cardlabel.root $outputroot/particleGun_bjet_$cardlabel.root $outputroot/particleGun_cjet_$cardlabel.root $outputroot/particleGun_taujet_$cardlabel.root $mainoutputdir/$output $version


# produce calo grid plots 
./CaloGrid $1 ECal
./CaloGrid $1 HCal
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=tmp/tmp1.pdf $outpdf ECal.pdf  
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$outpdf tmp/tmp1.pdf HCal.pdf
mv ECal.pdf $mainoutputdir/www/fig/img_ecal.pdf 
mv ECal.png $mainoutputdir/www/fig/img_ecal.png 
mv HCal.pdf $mainoutputdir/www/fig/img_hcal.pdf
mv HCal.png $mainoutputdir/www/fig/img_hcal.png




