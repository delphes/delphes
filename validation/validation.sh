#!/bin/bash
################################################################################
#
# This code produces a set of validation plots for a given detector card.
#
# In order to run this you need to compile Delphes with Pythia8 first, see:
#
# https://cp3.irmp.ucl.ac.be/projects/delphes/wiki/WorkBook/Pythia8
#
# After you (re-)compiled Delphes with Pythia8 you are ready to go, execute from Delphes main dir:
#
# ./validation/validation.sh [detector_card] [number_of_events]
#
#  e.g.
#
# ./validation/validation.sh cards/delphes_card_CMS.tcl 100000
#
# Note that the more events you specify, the more accurate the controls plots will be ...
# This said, 500k events should be ok for most cases.
#
################################################################################

EXPECTED_ARGS=2
E_BADARGS=65

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: ./validation/validation.sh [detector_card] [number_of_events]"
  echo "for instance: ./validation/validation.sh cards/delphes_card_CMS.tcl 10000"
  exit $E_BADARGS
fi

cardbase=$(basename $1)
carddir=$(dirname $1)
nEvents=$2
output=validation_${cardbase%.*}.root
mainoutputdir=report_${cardbase%.*}
outputrootdir=report_${cardbase%.*}/root
cardlabel=${cardbase%.*}
version=x.y.z
outpdf=$mainoutputdir/${output%.*}.pdf
cardsdir=validation/cards
samplesdir=validation/samples
validationcard=$carddir/validation_$cardbase

mkdir -p $cardsdir
mkdir -p $samplesdir
mkdir -p $outputrootdir
mkdir -p $mainoutputdir/www/fig

sed "s/delphes_card_CMS.tcl/$cardbase/g" cards/validation_card.tcl > $validationcard
sed -i "1i set MaxEvents $nEvents" $validationcard

function runParticleGun {
  name=$1
  pid=$2
  cmnd=$cardsdir/configParticleGun_$name.cmnd
  outputroot=particleGun_${name}_${cardlabel}.root
  sed "/Main:spareMode1/s|=[[:space:]]*[0-9]*|= $pid|; /Main:numberOfEvents/s|=[[:space:]]*[0-9]*|= $nEvents|" examples/Pythia8/configParticleGun.cmnd > $cmnd
  ./DelphesPythia8 $validationcard $cmnd $outputrootdir/$outputroot

}


function runJetsGun {
  name=$1
  pid=$2
  lhe=$samplesdir/events_$name.lhe
  cmnd=$cardsdir/configLHE_$name.cmnd
  inputroot=$samplesdir/$name.root
  outputroot=particleGun_${name}_${cardlabel}.root

  if [ ! -f $inputroot ]
  then
    python validation/flatGunLHEventProducer.py --pdg $pid  --guntype pt --pmin 1 --pmax 50000 --etamin -6 --etamax 6 --nevts $nEvents --seed 1 --output $lhe --log --ecm 100000
    gunzip $lhe.gz
    sed "/Beams:LHEF/s|=[[:space:]].*|= $lhe|; /Main:numberOfEvents/s/=[[:space:]]*[0-9]*/= $nEvents/" examples/Pythia8/configLHE.cmnd > $cmnd

    ./DelphesPythia8 cards/gen_card.tcl $cmnd $inputroot
  fi

  ./DelphesROOT $validationcard $outputrootdir/$outputroot $inputroot
}


runParticleGun pion 211 &
runParticleGun electron 11 &
runParticleGun muon 13 &
runParticleGun photon 22 &
runParticleGun neutron 2112 &
runParticleGun taujet 15 &
runJetsGun jet 1 &
runJetsGun bjet 5 &
runJetsGun cjet 4 &

wait
echo all particle guns complete ...

./DelphesValidation $outputrootdir/particleGun_pion_$cardlabel.root $outputrootdir/particleGun_electron_$cardlabel.root $outputrootdir/particleGun_muon_$cardlabel.root $outputrootdir/particleGun_photon_$cardlabel.root $outputrootdir/particleGun_neutron_$cardlabel.root $outputrootdir/particleGun_jet_$cardlabel.root $outputrootdir/particleGun_bjet_$cardlabel.root $outputrootdir/particleGun_cjet_$cardlabel.root $outputrootdir/particleGun_taujet_$cardlabel.root $mainoutputdir/$output $version


# produce calo grid plots
./CaloGrid $1 ECal
./CaloGrid $1 HCal
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=tmp/tmp1.pdf $outpdf ECal.pdf
gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=$outpdf tmp/tmp1.pdf HCal.pdf
mv ECal.pdf $mainoutputdir/www/fig/img_ecal.pdf
mv ECal.png $mainoutputdir/www/fig/img_ecal.png
mv HCal.pdf $mainoutputdir/www/fig/img_hcal.pdf
mv HCal.png $mainoutputdir/www/fig/img_hcal.png

