#! /bin/sh
function runParticleGun {
  name=$1
  pid=$2
  cmnd="examples/Pythia8/configParticleGun_$name.cmnd"
  sed '/Main:spareMode1/s/=[[:space:]]*[0-9]*/= '$pid'/' examples/Pythia8/configParticleGun.cmnd > $cmnd
  ./DelphesPythia8 cards/delphes_card_CMS.tcl $cmnd delphes_ParticleGun_$name.root
}

runParticleGun electron 11
runParticleGun muon 13
runParticleGun photon 22
runParticleGun jet 1
runParticleGun bjet 5
runParticleGun taujet 15

./Validation delphes_ParticleGun_electron.root delphes_ParticleGun_muon.root delphes_ParticleGun_photon.root delphes_ParticleGun_jet.root delphes_ParticleGun_bjet.root delphes_ParticleGun_taujet.root delphes_validation.root
