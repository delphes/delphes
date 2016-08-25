#!/bin/bash
function runParticleGun {
    in=$1
    ID=$2
    outcmnd="examples/Pythia8/configParticleGun$in.cmnd"
    sed 's/ID/'$ID'/g' examples/Pythia8/configParticleGun.cmnd > $outcmnd
    ./DelphesPythia8 cards/delphes_card_CMS.tcl $outcmnd delphes_ParticleGun$in.root
}
runParticleGun el 11
runParticleGun mu 13
runParticleGun ph 22
runParticleGun b 5
runParticleGun jet 1
runParticleGun tau 15

./Validation delphes_ParticleGunel.root delphes_ParticleGunmu.root delphes_ParticleGunph.root delphes_ParticleGunjet.root delphes_ParticleGunb.root delphes_ParticleGuntau.root validation.root

