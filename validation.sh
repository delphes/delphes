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
runParticleGun b 5
runParticleGun jet 1
hadd delphes_all.root delphes_ParticleGunmu.root delphes_ParticleGunb.root delphes_ParticleGunjet.root



