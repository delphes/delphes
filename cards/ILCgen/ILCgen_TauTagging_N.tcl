#####################################
# tau-tagging for inclusive jets, N=2
#####################################
module TauTagging TauTagging_N2 {
    set ParticleInputArray Delphes/allParticles
    set PartonInputArray Delphes/partons
    set JetInputArray JetFinder_N2/jets

    source ILCgen/ILCgen_TauTagging.tcl
}
#####################################
# tau-tagging for inclusive jets, N=3
#####################################
module TauTagging TauTagging_N3 {
    set ParticleInputArray Delphes/allParticles
    set PartonInputArray Delphes/partons
    set JetInputArray JetFinder_N3/jets

    source ILCgen/ILCgen_TauTagging.tcl
}
#####################################
# tau-tagging for inclusive jets, N=4
#####################################
module TauTagging TauTagging_N4 {
    set ParticleInputArray Delphes/allParticles
    set PartonInputArray Delphes/partons
    set JetInputArray JetFinder_N4/jets

    source ILCgen/ILCgen_TauTagging.tcl
}
#####################################
# tau-tagging for inclusive jets, N=5
#####################################
module TauTagging TauTagging_N5 {
    set ParticleInputArray Delphes/allParticles
    set PartonInputArray Delphes/partons
    set JetInputArray JetFinder_N5/jets

    source ILCgen/ILCgen_TauTagging.tcl
}
#####################################
# tau-tagging for inclusive jets, N=6
#####################################
module TauTagging TauTagging_N6 {
    set ParticleInputArray Delphes/allParticles
    set PartonInputArray Delphes/partons
    set JetInputArray JetFinder_N6/jets

    source ILCgen/ILCgen_TauTagging.tcl
}
