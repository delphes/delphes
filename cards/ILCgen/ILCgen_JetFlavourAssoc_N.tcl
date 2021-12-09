################################################
# Jet flavor association for inclusive jets, N=2
################################################
module JetFlavorAssociation JetFlavor_N2 {
    set PartonInputArray Delphes/partons
    set ParticleInputArray Delphes/allParticles
    set ParticleLHEFInputArray Delphes/allParticlesLHEF
    set JetInputArray JetFinder_N2/jets

    set DeltaR 0.5
    set PartonPTMin 1.0
    set PartonEtaMax 4.0
}
################################################
# Jet flavor association for inclusive jets, N=3
################################################
module JetFlavorAssociation JetFlavor_N3 {
    set PartonInputArray Delphes/partons
    set ParticleInputArray Delphes/allParticles
    set ParticleLHEFInputArray Delphes/allParticlesLHEF
    set JetInputArray JetFinder_N3/jets

    set DeltaR 0.5
    set PartonPTMin 1.0
    set PartonEtaMax 4.0
}
################################################
# Jet flavor association for inclusive jets, N=4
################################################
module JetFlavorAssociation JetFlavor_N4 {
    set PartonInputArray Delphes/partons
    set ParticleInputArray Delphes/allParticles
    set ParticleLHEFInputArray Delphes/allParticlesLHEF
    set JetInputArray JetFinder_N4/jets

    set DeltaR 0.5
    set PartonPTMin 1.0
    set PartonEtaMax 4.0
}
################################################
# Jet flavor association for inclusive jets, N=5
################################################
module JetFlavorAssociation JetFlavor_N5 {
    set PartonInputArray Delphes/partons
    set ParticleInputArray Delphes/allParticles
    set ParticleLHEFInputArray Delphes/allParticlesLHEF
    set JetInputArray JetFinder_N5/jets

    set DeltaR 0.5
    set PartonPTMin 1.0
    set PartonEtaMax 4.0
}
################################################
# Jet flavor association for inclusive jets, N=6
################################################
module JetFlavorAssociation JetFlavor_N6 {
    set PartonInputArray Delphes/partons
    set ParticleInputArray Delphes/allParticles
    set ParticleLHEFInputArray Delphes/allParticlesLHEF
    set JetInputArray JetFinder_N6/jets

    set DeltaR 0.5
    set PartonPTMin 1.0
    set PartonEtaMax 4.0
}
