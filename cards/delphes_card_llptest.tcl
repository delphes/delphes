set MaxEvents 100

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  ParticlePropagator

  ElectronTrackingEfficiency
  ElectronMomentumSmearing
  TreeWriter
}

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.00

  # magnetic field
  set Bz 3.8
}


##############################
# Electron tracking efficiency
##############################

module TrackEfficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons


  # tracking efficiency formula (efficiency formula as a function of  pt, eta,  phi,  energy, d0 (mm), dz (mm), ctgTheta)
  set EfficiencyFormula {
                          (d0 < 1e-3) *                                                     (pt <= 0.1)   * (0.00) +
                          (d0 < 1e-3) *                (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)      * (0.73) +
                          (d0 < 1e-3) *                (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e2)    * (0.95) +
                          (d0 < 1e-3) *                 (abs(eta) <= 1.5) * (pt > 1.0e2)                  * (0.99) +
                          (d0 < 1e-3) * (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.50) +
                          (d0 < 1e-3) * (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e2) * (0.83) +
                          (d0 < 1e-3) * (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e2)                * (0.90) +
                          (d0 < 1e-3) * (abs(eta) > 2.5)                                                  * (0.00) +
                          (d0 > 1e-3) *                                                                   * (0.00)
                         }
}

###################################
# Momentum resolution for electrons
###################################

module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

  # resolution formula for electrons
  # based on arXiv:1502.02701
  set ResolutionFormula {                  (abs(eta) <= 0.5) * (pt > 0.1) * sqrt(0.03^2 + pt^2*1.3e-3^2) +
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 0.1) * sqrt(0.05^2 + pt^2*1.7e-3^2) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1) * sqrt(0.15^2 + pt^2*3.1e-3^2)}
}

##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle
  add Branch ElectronMomentumSmearing/electrons Electron Electron
}
