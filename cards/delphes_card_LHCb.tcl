#set MaxEvents 1000
#set RandomSeed 123


#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  ParticlePropagator

  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing

  TrackMerger
  ImpactParameterSmearing

  IdentificationMap

  ECal
  HCal

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

  set Radius 3.31


  # half-length of the magnetic field coverage, in m
  set HalfLength 12.0

  # magnetic field
  set Bz 1.1

  # Need to veto anything with theta > 0.269 rad  -> eta = 2
  #                            theta < 0.0135 rad -> eta = 5

  # tracker and calos are at approx 0.269 rad, R = 12*tan(0.269)

}


########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for charged hadrons
  set ResolutionFormula {(eta > 2.0  && eta <= 5.0)      * (pt > 0.5) * (0.005)}
}

#################################
# Energy resolution for electrons
#################################

module EnergySmearing ElectronEnergySmearing {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

  # resolution formula for electrons
  set ResolutionFormula { (eta > 2.0  && eta <= 5.0) * (energy > 0.1   && energy <= 8.0) * (energy*0.05) +
                          (eta > 2.0  && eta <= 5.0) * (energy > 8.0)                    *  sqrt(energy^2*0.015^2 + energy*0.10^2)}
  }

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  # resolution formula for muons
  set ResolutionFormula {(eta > 2.0  && eta <= 5.0)      * (pt > 0.5)* (0.005)}
}


##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronEnergySmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}



################################
# Track impact parameter smearing
################################

module ImpactParameterSmearing ImpactParameterSmearing {
  set InputArray TrackMerger/tracks
  set OutputArray tracks


# absolute impact parameter smearing formula (in mm) as a function of pt and eta
set ResolutionFormula {0.0116 + 0.0234/pt}

}

###### ref. JINST 9 C01065 #############


####################################
# Charged hadron PID
####################################

module IdentificationMap IdentificationMap {
  set InputArray ImpactParameterSmearing/tracks
  set OutputArray tracks

  # {PID in} {PID out} {formula}
  # make sure "PID in" and "PID out" have the same charge (e.g {-13} {211} or {-321} {211})
  # {211} {-13} is equivalent to {-211} {13} (and needs to be written once only...)






  # --- pions ---

  add EfficiencyFormula {211} {211} {      (eta <= 2.0)                                  * (0.00) +
                                           (eta > 2.0  && eta <= 5.0) *       (pt < 0.8) * (0.00) +
					   (eta > 2.0  && eta <= 5.0) *       (pt >= 0.8)* (0.95) +
					   (eta > 5.0)                                   * (0.00)}

  add EfficiencyFormula {211} {-13} {      (eta <= 2.0)                                 * (0.00) +
                                           (eta > 2.0  && eta <= 5.0) *      (pt < 0.8) * (0.00) +
					   (eta > 2.0  && eta <= 5.0) *      (pt >= 0.8)* (0.005 + 0.0663*exp(-0.13*pt*cosh(eta))) +
					   (eta > 5.0)                                  * (0.00)}


 # --- kaons ---


  add EfficiencyFormula {321} {321} {      (eta <= 2.0)                                  * (0.00) +
                                           (eta > 2.0  && eta <= 5.0) *       (pt < 0.8) * (0.00) +
					   (eta > 2.0  && eta <= 5.0) *       (pt >= 0.8)* (0.95) +
					   (eta > 5.0)                                   * (0.00)}

  add EfficiencyFormula {321} {-13} {      (eta <= 2.0)                                 * (0.00) +
                                           (eta > 2.0  && eta <= 5.0) *      (pt < 0.8) * (0.00) +
					   (eta > 2.0  && eta <= 5.0) *      (pt >= 0.8)* (0.005 + 0.086*exp(-0.11*pt*cosh(eta))) +
					   (eta > 5.0)                                  * (0.00)}


 # --- protons ---


  add EfficiencyFormula {2212} {2212} {    (eta <= 2.0)                                  * (0.00) +
                                           (eta > 2.0  && eta <= 5.0) *       (pt < 0.8) * (0.00) +
					   (eta > 2.0  && eta <= 5.0) *       (pt >= 0.8)* (0.95) +
					   (eta > 5.0)                                   * (0.00)}

  add EfficiencyFormula {2212} {-13} {     (eta <= 2.0)                                 * (0.00) +
                                           (eta > 2.0  && eta <= 5.0) *      (pt < 0.8) * (0.00) +
					   (eta > 2.0  && eta <= 5.0) *      (pt >= 0.8)* (0.002) +
					   (eta > 5.0)                                  * (0.00)}



 # --- muons ---



  add EfficiencyFormula {-13} {-13} {      (eta <= 2.0)                                * (0.00) +
                                           (eta > 2.0  && eta <= 5.0) *      (pt < 0.8)* (0.00) +
                                           (eta > 2.0  && eta <= 5.0) *     (pt >= 0.8)* (0.97) +
                                           (eta > 5.0)                                 * (0.00)}


 # efficiency for other charged particles (should be always {0} {0} {formula})

  add EfficiencyFormula {0} {0}     {      (eta <= 2.0)                                * (0.00) +
                                           (eta > 2.0  && eta <= 5.0) *     (pt < 0.8) * (0.00) +
					   (eta > 2.0  && eta <= 5.0) *     (pt > 0.8) * (0.95) +
                                           (eta > 5.0)                                 * (0.00)}

}

###### ref. JINST 8 P10020 #############


#############
#   ECAL
#############

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray IdentificationMap/tracks

  set TowerOutputArray ecalTowers
  set EFlowTowerOutputArray eflowPhotons

  set EnergyMin 0.0
  set EnergySignificanceMin 0.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

   # 1 degree towers
  set PhiBins {}
  for {set i -180} {$i <= 180} {incr i} {
    add PhiBins [expr {$i * $pi/180.0}]
  }

  # 0.02 unit in eta  from eta = 3.2 to eta = 5.0
  for {set i 1} {$i <= 90} {incr i} {
    set eta [expr {3.2 + $i * 0.02}]
    add EtaPhiBins $eta $PhiBins
  }

    # 1.25 degree towers
  set PhiBins {}
  for {set i -135} {$i <= 135} {incr i} {
    add PhiBins [expr {$i * $pi/135.0}]
  }

  # 0.025 unit in eta  from eta = 2.6 to eta = 3.2
  for {set i 1} {$i <= 24} {incr i} {
    set eta [expr {2.6 + $i * 0.025}]
    add EtaPhiBins $eta $PhiBins
  }

     # 1.25 degree towers
  set PhiBins {}
  for {set i -100} {$i <= 100} {incr i} {
    add PhiBins [expr {$i * $pi/100.0}]
  }

  # 0.04 unit in eta  from eta = 2.0 to eta = 2.6
  for {set i 0} {$i <= 24} {incr i} {
    set eta [expr {2.0 + $i * 0.04}]
    add EtaPhiBins $eta $PhiBins
  }

  add EnergyFraction {0} {0.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0}
  add EnergyFraction {22} {1.0}
  add EnergyFraction {111} {1.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3}
  add EnergyFraction {3122} {0.3}

  set ResolutionFormula {(eta <= 5.0 && eta > 2.0) * sqrt(energy^2*0.015^2 + energy*0.10^2)}
}

#############
#   HCAL
#############

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray IdentificationMap/tracks

  set TowerOutputArray hcalTowers
  set EFlowTowerOutputArray eflowNeutralHadrons

  set EnergyMin 0.0
  set EnergySignificanceMin 0.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

    # 1 degree towers
  set PhiBins {}
  for {set i -16} {$i <= 16} {incr i} {
    add PhiBins [expr {$i * $pi/16.0}]
  }

  # 0.20 unit in eta  from eta = 2.6 to eta = 5.0
  for {set i 1} {$i <= 12} {incr i} {
    set eta [expr {2.6 + $i * 0.2}]
    add EtaPhiBins $eta $PhiBins
  }

    # 1 degree towers
  set PhiBins {}
  for {set i -32} {$i <= 32} {incr i} {
    add PhiBins [expr {$i * $pi/32.0}]
  }

  # 0.1 unit in eta  from eta = 2 to eta = 2.6
  for {set i 0} {$i <= 6} {incr i} {
    set eta [expr {2.0 + $i * 0.1}]
    add EtaPhiBins $eta $PhiBins
  }



  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
  add EnergyFraction {0} {1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {0.0}
  add EnergyFraction {22} {0.0}
  add EnergyFraction {111} {0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.7}
  add EnergyFraction {3122} {0.7}

  set ResolutionFormula { (eta <= 5.0 && eta > 2.0) * sqrt(energy^2*0.05^2 + energy*0.80^2)}
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

  add Branch IdentificationMap/tracks Track Track
  add Branch HCal/eflowNeutralHadrons NeutralHadron Tower
  add Branch ECal/eflowPhotons Photon Photon

}

