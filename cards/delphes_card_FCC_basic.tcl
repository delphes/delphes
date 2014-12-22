#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing

  TrackMerger
  AngularSmearing
  ImpactParameterSmearing

  ECal
  HCal

  TowerMerger
  EFlowMerger

  MissingET

  NeutrinoFilter
  GenJetFinder
  FastJetFinder

  JetEnergyScale

  TrackCountingBTagging
  TauTagging

  ScalarHT

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
  set Radius 2.60
  # half-length of the magnetic field coverage, in m
  set HalfLength 6.00

  # magnetic field
  set Bz 6.0
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

  # add EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for charged hadrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.75) +
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.95) +
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.60) +
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 1.0)                  * (0.90) +
                         (abs(eta) > 4.0)                                                  * (0.00)}

}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for electrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.75) +
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.99) +
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.70) +
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 1.0)                  * (0.98) +
                         (abs(eta) > 4.0)                                                  * (0.00)}
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for muons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.75) +
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.99) +
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.70) +
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 1.0)                  * (0.98) +
                         (abs(eta) > 4.0)                                                  * (0.00)}
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for charged hadrons
  set ResolutionFormula {    (abs(eta) <= 1.5)                   * (pt > 0.1) * (0.01 + pt*2.e-5) +
                             (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 0.1) * (0.02 + pt*3.e-5)}


}

#################################
# Energy resolution for electrons
#################################

module EnergySmearing ElectronEnergySmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

  # resolution formula for electrons
  set ResolutionFormula {                  (abs(eta) <= 4.0) * (energy > 0.1   && energy <= 2.0e1) * (energy*0.007) +
                                           (abs(eta) <= 4.0) * (energy > 2.0e1)                    * sqrt(energy^2*0.005^2 + energy*0.02^2) +
                                           (abs(eta) > 4.0 && abs(eta) <= 6.0)                     * sqrt(energy^2*0.05^2 + energy*1.00^2)}

}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for muons
  set ResolutionFormula {    (abs(eta) <= 1.5)                   * (pt > 0.1) * (0.01 + pt*5.e-6) +
                             (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 0.1) * (0.02 + pt*1.e-5)}


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
# Track angular smearing
################################

module AngularSmearing AngularSmearing {
  set InputArray TrackMerger/tracks
  set OutputArray tracks


  # angular smearing  in eta formula as a function of pt and eta
  set EtaResolutionFormula { 0.001 }

  # angular smearing  in phi formula as a function of pt and eta
  set PhiResolutionFormula { 0.001 }

}

################################
# Track impact parameter smearing
################################

module ImpactParameterSmearing ImpactParameterSmearing {
  set InputArray AngularSmearing/tracks
  set OutputArray tracks


  # absolute impact parameter smearing formula (in mm) as a function of pt and eta
  set ResolutionFormula {(pt > 0.1  && pt <= 5.0)   * (0.010) +
                         (pt > 5.0)                 * (0.005)}

}

#############
#   ECAL
#############

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ImpactParameterSmearing/tracks

  set TowerOutputArray ecalTowers
  set EFlowTowerOutputArray eflowPhotons

  set EnergyMin 0.5
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 0.5 degree towers
  set PhiBins {}
  for {set i -360} {$i <= 360} {incr i} {
    add PhiBins [expr {$i * $pi/360.0}]
  }

  # 0.01 unit in eta up to eta = 3
  for {set i -300} {$i <= 300} {incr i} {
    set eta [expr {$i * 0.01}]
    add EtaPhiBins $eta $PhiBins
  }

  # default energy fractions {abs(PDG code)} {fraction of energy deposited in ECAL}

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

  # set ECalResolutionFormula {resolution formula as a function of eta and energy}

  # This is the CMS ECAL resolution, extended up eta = 6.0
  set ResolutionFormula { (abs(eta) <= 3.0)                   * sqrt(energy^2*0.003^2 + energy*0.029^2 + 0.125^2)  +
                 	  (abs(eta) > 3.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.107^2 + energy*2.08^2)}


}

#############
#   HCAL
#############

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ImpactParameterSmearing/tracks

  set TowerOutputArray hcalTowers
  set EFlowTowerOutputArray eflowNeutralHadrons

  set EnergyMin 1.0
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower


  # 5 degree towers
  set PhiBins {}
  for {set i -72} {$i <= 72} {incr i} {
    add PhiBins [expr {$i * $pi/75.0}]
  }

  # 0.05 unit in eta up to eta = 3
  for {set i -60} {$i <= 60} {incr i} {
    set eta [expr {$i * 0.05}]
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

  # set HCalResolutionFormula {resolution formula as a function of eta and energy}

  # This is the ATLAS HCAL resolution, extended up eta = 6.0
  set HCalResolutionFormula {                  (abs(eta) <= 1.7) * sqrt(energy^2*0.0302^2 + energy*0.5205^2 + 1.59^2) +
                             (abs(eta) > 1.7 && abs(eta) <= 3.2) * sqrt(energy^2*0.0500^2 + energy*0.706^2) +
                             (abs(eta) > 3.2 && abs(eta) <= 6.0) * sqrt(energy^2*0.09420^2 + energy*1.00^2)}

}

####################
# Tower Merger (in case not using e-flow algorithm)
####################

module Merger TowerMerger {
# add InputArray InputArray
  add InputArray ECal/ecalTowers
  add InputArray HCal/hcalTowers
  set OutputArray towers
}

####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray ImpactParameterSmearing/tracks
  add InputArray ECal/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflow
}


###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}


##################
# Scalar HT merger
##################

module Merger ScalarHT {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set EnergyOutputArray energy
}

#####################
# Neutrino Filter
#####################

module PdgCodeFilter NeutrinoFilter {

  set InputArray Delphes/stableParticles
  set OutputArray filteredParticles

  set PTMin 0.0

  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}

}


#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.5

  set JetPTMin 100.0
}

############
# Jet finder
############

module FastJetFinder FastJetFinder {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.5

  set JetPTMin 100.0
}

##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.00}
}

##########################
# Track Counting b-tagging
##########################

module TrackCountingBTagging TrackCountingBTagging {
  set TrackInputArray ImpactParameterSmearing/tracks
  set JetInputArray JetEnergyScale/jets

  set BitNumber 0

  # maximum distance between jet and track
  set DeltaR 0.3

  # minimum pt of tracks
  set TrackPTMin 1.0

  # minimum transverse impact parameter (in mm)
  set TrackIPMax 2.0

  # minimum ip significance for the track to be counted
  set SigMin 6.5

  # minimum number of tracks (high efficiency n=2, high purity n=3)
  set Ntracks 3
}


##########################
# tau-tagging
##########################


module TauTagging TauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.4}
}

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle
  add Branch GenJetFinder/jets GenJet Jet

  add Branch ChargedHadronMomentumSmearing/chargedHadrons ChargedHadron Track
  add Branch HCal/eflowNeutralHadrons NeutralHadron Tower
  add Branch ECal/eflowPhotons Photon Photon

  add Branch ElectronEnergySmearing/electrons Electron Electron
  add Branch MuonMomentumSmearing/muons Muon Muon
  add Branch JetEnergyScale/jets Jet Jet
  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}
