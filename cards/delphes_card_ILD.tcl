# based on arXiv:1306.6329

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  MuonMomentumSmearing

  TrackMerger
  AngularSmearing
  ImpactParameterSmearing

  ECal
  HCal

  ElectronFilter

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
  set Radius 1.8
  # half-length of the magnetic field coverage, in m
  set HalfLength 2.4

  # magnetic field
  set Bz 3.5
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
                                           (abs(eta) <= 2.4)               * (pt > 0.1)    * (0.99) +
                                           (abs(eta) >  2.4)                               * (0.00)}
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
                                           (abs(eta) <= 2.4)               * (pt > 0.1)    * (0.99) +
                                           (abs(eta) >  2.4)                               * (0.00)}
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
                                           (abs(eta) <= 2.4)               * (pt > 0.1)    * (0.99) +
                                           (abs(eta) >  2.4)                               * (0.00)}
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

  # resolution formula for charged hadrons
  set ResolutionFormula {    (abs(eta) <= 1.0)                   * sqrt(0.001^2 + pt^2*1.e-5^2) +
                             (abs(eta) > 1.0 && abs(eta) <= 2.4) * sqrt(0.01^2 + pt^2*1.e-4^2)}


}

###################################
# Momentum resolution for electrons
###################################

module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

   # resolution formula for charged hadrons
  set ResolutionFormula {    (abs(eta) <= 1.0)                   * sqrt(0.001^2 + pt^2*1.e-5^2) +
                             (abs(eta) > 1.0 && abs(eta) <= 2.4) * sqrt(0.01^2 + pt^2*1.e-4^2)}
}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

   # resolution formula for charged hadrons
  set ResolutionFormula {    (abs(eta) <= 1.0)                   * sqrt(0.001^2 + pt^2*1.e-5^2) +
                             (abs(eta) > 1.0 && abs(eta) <= 2.4) * sqrt(0.01^2 + pt^2*1.e-4^2)}

}

##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronMomentumSmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}


########################
# Track angular smearing
########################

module AngularSmearing AngularSmearing {
  set InputArray TrackMerger/tracks
  set OutputArray tracks


  # angular smearing  in eta formula as a function of pt and eta
  set EtaResolutionFormula { 0.001 }

  # angular smearing  in phi formula as a function of pt and eta
  set PhiResolutionFormula { 0.001 }

}

#################################
# Track impact parameter smearing
#################################

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
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true 
 
  set EnergyMin 0.5
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 0.5 degree towers (5x5 mm^2)
  set PhiBins {}
  for {set i -360} {$i <= 360} {incr i} {
    add PhiBins [expr {$i * $pi/360.0}]
  }

  # 0.01 unit in eta up to eta = 2.5
  for {set i -500} {$i <= 500} {incr i} {
    set eta [expr {$i * 0.005}]
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

  set ResolutionFormula { (abs(eta) <= 3.0)                   * sqrt(energy^2*0.01^2 + energy*0.15^2) }

}

#############
#   HCAL
#############

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ECal/eflowTracks

  set TowerOutputArray hcalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowNeutralHadrons

  set IsEcal false 
 
  set EnergyMin 1.0
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower


  # 6 degree towers
  set PhiBins {}
  for {set i -60} {$i <= 60} {incr i} {
    add PhiBins [expr {$i * $pi/60.0}]
  }

  # 0.5 unit in eta up to eta = 3
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

  set ResolutionFormula {                  (abs(eta) <= 3.0) * sqrt(energy^2*0.015^2 + energy*0.50^2)}

}

#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}

###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

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
  add InputArray HCal/eflowTracks
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

#################
# Neutrino Filter
#################

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

  set JetPTMin 20.0
}

############
# Jet finder
############

module FastJetFinder FastJetFinder {
#  set InputArray TowerMerger/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.5

  set JetPTMin 20.0
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


#############
# tau-tagging
#############


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

  add Branch ElectronFilter/electrons Electron Electron
  add Branch MuonMomentumSmearing/muons Muon Muon
  add Branch JetEnergyScale/jets Jet Jet
  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}

