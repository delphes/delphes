#
# Official Delphes card prepared by FCC-hh collaboration
#
#  Main authors:  Michele Selvaggi (CERN)
#
#  Released on: October 14th, 2020
#
#  - fix muon resolution at high pT
#  - updated btagging, tau tagging and photon ID
#
#
#  Configuration: FCC-hh baseline detector
#
#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  
  ParticlePropagator
  TrackMergerProp

  DenseProp
  DenseMergeTracks
  DenseTrackFilter

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  MuonMomentumSmearing

  TrackMerger

  ECal
  HCal

  Calorimeter
  EFlowMerger
  EFlowFilter
  
  PhotonEfficiency
  PhotonIsolation

  ElectronFilter
  ElectronIsolation

  ChargedHadronFilter

  MuonIsolation

  NeutrinoFilter

  MissingET
  GenMissingET

  GenJetFinder02
  GenJetFinder04
  GenJetFinder08
  GenJetFinder15

  FastJetFinder02
  FastJetFinder04
  FastJetFinder08
  FastJetFinder15

  CaloJetFinder02
  CaloJetFinder04
  CaloJetFinder08
  CaloJetFinder15

  TrackJetFinder02
  TrackJetFinder04
  TrackJetFinder08
  TrackJetFinder15

  JetEnergyScale

  JetFlavorAssociation

  BTagging
  CTagging
  TauTagging

  ScalarHT

  UniqueObjectFinder

  TreeWriter
}

#####################################
# Track propagation to calorimeters
#####################################

module ParticlePropagator ParticlePropagator {
  set InputArray Delphes/stableParticles
  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.5
  # half-length of the magnetic field coverage, in m
  set HalfLength 5

  # magnetic field
  set Bz 4.0
}


##############
# Track merger
##############

module Merger TrackMergerProp {
# add InputArray InputArray
  add InputArray ParticlePropagator/chargedHadrons
  add InputArray ParticlePropagator/electrons
  add InputArray ParticlePropagator/muons
  set OutputArray tracks
}


####################################
# Track propagation to pseudo-pixel
####################################

module ParticlePropagator DenseProp {

  set InputArray TrackMergerProp/tracks

  # radius of the magnetic field coverage, in m
  set Radius 0.45
  set RadiusMax 1.5
  # half-length of the magnetic field coverage, in m
  set HalfLength 0.8
  set HalfLengthMax 5

  # magnetic field
  set Bz 4.0
}

#####################
# Dense Track merger
#####################

module Merger DenseMergeTracks {
# add InputArray InputArray
  add InputArray DenseProp/chargedHadrons
  add InputArray DenseProp/electrons
  add InputArray DenseProp/muons
  set OutputArray tracks
}


######################
#   Dense Track Filter
######################

module DenseTrackFilter DenseTrackFilter {

  set TrackInputArray DenseMergeTracks/tracks

  set TrackOutputArray tracks
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  set EtaPhiRes 0.001
  set EtaMax 6.0

  set pi [expr {acos(-1)}]

  set nbins_phi [expr {$pi/$EtaPhiRes} ]
  set nbins_phi [expr {int($nbins_phi)} ]

  set PhiBins {}
  for {set i -$nbins_phi} {$i <= $nbins_phi} {incr i} {
    add PhiBins [expr {$i * $pi/$nbins_phi}]
  }

  set nbins_eta [expr {$EtaMax/$EtaPhiRes} ]
  set nbins_eta [expr {int($nbins_eta)} ]

  for {set i -$nbins_eta} {$i <= $nbins_eta} {incr i} {
    set eta [expr {$i * $EtaPhiRes}]
    add EtaPhiBins $eta $PhiBins
  }
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray DenseTrackFilter/chargedHadrons
  set OutputArray chargedHadrons

 # TBC (which eta_max ? which pT min?)

 # tracking efficiency formula for charged hadrons

  set EfficiencyFormula { (pt <= 0.5) * (0.00) + \
(abs(eta) <= 2.5) * (pt > 0.5 && pt <= 1) * (0.90) + \
(abs(eta) <= 2.5) * (pt > 1) * (0.95) + \
(abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.5 && pt <= 1) * (0.85) + \
(abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 1) * (0.90) + \
(abs(eta) > 4 && abs(eta) <= 6) * (pt > 0.5 && pt <= 1) * (0.80) + \
(abs(eta) > 4 && abs(eta) <= 6) * (pt > 1.0) * (0.85) + \
(abs(eta) > 6.0) * (0.00)}

}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray DenseTrackFilter/electrons
  set OutputArray electrons

# TBC (which eta_max ?)
# putting same as charged hadrons for now...

  set EfficiencyFormula { (pt <= 0.5) * (0.00) + \
  (abs(eta) <= 2.5) * (pt > 0.5 && pt <= 1) * (0.90) + \
  (abs(eta) <= 2.5) * (pt > 1) * (0.95) + \
  (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.5 && pt <= 1) * (0.85) + \
  (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 1) * (0.90) + \
  (abs(eta) > 4 && abs(eta) <= 6) * (pt > 0.5 && pt <= 1) * (0.80) + \
  (abs(eta) > 4 && abs(eta) <= 6) * (pt > 1.0) * (0.85) + \
  (abs(eta) > 6.0) * (0.00)}

}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray DenseTrackFilter/muons
  set OutputArray muons

# TBC (which eta_max ? why eff = 0 for 4 < eta < 6 ? for now put the same as central)
# what about high pT ?
 # tracking efficiency formula for muons
  set EfficiencyFormula { (pt <= 0.5) * (0.00) + \
  (abs(eta) <= 6.0) * (pt > 0.5 && pt <= 1) * (0.90) + \
  (abs(eta) <= 6.0) * (pt > 1) * (0.99) + \
  (abs(eta) > 6.0) * (0.00)}

}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  source momentumResolutionVsP.tcl 
}


###################################
# Momentum resolution for electrons
###################################

module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  source momentumResolutionVsP.tcl 
}


###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # TBC for just putting tracker resolution/ need to add improvement at high pT

  source muonMomentumResolutionVsP.tcl 
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


#############
#   ECAL
#############

# TBC : calos seems ok, check eta max value though.

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true

  set EnergyMin 0.5
  set EnergySignificanceMin 2.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 0.012 rad towers up to eta = 2.5 (barrel)
   set PhiBins {}
  for {set i -256} {$i <= 256} {incr i} {
    add PhiBins [expr {$i * $pi/256.0}]
  } 

  # 0.01 unit in eta up to eta = 2.5 (barrel)
  for {set i -249} {$i <= 250} {incr i} {
    set eta [expr {$i * 0.01}]
    add EtaPhiBins $eta $PhiBins
  }

  # 0.025 rad between 2.5 and 6.0
  set PhiBins {}
  for {set i -128} {$i <= 128} {incr i} {
    add PhiBins [expr {$i * $pi/128.0}]
  }

  # 0.025 unit in eta between 2.5 and 6.0
  for {set i 0} {$i <= 140} {incr i} {
    set eta [expr { -6.00 + $i * 0.025}]
    add EtaPhiBins $eta $PhiBins
  }

  for {set i 0} {$i <= 139} {incr i} {
    set eta [expr { 2.525 + $i * 0.025}]
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
  # add EnergyFraction {310} {0.3}
  # add EnergyFraction {3122} {0.3}

  # set ECalResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {                     (abs(eta) <= 4.0) * sqrt(energy^2*0.007^2 + energy*0.10^2) + \
                            (abs(eta) > 4.0 && abs(eta) <= 6.0) * sqrt(energy^2*0.035^2  + energy*0.30^2)}


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
  set EnergySignificanceMin 2.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # 0.025 rad towers up to eta = 2.5 (barrel)
   set PhiBins {}
  for {set i -128} {$i <= 128} {incr i} {
    add PhiBins [expr {$i * $pi/128.0}]
  } 

  # 0.025 unit in eta up to eta = 2.5 (barrel)
  for {set i -99} {$i <= 100} {incr i} {
    set eta [expr {$i * 0.025}]
    add EtaPhiBins $eta $PhiBins
  }

  # 0.05 rad between 2.5 and 6.0
  set PhiBins {}
  for {set i -64} {$i <= 64} {incr i} {
    add PhiBins [expr {$i * $pi/64.0}]
  }

  # 0.05 unit in eta between 2.5 and 6.0
  for {set i 0} {$i <= 70} {incr i} {
    set eta [expr { -6.00 + $i * 0.05}]
    add EtaPhiBins $eta $PhiBins
  }

  for {set i 0} {$i <= 69} {incr i} {
    set eta [expr { 2.55 + $i * 0.05}]
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
  # add EnergyFraction {310} {0.7}
  # add EnergyFraction {3122} {0.7}

   # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {                     (abs(eta) <= 1.7) * sqrt(energy^2*0.03^2 + energy*0.50^2) + \
                            (abs(eta) > 1.7 && abs(eta) <= 4.0) * sqrt(energy^2*0.03^2 + energy*0.60^2) + \
                            (abs(eta) > 4.0 && abs(eta) <= 6.0) * sqrt(energy^2*0.10^2 + energy*1.00^2)}
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


######################
# ChargedHadronFilter
######################

module PdgCodeFilter ChargedHadronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray chargedHadrons
  
  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}


###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

module Merger Calorimeter {
# add InputArray InputArray
  add InputArray ECal/ecalTowers
  add InputArray HCal/hcalTowers
  add InputArray MuonMomentumSmearing/muons
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

######################
# EFlowFilter
######################

module PdgCodeFilter EFlowFilter {
  set InputArray EFlowMerger/eflow
  set OutputArray eflow
  
  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
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


#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {

# add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}



#####################
# MC truth jet finder
#####################

# TBC: is jet radius fine?

module FastJetFinder GenJetFinder02 {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.2

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.2

  set JetPTMin 25.0
}


#####################
# MC truth jet finder
#####################

# TBC: is jet radius fine?

module FastJetFinder GenJetFinder04 {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.4

  set JetPTMin 25.0
}
#####################
# MC truth jet finder
#####################

# TBC: is jet radius fine?

module FastJetFinder GenJetFinder08 {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.8

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

  set JetPTMin 25.0
}

#####################
# MC truth jet finder
#####################

# TBC: is jet radius fine?

module FastJetFinder GenJetFinder15 {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.5

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 1.5

  set JetPTMin 25.0
}


##################
# Fast Jet finder
##################

module FastJetFinder FastJetFinder02 {
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.2

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.2

  set JetPTMin 25.0
}

##################
# Fast Jet finder
##################

module FastJetFinder FastJetFinder04 {
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.4

  set JetPTMin 25.0
}


##################
# Fast Jet finder
##################

module FastJetFinder FastJetFinder08 {
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.8

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

  set JetPTMin 25.0
}

##################
# Fast Jet finder
##################

module FastJetFinder FastJetFinder15 {
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.5

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 1.5

  set JetPTMin 25.0
}


##################
# Fast Jet finder
##################

module FastJetFinder CaloJetFinder02 {
  set InputArray Calorimeter/towers

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.2

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.2

  set JetPTMin 25.0
}

##################
# Fast Jet finder
##################

module FastJetFinder CaloJetFinder04 {
  set InputArray Calorimeter/towers

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.4

  set JetPTMin 25.0
}


##################
# Fast Jet finder
##################

module FastJetFinder CaloJetFinder08 {
  set InputArray Calorimeter/towers

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.8

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

  set JetPTMin 25.0
}

##################
# Fast Jet finder
##################

module FastJetFinder CaloJetFinder15 {
  set InputArray Calorimeter/towers

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.5

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 1.5

  set JetPTMin 25.0
}


##################
# Fast Jet finder
##################

module FastJetFinder TrackJetFinder02 {
  set InputArray TrackMerger/tracks

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.2

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.2

  set JetPTMin 25.0
}

##################
# Fast Jet finder
##################

module FastJetFinder TrackJetFinder04 {
  set InputArray TrackMerger/tracks

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.4

  set JetPTMin 25.0
}


##################
# Fast Jet finder
##################

module FastJetFinder TrackJetFinder08 {
  set InputArray TrackMerger/tracks

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.8

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

  set JetPTMin 25.0
}

##################
# Fast Jet finder
##################

module FastJetFinder TrackJetFinder15 {
  set InputArray TrackMerger/tracks

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.5

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 1.5

  set JetPTMin 25.0
}


##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder04/jets
  set OutputArray jets

 # scale formula for jets
  set ScaleFormula {1.00}
}


########################
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5
  set PartonPTMin 5.0
  set PartonEtaMax 6.0

}

###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray ECal/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  set EfficiencyFormula {
  (pt <= 1.0) * (0.00) + \
  (abs(eta) <= 2.5) * (pt > 1.0 && pt < 5.0)  * (0.70) +
  (abs(eta) <= 2.5) * (pt > 5.0 && pt < 10.0) * (0.85) +
  (abs(eta) <= 2.5) * (pt > 10.0)             * (0.90) +
    
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt < 5.0)  * (0.60) + 
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 5.0 && pt < 10.0) * (0.80) + 
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0)		* (0.90) +   
  
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 1.0 && pt < 5.0)  * (0.50) + \
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 5.0 && pt < 10.0) * (0.70) + \
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 10.0)             * (0.80) + \
  (abs(eta) > 6.0) * (0.00)}

}

##################
# Photon isolation
##################

# TBC: check values for iso cuts

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray photons

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.1
}


####################
# Electron isolation
####################

# TBC: check values for iso cuts

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronFilter/electrons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray electrons

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.1
}


################
# Muon isolation
################

# TBC: check values for iso cuts

module Isolation MuonIsolation {
  set CandidateInputArray MuonMomentumSmearing/muons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray muons

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.2
}



###########
# b-tagging
###########

module BTagging BTagging {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 0

  add EfficiencyFormula {0} {

  (pt <= 10.0)                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 10.0 && pt < 500)      * (0.01) + \
  (abs(eta) < 2.5)                   * (pt > 500.0 && pt < 15000.0) * (0.01)*(1.0 - pt/15000.) + \
  (abs(eta) < 2.5)                   * (pt > 15000.0)               * (0.00) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 500)      * (0.01) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 15000.0) * (0.01)*(1.0 - pt/15000.) + \
  (abs(eta) < 2.5 && abs(eta) < 4.0) * (pt > 15000.0)               * (0.000) + \
  (abs(eta) > 4.0) * (0.00)}

  add EfficiencyFormula {4} {

  (pt <= 10.0)                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 10.0 && pt < 500)      * (0.15) + \
  (abs(eta) < 2.5)                   * (pt > 500.0 && pt < 15000.0) * (0.15)*(1.0 - pt/15000.) + \
  (abs(eta) < 2.5)                   * (pt > 15000.0)               * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 500)      * (0.10) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 15000.0) * (0.10)*(1.0 - pt/15000.) + \
  (abs(eta) < 2.5 && abs(eta) < 4.0) * (pt > 15000.0)               * (0.000) + \
  (abs(eta) > 4.0) * (0.00)}

  add EfficiencyFormula {5} {

  (pt <= 10.0)                                                       * (0.00) +
  (abs(eta) < 2.5)                    * (pt > 10.0 && pt < 500)      * (0.82) + 
  (abs(eta) < 2.5)                    * (pt > 500.0 && pt < 15000.0) * (0.82)*(1.0 - pt/15000.) + 
  (abs(eta) < 2.5)                    * (pt > 15000.0)               * (0.000) + 
  (abs(eta) >= 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 500)      * (0.64) + 
  (abs(eta) >= 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 15000.0) * (0.64)*(1.0 - pt/15000.) + 
  (abs(eta) <= 2.5 && abs(eta) < 4.0) * (pt > 15000.0)               * (0.000) + 
  (abs(eta) >= 4.0) * (0.00)}

}

###########
# c-tagging
###########

module BTagging CTagging {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 1

  add EfficiencyFormula {0} {

  (pt <= 10.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 10.0) * (0.01) + \
  (abs(eta) > 4.0) * (pt > 10.0) * (0.00)}

  add EfficiencyFormula {4} {

  (pt <= 10.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 10.0) * (0.25) + \
  (abs(eta) > 4.0) * (pt > 10.0) * (0.00)}

  add EfficiencyFormula {5} {

  (pt <= 10.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 10.0) * (0.03) + \
  (abs(eta) > 4.0) * (pt > 10.0) * (0.00)}

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
  add EfficiencyFormula {0} {

  (pt <= 10.0)                                                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 10.0 && pt < 5000.0)    * (0.02) + \
  (abs(eta) < 2.5)                   * (pt > 5000.0 && pt < 34000.0) * (0.02)  *(8./9. - pt/30000.) + \
  (abs(eta) < 2.5)                   * (pt > 34000.0)                * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 5000.0)    * (0.0075) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 5000.0 && pt < 34000.0) * (0.0075)*(8./9. - pt/30000.) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 34000.0)                * (0.00) + \
  (abs(eta) > 4.0)                   * (0.00)}

  add EfficiencyFormula {11} {

  (pt <= 10.0)                                                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 10.0 && pt < 5000.0)    * (0.001) + \
  (abs(eta) < 2.5)                   * (pt > 5000.0 && pt < 34000.0) * (0.001)  *(8./9. - pt/30000.) + \
  (abs(eta) < 2.5)                   * (pt > 34000.0)                * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 5000.0)    * (0.001) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 5000.0 && pt < 34000.0) * (0.001)*(8./9. - pt/30000.) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 34000.0)                * (0.00) + \
  (abs(eta) > 4.0)                   * (0.00)}

  add EfficiencyFormula {15} {

  (pt <= 10.0)                                                       * (0.00) +
  (abs(eta) < 2.5) * (pt > 10.0 && pt < 5000.0)                      * (0.8)                      + \
  (abs(eta) < 2.5) * (pt > 5000.0 && pt < 34000.0)                   * (0.8) *(8./9. - pt/30000.) + \
  (abs(eta) < 2.5)                   * (pt > 34000.0)                * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 10.0 && pt < 5000.0)    * (0.65) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 5000.0 && pt < 34000.0) * (0.65)*(8./9. - pt/30000.) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 34000.0)                * (0.00) + \
  (abs(eta) > 4.0)                                                   * (0.00)}

}

#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

module UniqueObjectFinder UniqueObjectFinder {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray MuonIsolation/muons muons
  add InputArray JetEnergyScale/jets jets

}


##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle

  add Branch GenMissingET/momentum GenMissingET MissingET
  add Branch GenJetFinder02/jets GenJet Jet

  add Branch TrackMerger/tracks Track Track
  add Branch Calorimeter/towers Tower Tower

  add Branch HCal/eflowTracks EFlowTrack Track
  add Branch ECal/eflowPhotons EFlowPhoton Tower
  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch UniqueObjectFinder/photons Photon Photon
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/muons Muon Muon
  add Branch UniqueObjectFinder/jets Jet Jet

  add Branch GenJetFinder02/jets GenJet02 Jet
  add Branch GenJetFinder04/jets GenJet04 Jet
  add Branch GenJetFinder08/jets GenJet08 Jet
  add Branch GenJetFinder15/jets GenJet15 Jet

  add Branch FastJetFinder02/jets ParticleFlowJet02 Jet
  add Branch FastJetFinder04/jets ParticleFlowJet04 Jet
  add Branch FastJetFinder08/jets ParticleFlowJet08 Jet
  add Branch FastJetFinder15/jets ParticleFlowJet15 Jet

  add Branch CaloJetFinder02/jets CaloJet02 Jet
  add Branch CaloJetFinder04/jets CaloJet04 Jet
  add Branch CaloJetFinder08/jets CaloJet08 Jet
  add Branch CaloJetFinder15/jets CaloJet15 Jet

  add Branch TrackJetFinder02/jets TrackJet02 Jet
  add Branch TrackJetFinder04/jets TrackJet04 Jet
  add Branch TrackJetFinder08/jets TrackJet08 Jet
  add Branch TrackJetFinder15/jets TrackJet15 Jet

  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}

