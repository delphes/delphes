set RandomSeed 123

#
# Official Delphes card prepared by FCC-hh collaboration
#
#  Main authors:  Michele Selvaggi (CERN)
#
#  Released on: July 27th, 2022
#
#
#  Configuration: FCC-hh Scenario I (optimistic)
#
#
#   - removed DenseTrackFilter (not relevant only for pT> 5 TeV jet studies)
#   - nominal CLD detector tracking efficiency
#       -  https://indico.cern.ch/event/1064327/contributions/4893180/attachments/2452299/4204992/220601_FCCWeek_CLD_sailer.pdf
#   - nominal FCChh detector tracking/muon resolution for scenario I
#   - use DualRedoutCalorimeter module for simplicity:
#      - ECAL is CMS crystals (test beam)
#      - HCAL is DualReadout
#   - BTagging (3WPs) shape based on latest ParticleNet Run II (still private):
#      - nominal PNet will be scenario II, for scenario I(III) improve(degrade) btagging eff by +/-5%
#   - TauTagging (3WPs) based on DeepJet CMS:
#      - https://arxiv.org/pdf/2201.08458.pdf
#
#   TODO:
#
#
#   - optimise PF parameters using mbb and maa and eventually on boosted jets
#   - check that PF parameters optimisation does not affect electron isolation at high pT
#   - check if can propagate correctly BTag WPs in FCCSW
#   - FatJet tagging?
#
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
  ForwardLooperTracks

  Calorimeter
  EFlowTrackMerger
  EFlowMerger
  EFlowFilter

  PhotonEfficiency
  PhotonIsolation

  ElectronFilter
  ElectronEfficiency
  ElectronIsolation

  MuonFilter
  MuonEfficiency
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

  TowerMerger
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

  BTaggingL
  BTaggingM
  BTaggingT

  CTaggingL
  CTaggingM
  CTaggingT

  TauTaggingL
  TauTaggingM
  TauTaggingT

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



####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  set UseMomentumVector true

 # TBC (which eta_max ? which pT min?)

 # tracking efficiency formula for charged hadrons

  set EfficiencyFormula {

  (pt <= 0.1) * (0.00) +

  (abs(eta) <= 2.5) * (pt > 0.1 && pt <= 0.2) * (0.93) +
  (abs(eta) <= 2.5) * (pt > 0.2 && pt <= 0.3) * (0.95) +
  (abs(eta) <= 2.5) * (pt > 0.3 && pt <= 0.5) * (0.97) +
  (abs(eta) <= 2.5) * (pt > 0.5 && pt <= 1) * (0.99) +
  (abs(eta) <= 2.5) * (pt > 1) * (1.00) +

  0.95 * (
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.1 && pt <= 0.2) * (0.93) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.2 && pt <= 0.3) * (0.95) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.3 && pt <= 0.5) * (0.97) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.5 && pt <= 1) * (0.99) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 1) * (1.00)
  )+
  0.90 * (
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.1 && pt <= 0.2) * (0.93) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.2 && pt <= 0.3) * (0.95) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.3 && pt <= 0.5) * (0.97) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.5 && pt <= 1) * (0.99) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 1) * (1.00)
  )+
  (abs(eta) > 6.0) * (0.00)
  }

}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons
  set UseMomentumVector true

  set EfficiencyFormula {

  (pt <= 0.1) * (0.00) +

  (abs(eta) <= 2.5) * (pt > 0.1 && pt <= 0.2) * (0.93) +
  (abs(eta) <= 2.5) * (pt > 0.2 && pt <= 0.3) * (0.95) +
  (abs(eta) <= 2.5) * (pt > 0.3 && pt <= 0.5) * (0.97) +
  (abs(eta) <= 2.5) * (pt > 0.5 && pt <= 1) * (0.99) +
  (abs(eta) <= 2.5) * (pt > 1) * (1.00) +

  0.95 * (
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.1 && pt <= 0.2) * (0.93) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.2 && pt <= 0.3) * (0.95) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.3 && pt <= 0.5) * (0.97) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.5 && pt <= 1) * (0.99) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 1) * (1.00)
  )+
  0.90 * (
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.1 && pt <= 0.2) * (0.93) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.2 && pt <= 0.3) * (0.95) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.3 && pt <= 0.5) * (0.97) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.5 && pt <= 1) * (0.99) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 1) * (1.00)
  )+
  (abs(eta) > 6.0) * (0.00)
  }


}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons
  set UseMomentumVector true

  set EfficiencyFormula {

  (pt <= 0.1) * (0.00) +

  (abs(eta) <= 2.5) * (pt > 0.1 && pt <= 0.2) * (0.99) +
  (abs(eta) <= 2.5) * (pt > 0.2 && pt <= 0.3) * (0.99) +
  (abs(eta) <= 2.5) * (pt > 0.3 && pt <= 0.5) * (0.99) +
  (abs(eta) <= 2.5) * (pt > 0.5 && pt <= 1) * (0.99) +
  (abs(eta) <= 2.5) * (pt > 1) * (1.00) +

  0.975 * (
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.1 && pt <= 0.2) * (0.99) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.2 && pt <= 0.3) * (0.99) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.3 && pt <= 0.5) * (0.99) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 0.5 && pt <= 1) * (0.99) +
    (abs(eta) > 2.5 && abs(eta) <= 4) * (pt > 1) * (1.00)
  )+
  0.95 * (
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.1 && pt <= 0.2) * (0.99) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.2 && pt <= 0.3) * (0.99) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.3 && pt <= 0.5) * (0.99) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 0.5 && pt <= 1) * (0.99) +
    (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 1) * (1.00)
  )+
  (abs(eta) > 6.0) * (0.00)
  }

}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  source trackMomentumResolution_I.tcl
}


###################################
# Momentum resolution for electrons
###################################

module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  source electronMomentumResolution_I.tcl
}


###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  source muonMomentumResolution_I.tcl
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

######################
# Looper Selection
######################

module Efficiency ForwardLooperTracks  {
  set InputArray TrackMerger/tracks
  set OutputArray tracks
  set UseMomentumVector False

  ## select looping tracks that end up in position |eta| > 6.000 (lost by calo)
  set EfficiencyFormula {
    (abs(eta) > 6.0 )                                 * (1.000) +
    (abs(eta) <= 6.0 )                                * (0.000)
  }

}

#############
# Calorimeter
#############
module DualReadoutCalorimeter Calorimeter {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray towers
  set PhotonOutputArray photons

  set EFlowTrackOutputArray eflowTracks
  set EFlowPhotonOutputArray eflowPhotons
  set EFlowNeutralHadronOutputArray eflowNeutralHadrons

  set ECalMinSignificance 2.0
  set HCalMinSignificance 2.5

  set SmearLogNormal false

  set SmearTowerCenter true
  #set SmearTowerCenter false
    set pi [expr {acos(-1)}]

    # Lists of the edges of each tower in eta and phi;
    # each list starts with the lower edge of the first tower;
    # the list ends with the higher edged of the last tower.
    # Barrel:  deta=0.02 towers up to |eta| <= 0.88 ( up to 45°)
    # Endcaps: deta=0.02 towers up to |eta| <= 3.0 (8.6° = 100 mrad)
    # Cell size: about 6 cm x 6 cm

    set EtaPhiRes 0.02
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

    # default energy fractions {abs(PDG code)} {Fecal Fhcal}
    add EnergyFraction {0} {0.0 1.0}
    # energy fractions for e, gamma and pi0
    add EnergyFraction {11} {1.0 0.0}
    add EnergyFraction {22} {1.0 0.0}
    add EnergyFraction {111} {1.0 0.0}
    # energy fractions for muon, neutrinos and neutralinos
    add EnergyFraction {12} {0.0 0.0}
    add EnergyFraction {13} {0.0 0.0}
    add EnergyFraction {14} {0.0 0.0}
    add EnergyFraction {16} {0.0 0.0}
    add EnergyFraction {1000022} {0.0 0.0}
    add EnergyFraction {1000023} {0.0 0.0}
    add EnergyFraction {1000025} {0.0 0.0}
    add EnergyFraction {1000035} {0.0 0.0}
    add EnergyFraction {1000045} {0.0 0.0}
    # energy fractions for K0short and Lambda
    add EnergyFraction {310} {0.3 0.7}
    add EnergyFraction {130} {0.3 0.7}
    add EnergyFraction {3122} {0.3 0.7}


    ## ECAL crystals for the EM part from 2008.00338
    # set ECalResolutionFormula {resolution formula as a function of eta and energy}
    set ECalResolutionFormula {
    (abs(eta) <= 0.88 )                     * sqrt(energy^2*0.005^2 + energy*0.03^2 + 0.002^2)+
    (abs(eta) > 0.88 && abs(eta) <= 6.0)    * sqrt(energy^2*0.005^2 + energy*0.03^2 + 0.002^2)
    }


    # Dual Readout
    # set HCalResolutionFormula {resolution formula as a function of eta and energy}
    set HCalResolutionFormula {
    (abs(eta) <= 0.88 )                     * sqrt(energy^2*0.01^2 + energy*0.3^2 + 0.05^2)+
    (abs(eta) > 0.88 && abs(eta) <= 6.0)    * sqrt(energy^2*0.01^2 + energy*0.3^2 + 0.05^2)
    }
}



#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray EFlowTrackMerger/eflowTracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}


#################
# Muon filter
#################

module PdgCodeFilter MuonFilter {
  set InputArray EFlowTrackMerger/eflowTracks
  set OutputArray muons
  set Invert true
  add PdgCode {13}
  add PdgCode {-13}
}



###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

module Merger TowerMerger {
# add InputArray InputArray
  add InputArray Calorimeter/towers
  add InputArray MuonFilter/muons
  set OutputArray towers
}



############################
# Energy flow track merger
############################

module Merger EFlowTrackMerger {
# add InputArray InputArray
  add InputArray Calorimeter/eflowTracks
  add InputArray ForwardLooperTracks/tracks
  set OutputArray eflowTracks
}


####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray EFlowTrackMerger/eflowTracks
  add InputArray Calorimeter/eflowPhotons
  add InputArray Calorimeter/eflowNeutralHadrons
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
  set InputArray TowerMerger/towers

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
  set InputArray TowerMerger/towers

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
  set InputArray TowerMerger/towers

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
  set InputArray TowerMerger/towers

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
  set InputArray Calorimeter/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  set EfficiencyFormula {
  (pt <= 1.0) * (0.00) + \
  (abs(eta) <= 2.5) * (pt > 1.0 && pt < 5.0)  * (0.80) +
  (abs(eta) <= 2.5) * (pt > 5.0 && pt < 10.0) * (0.90) +
  (abs(eta) <= 2.5) * (pt > 10.0)             * (0.95) +

  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt < 5.0)  * (0.80) +
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 5.0 && pt < 10.0) * (0.85) +
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0)		* (0.90) +

  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 1.0 && pt < 5.0)  * (0.75) + \
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 5.0 && pt < 10.0) * (0.80) + \
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 10.0)             * (0.85) + \
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

  # veto isolation cand. based on proximity to input cand. -> turned out not to be needed
  # set DeltaRMin 0.02
  # set UseMiniCone true

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.1
}


###################
# Electron efficiency
###################

module Efficiency ElectronEfficiency {
  set InputArray ElectronFilter/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  set EfficiencyFormula {
  (pt <= 1.0) * (0.00) + \
  (abs(eta) <= 2.5) * (pt > 1.0 && pt < 5.0)  * (0.80) +
  (abs(eta) <= 2.5) * (pt > 5.0 && pt < 10.0) * (0.90) +
  (abs(eta) <= 2.5) * (pt > 10.0)             * (0.95) +

  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt < 5.0)  * (0.80) +
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 5.0 && pt < 10.0) * (0.85) +
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0)		* (0.90) +

  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 1.0 && pt < 5.0)  * (0.75) + \
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 5.0 && pt < 10.0) * (0.80) + \
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 10.0)             * (0.85) + \
  (abs(eta) > 6.0) * (0.00)}

}


####################
# Electron isolation
####################

# TBC: check values for iso cuts

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray electrons

  # veto isolation cand. based on proximity to input cand. -> turned out not to be needed
  # set DeltaRMin 0.02
  # set UseMiniCone true

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.1
}




###################
# Muon efficiency
###################

module Efficiency MuonEfficiency {
  set InputArray MuonFilter/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  set EfficiencyFormula {
  (pt <= 1.0) * (0.00) + \
  (abs(eta) <= 2.5) * (pt > 1.0 && pt < 5.0)  * (0.95) +
  (abs(eta) <= 2.5) * (pt > 5.0 && pt < 10.0) * (0.97) +
  (abs(eta) <= 2.5) * (pt > 10.0)             * (0.999) +

  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt < 5.0)  * (0.90) +
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 5.0 && pt < 10.0) * (0.95) +
  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0)		* (0.99) +

  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 1.0 && pt < 5.0)  * (0.85) + \
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 5.0 && pt < 10.0) * (0.90) + \
  (abs(eta) > 4.0 && abs(eta) <= 6.0) * (pt > 10.0)             * (0.95) + \
  (abs(eta) > 6.0) * (0.00)}

}



################
# Muon isolation
################

# TBC: check values for iso cuts

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray muons

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.2
}



###############
# b-tagging L
###############

module BTagging BTaggingL {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 0

  add EfficiencyFormula {5} {

  1.05 * (
  (pt <= 20.0)                                                        * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0 && pt < 50)       * (0.88) +
  (abs(eta) < 4.0)                     * (pt > 50.0 && pt < 100)      * (0.92) +
  (abs(eta) < 4.0)                     * (pt > 100.0 && pt < 200)     * (0.93) +
  (abs(eta) < 4.0)                     * (pt > 200.0 && pt < 500)     * (0.90) +
  (abs(eta) < 4.0)                     * (pt > 500.0 && pt < 1000)    * (0.87) +
  (abs(eta) < 4.0)                     * (pt > 1000.0)                * (0.85) +
  (abs(eta) >= 4.0) * (0.00)
  )

  }

  add EfficiencyFormula {0} {

  (pt <= 20.0)                                                        * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0 && pt < 50)       * (0.15) +
  (abs(eta) < 4.0)                     * (pt > 50.0 && pt < 100)      * (0.08) +
  (abs(eta) < 4.0)                     * (pt > 100.0 && pt < 200)     * (0.06) +
  (abs(eta) < 4.0)                     * (pt > 200.0 && pt < 500)     * (0.08) +
  (abs(eta) < 4.0)                     * (pt > 500.0 && pt < 1000)    * (0.10) +
  (abs(eta) < 4.0)                     * (pt > 1000.0)                * (0.12) +
  (abs(eta) >= 4.0) * (0.00)

  }

  add EfficiencyFormula {4} {

  (pt <= 20.0)                                                        * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0 && pt < 50)       * (0.50) +
  (abs(eta) < 4.0)                     * (pt > 50.0 && pt < 100)      * (0.45) +
  (abs(eta) < 4.0)                     * (pt > 100.0 && pt < 200)     * (0.40) +
  (abs(eta) < 4.0)                     * (pt > 200.0 && pt < 500)     * (0.40) +
  (abs(eta) < 4.0)                     * (pt > 500.0 && pt < 1000)    * (0.40) +
  (abs(eta) < 4.0)                     * (pt > 1000.0)                * (0.40) +
  (abs(eta) >= 4.0) * (0.00)

  }
}

###############
# b-tagging M
###############

module BTagging BTaggingM {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 1

  add EfficiencyFormula {5} {

  1.05 * (
  (pt <= 20.0)                                                        * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0 && pt < 50)       * (0.75) +
  (abs(eta) < 4.0)                     * (pt > 50.0 && pt < 100)      * (0.84) +
  (abs(eta) < 4.0)                     * (pt > 100.0 && pt < 200)     * (0.85) +
  (abs(eta) < 4.0)                     * (pt > 200.0 && pt < 500)     * (0.83) +
  (abs(eta) < 4.0)                     * (pt > 500.0 && pt < 1000)    * (0.80) +
  (abs(eta) < 4.0)                     * (pt > 1000.0)                * (0.80) +
  (abs(eta) >= 4.0) * (0.00)
  )

  }

  add EfficiencyFormula {0} {

  (pt <= 20.0)                                                        * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0 && pt < 50)       * (0.015) +
  (abs(eta) < 4.0)                     * (pt > 50.0 && pt < 100)      * (0.008) +
  (abs(eta) < 4.0)                     * (pt > 100.0 && pt < 200)     * (0.006) +
  (abs(eta) < 4.0)                     * (pt > 200.0 && pt < 500)     * (0.008) +
  (abs(eta) < 4.0)                     * (pt > 500.0 && pt < 1000)    * (0.010) +
  (abs(eta) < 4.0)                     * (pt > 1000.0)                * (0.012) +
  (abs(eta) >= 4.0) * (0.00)

  }

  add EfficiencyFormula {4} {

  (pt <= 20.0)                                                        * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0 && pt < 50)       * (0.20) +
  (abs(eta) < 4.0)                     * (pt > 50.0 && pt < 100)      * (0.15) +
  (abs(eta) < 4.0)                     * (pt > 100.0 && pt < 200)     * (0.13) +
  (abs(eta) < 4.0)                     * (pt > 200.0 && pt < 500)     * (0.12) +
  (abs(eta) < 4.0)                     * (pt > 500.0 && pt < 1000)    * (0.11) +
  (abs(eta) < 4.0)                     * (pt > 1000.0)                * (0.10) +
  (abs(eta) >= 4.0) * (0.00)

  }
}


###############
# b-tagging T
###############

module BTagging BTaggingT {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 2

  add EfficiencyFormula {5} {

  1.05 * (
  (pt <= 20.0)                                                        * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0 && pt < 50)       * (0.55) +
  (abs(eta) < 4.0)                     * (pt > 50.0 && pt < 100)      * (0.65) +
  (abs(eta) < 4.0)                     * (pt > 100.0 && pt < 200)     * (0.70) +
  (abs(eta) < 4.0)                     * (pt > 200.0 && pt < 500)     * (0.65) +
  (abs(eta) < 4.0)                     * (pt > 500.0 && pt < 1000)    * (0.60) +
  (abs(eta) < 4.0)                     * (pt > 1000.0)                * (0.55) +
  (abs(eta) >= 4.0) * (0.00)
  )

  }

  add EfficiencyFormula {0} {

  (pt <= 20.0)                                                        * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0 && pt < 50)       * (0.005) +
  (abs(eta) < 4.0)                     * (pt > 50.0 && pt < 100)      * (0.005) +
  (abs(eta) < 4.0)                     * (pt > 100.0 && pt < 200)     * (0.006) +
  (abs(eta) < 4.0)                     * (pt > 200.0 && pt < 500)     * (0.007) +
  (abs(eta) < 4.0)                     * (pt > 500.0 && pt < 1000)    * (0.007) +
  (abs(eta) < 4.0)                     * (pt > 1000.0)                * (0.008) +
  (abs(eta) >= 4.0) * (0.00)

  }

  add EfficiencyFormula {4} {

  (pt <= 20.0)                                                        * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0 && pt < 50)       * (0.025) +
  (abs(eta) < 4.0)                     * (pt > 50.0 && pt < 100)      * (0.025) +
  (abs(eta) < 4.0)                     * (pt > 100.0 && pt < 200)     * (0.025) +
  (abs(eta) < 4.0)                     * (pt > 200.0 && pt < 500)     * (0.025) +
  (abs(eta) < 4.0)                     * (pt > 500.0 && pt < 1000)    * (0.02) +
  (abs(eta) < 4.0)                     * (pt > 1000.0)                * (0.02) +
  (abs(eta) >= 4.0) * (0.00)

  }
}


###############
# c-tagging L
###############

module BTagging CTaggingL {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 4

  add EfficiencyFormula {4} {

  1.05 * (
  (pt <= 20.0)                                             * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0)       * (0.95)
  )

  }

  add EfficiencyFormula {5} {

  (pt <= 20.0)                                             * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0)       * (0.25)

  }

  add EfficiencyFormula {0} {

  (pt <= 20.0)                                             * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0)       * (0.80)

  }
}


###############
# c-tagging M
###############

module BTagging CTaggingM {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 5

  add EfficiencyFormula {4} {

  1.05 * (
  (pt <= 20.0)                                             * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0)       * (0.65)
  )

  }

  add EfficiencyFormula {5} {

  (pt <= 20.0)                                             * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0)       * (0.07)

  }

  add EfficiencyFormula {0} {

  (pt <= 20.0)                                             * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0)       * (0.14)

  }
}



###############
# c-tagging T
###############

module BTagging CTaggingT {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 6

  add EfficiencyFormula {4} {

  1.05 * (
  (pt <= 20.0)                                             * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0)       * (0.40)
  )

  }

  add EfficiencyFormula {5} {

  (pt <= 20.0)                                             * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0)       * (0.025)

  }

  add EfficiencyFormula {0} {

  (pt <= 20.0)                                             * (0.00) +
  (abs(eta) < 4.0)                     * (pt > 20.0)       * (0.01)

  }
}



#################
# tau-tagging L
#################


module TauTagging TauTaggingL {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set BitNumber 0

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  add EfficiencyFormula {15} {
  1.05 *
  ((pt <= 20.0)                                                       * (0.00) +
  (abs(eta) < 4.0)                   * (pt > 20.0 && pt <= 100)      * (0.90) +
  (abs(eta) < 4.0)                   * (pt > 100)                    * (0.94) +
  (abs(eta) > 4.0)                   * (0.00)
  )
  }

  add EfficiencyFormula {0} {

  (pt <= 20.0)                                                       * (0.00) +
  (abs(eta) < 4.0)                   * (pt > 20.0 && pt <= 100)      * (0.035) +
  (abs(eta) < 4.0)                   * (pt > 100)                    * (0.02) +
  (abs(eta) > 4.0)                   * (0.00)
  }

}



#################
# tau-tagging M
#################


module TauTagging TauTaggingM {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set BitNumber 1

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  add EfficiencyFormula {15} {

  1.05 * (
  (pt <= 20.0)                                                       * (0.00) +
  (abs(eta) < 4.0)                   * (pt > 20.0 && pt <= 100)      * (0.80) +
  (abs(eta) < 4.0)                   * (pt > 100)                    * (0.88) +
  (abs(eta) > 4.0)                   * (0.00)
  )
  }

  add EfficiencyFormula {0} {

  (pt <= 20.0)                                                       * (0.00) +
  (abs(eta) < 4.0)                   * (pt > 20.0 && pt <= 100)      * (0.02) +
  (abs(eta) < 4.0)                   * (pt > 100)                    * (0.01) +
  (abs(eta) > 4.0)                   * (0.00)
  }

}


#################
# tau-tagging T
#################


module TauTagging TauTaggingT {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set BitNumber 2

  set DeltaR 0.5

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  add EfficiencyFormula {15} {

  1.05 * (
  (pt <= 20.0)                                                       * (0.00) +
  (abs(eta) < 4.0)                   * (pt > 20.0 && pt <= 100)      * (0.70) +
  (abs(eta) < 4.0)                   * (pt > 100)                    * (0.82) +
  (abs(eta) > 4.0)                   * (0.00)
  )

  }

  add EfficiencyFormula {0} {

  (pt <= 20.0)                                                       * (0.00) +
  (abs(eta) < 4.0)                   * (pt > 20.0 && pt <= 100)      * (0.01) +
  (abs(eta) < 4.0)                   * (pt > 100)                    * (0.008) +
  (abs(eta) > 4.0)                   * (0.00)
  }

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

  #add Branch TrackMerger/tracks Track Track
  #add Branch TowerMerger/towers Tower Tower

  #Temporary addition for manual isoVar validation/recalculation
  # add Branch EFlowFilter/eflow ParticleFlowCandidates ParticleFlowCandidate 

  add Branch EFlowTrackMerger/eflowTracks EFlowTrack Track
  add Branch Calorimeter/eflowPhotons EFlowPhoton Tower
  add Branch Calorimeter/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch UniqueObjectFinder/photons Photon Photon
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/muons Muon Muon
  add Branch UniqueObjectFinder/jets Jet Jet 

  #collections for objects before isolation and uniqueobjectfinder
  add Branch PhotonEfficiency/photons PhotonNoIso Photon
  add Branch ElectronEfficiency/electrons ElectronNoIso Electron
  add Branch MuonEfficiency/muons MuonNoIso Muon
  add Branch JetEnergyScale/jets JetNoIso Jet

  #collections for objects after isolation but before uniqueobjectfinder -> not needed?
  # add Branch PhotonIsolation/photons PhotonNoOR Photon
  # add Branch ElectronIsolation/electrons ElectronNoOR Electron
  # add Branch MuonIsolation/muons MuonNoOR Muon

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
