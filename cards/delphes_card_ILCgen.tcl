################################################
#
# ILCgen model
#
# Generic ILC detector model for Delphes,
# developed for the Snowmass 2021 study.
# While it has been mainly based on the ILD
# detector concept, as presented in ILD IDR,
# it can be considered a generic ILC detector
# model, as expected performances of both ILD
# and SiD are very similar and details of the
# detector design are not taken into account.
#
# For more details and references see:
#   https://github.com/iLCSoft/ILCDelphes
#
# Questions and comments can be send to:
# Aleksander Filip Zarnecki <zarnecki@fuw.edu.pl>
#
################################################

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  NeutrinoFilter

  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  MuonMomentumSmearing

  TrackMerger


    
  ECal
  LumiCalF
  LumiCalR

  HCal
  LHCalF
  LHCalR
    
  Calorimeter

  
  
  BeamCalF
  BeamCalR

  BCalTowers
  BCalMerger
  BCalEfficiency
   
  HCalMerger  
  ElectronFilter

    
  EFlowMerger
  PhotonMerger
  NeutralMerger

  PhotonEfficiency
  PhotonIsolation
    
  ElectronEfficiency
  ElectronIsolation

  MuonEfficiency
  MuonIsolation

  
  EFlowFilter

  GenJetFinder

  JetFinder
  JetFlavorAssociation
  BTagging80
  BTagging70
  BTagging50
  CTagging55
  CTagging30
  CTagging20
  TauTagging

  
  JetFinder_N2
  JetFlavor_N2
  BTagging80_N2
  BTagging70_N2
  BTagging50_N2
  CTagging55_N2
  CTagging30_N2
  CTagging20_N2
  TauTagging_N2
    

  JetFinder_N3
  JetFlavor_N3
  BTagging80_N3
  BTagging70_N3
  BTagging50_N3
  CTagging55_N3
  CTagging30_N3
  CTagging20_N3
  TauTagging_N3
    

  JetFinder_N4
  JetFlavor_N4
  BTagging80_N4
  BTagging70_N4
  BTagging50_N4
  CTagging55_N4
  CTagging30_N4
  CTagging20_N4
  TauTagging_N4
    

  JetFinder_N5
  JetFlavor_N5
  BTagging80_N5
  BTagging70_N5
  BTagging50_N5
  CTagging55_N5
  CTagging30_N5
  CTagging20_N5
  TauTagging_N5
    

  JetFinder_N6
  JetFlavor_N6
  BTagging80_N6
  BTagging70_N6
  BTagging50_N6
  CTagging55_N6
  CTagging30_N6
  CTagging20_N6
  TauTagging_N6
    

  MissingET
  GenMissingET
  ScalarHT


    
  MainCalorimeter

  EFlowMerger_MainCal

  PhotonEfficiency_MainCal
  PhotonIsolation_MainCal

  EFlowFilter_MainCal

  JetFinder_MainCal
  JetFlavorAssociation_MainCal
  BTagging80_MainCal
  BTagging70_MainCal
  BTagging50_MainCal
  CTagging55_MainCal
  CTagging30_MainCal
  CTagging20_MainCal
  TauTagging_MainCal

  MissingET_MainCal
  ScalarHT_MainCal

    
    
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
    
    source ILCgen/ILCgen_Propagator.tcl
    
}

####################################
# Charged hadron tracking efficiency
####################################
module Efficiency ChargedHadronTrackingEfficiency {
    set InputArray ParticlePropagator/chargedHadrons
    set OutputArray chargedHadrons


    source ILCgen/ILCgen_ChrgHadTrackingEff.tcl 
}

##############################
# Electron tracking efficiency
##############################
module Efficiency ElectronTrackingEfficiency {
    set InputArray ParticlePropagator/electrons
    set OutputArray electrons

    source ILCgen/ILCgen_ElectronTrackingEff.tcl 
}

##########################
# Muon tracking efficiency
##########################
module Efficiency MuonTrackingEfficiency {
    set InputArray ParticlePropagator/muons
    set OutputArray muons

    source ILCgen/ILCgen_MuonTrackingEff.tcl 
}

########################################
# Momentum resolution for charged tracks
########################################
module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

    source ILCgen/ILCgen_ChrgHadMomentumSmearing.tcl

}

###################################
# Momentum resolution for electrons
###################################
module MomentumSmearing ElectronMomentumSmearing {
    set InputArray ElectronTrackingEfficiency/electrons
    set OutputArray electrons

    source ILCgen/ILCgen_ElectronMomentumSmearing.tcl 
}

###############################
# Momentum resolution for muons
###############################
module MomentumSmearing MuonMomentumSmearing {
    set InputArray MuonTrackingEfficiency/muons
    set OutputArray muons

    source ILCgen/ILCgen_MuonMomentumSmearing.tcl 
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
module SimpleCalorimeter ECal {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray TrackMerger/tracks

    set TowerOutputArray ecalTowers
    
    set EFlowTrackOutputArray eflowTracks
    set EFlowTowerOutputArray eflowPhotons
    
    set IsEcal true
    
    set EnergyMin 0.1
    set EnergySignificanceMin 1.0
    
    set SmearTowerCenter true
    
    source ILCgen/ILCgen_ECAL_Binning.tcl
    source ILCgen/ILCgen_ECAL_EnergyFractions.tcl
    source ILCgen/ILCgen_ECAL_Resolution.tcl 
}

##############
# LumiCal
##############
module SimpleCalorimeter LumiCalF {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray TrackMerger/tracks
    
    set TowerOutputArray lumicalTowers
    
    set EFlowTrackOutputArray eflowTracks
    set EFlowTowerOutputArray eflowPhotons
    
    set IsEcal true 
    
    set EnergyMin 2.0
    set EnergySignificanceMin 1.0
    
    set SmearTowerCenter true
    
    source ILCgen/ILCgen_LumiCalF_Binning.tcl
    source ILCgen/ILCgen_ECAL_EnergyFractions.tcl
    source ILCgen/ILCgen_ECAL_Resolution.tcl
}

module SimpleCalorimeter LumiCalR {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray TrackMerger/tracks
    
    set TowerOutputArray lumicalTowers
    
    set EFlowTrackOutputArray eflowTracks
    set EFlowTowerOutputArray eflowPhotons
    
    set IsEcal true 
    
    set EnergyMin 2.0
    set EnergySignificanceMin 1.0
    
    set SmearTowerCenter true
    
    source ILCgen/ILCgen_LumiCalR_Binning.tcl
    source ILCgen/ILCgen_ECAL_EnergyFractions.tcl
    source ILCgen/ILCgen_ECAL_Resolution.tcl
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
    
    set EnergyMin 0.25
    set EnergySignificanceMin 1.0
    
    set SmearTowerCenter true
    
    source ILCgen/ILCgen_HCAL_Binning.tcl
    source ILCgen/ILCgen_HCAL_EnergyFractions.tcl
    source ILCgen/ILCgen_HCAL_Resolution.tcl 
}

##############
# LHCal
##############
module SimpleCalorimeter LHCalR {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray ECal/eflowTracks
    
    set TowerOutputArray lhcalTowers
    set EFlowTrackOutputArray eflowTracks
    set EFlowTowerOutputArray eflowNeutralHadrons
    
    set IsEcal false 
    
    set EnergyMin 2.0
    set EnergySignificanceMin 1.0
    
    set SmearTowerCenter true
    
    source ILCgen/ILCgen_LHCalR_Binning.tcl
    source ILCgen/ILCgen_HCAL_EnergyFractions.tcl
    source ILCgen/ILCgen_HCAL_Resolution.tcl
}

module SimpleCalorimeter LHCalF {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray ECal/eflowTracks
    
    set TowerOutputArray lhcalTowers
    set EFlowTrackOutputArray eflowTracks
    set EFlowTowerOutputArray eflowNeutralHadrons
    
    set IsEcal false 
    
    set EnergyMin 2.0
    set EnergySignificanceMin 1.0
    
    set SmearTowerCenter true
    
    source ILCgen/ILCgen_LHCalF_Binning.tcl
    source ILCgen/ILCgen_HCAL_EnergyFractions.tcl
    source ILCgen/ILCgen_HCAL_Resolution.tcl
}

##############
# BeamCal
##############
module SimpleCalorimeter BeamCalR {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray TrackMerger/tracks
    
    set TowerOutputArray bcalTowers
    set EFlowTowerOutputArray bcalPhotons
    
    set IsEcal true 
    
    set EnergyMin 5.0
    set EnergySignificanceMin 1.0
    
    set SmearTowerCenter true
    
    source ILCgen/ILCgen_BeamCalR_Binning.tcl
    source ILCgen/ILCgen_BeamCal_EnergyFractions.tcl
    source ILCgen/ILCgen_BeamCal_Resolution.tcl
}

module SimpleCalorimeter BeamCalF {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray TrackMerger/tracks
    
    set TowerOutputArray bcalTowers
    set EFlowTowerOutputArray bcalPhotons
    
    set IsEcal true 
    
    set EnergyMin 5.0

    set EnergySignificanceMin 1.0
    
    set SmearTowerCenter true
    
    source ILCgen/ILCgen_BeamCalF_Binning.tcl
    source ILCgen/ILCgen_BeamCal_EnergyFractions.tcl
    source ILCgen/ILCgen_BeamCal_Resolution.tcl
}

#################
# Electron merger
#################
module Merger HCalMerger {
# add InputArray InputArray
  add InputArray HCal/eflowTracks
  add InputArray LHCalF/eflowTracks
  add InputArray LHCalR/eflowTracks
  set OutputArray eflowTracks
}

#################
# Electron filter
#################
module PdgCodeFilter ElectronFilter {
  set InputArray  HCalMerger/eflowTracks
  set OutputArray electrons
  set Invert true

  add PdgCode {11}
  add PdgCode {-11}
}

###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################
module Merger Calorimeter {
# add InputArray InputArray
  add InputArray ECal/ecalTowers
  add InputArray HCal/hcalTowers
  add InputArray LumiCalF/lumicalTowers
  add InputArray LumiCalR/lumicalTowers
  add InputArray LHCalF/lhcalTowers
  add InputArray LHCalR/lhcalTowers
  set OutputArray towers
}

############################################
# Tower Merger for central calorimeters only
############################################
module Merger MainCalorimeter {
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
  add InputArray HCalMerger/eflowTracks
  add InputArray ECal/eflowPhotons
  add InputArray LumiCalF/eflowPhotons
  add InputArray LumiCalR/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  add InputArray LHCalF/eflowNeutralHadrons
  add InputArray LHCalR/eflowNeutralHadrons
  set OutputArray eflow
}

##################################################
# Energy flow merger for central calorimeters only
##################################################
module Merger EFlowMerger_MainCal {
# add InputArray InputArray
  add InputArray HCal/eflowTracks
  add InputArray ECal/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflow
}

###############
# Photon merger
###############
module Merger PhotonMerger {
# add InputArray InputArray
  add InputArray ECal/eflowPhotons
  add InputArray LumiCalF/eflowPhotons
  add InputArray LumiCalR/eflowPhotons

  set OutputArray eflowPhotons
}

#######################
# Neutral hadron merger
#######################
module Merger NeutralMerger {
# add InputArray InputArray
  add InputArray HCal/eflowNeutralHadrons
  add InputArray LHCalF/eflowNeutralHadrons
  add InputArray LHCalR/eflowNeutralHadrons
  set OutputArray eflowNeutralHadrons
}

###############################
# BeamCal tower merger
###############################
module Merger BCalTowers {
# add InputArray InputArray
  add InputArray BeamCalF/bcalTowers
  add InputArray BeamCalR/bcalTowers
  set OutputArray bcalTowers
}
  
###############################
# BeamCal energy flow merger
###############################
module Merger BCalMerger {
# add InputArray InputArray
  add InputArray BeamCalF/bcalPhotons
  add InputArray BeamCalR/bcalPhotons
  set OutputArray bcalPhotons
}
  
##############################
# BeamCal photon efficiency
##############################
module Efficiency BCalEfficiency {
    set InputArray  BCalMerger/bcalPhotons
    set OutputArray bcalPhotons

    source ILCgen/ILCgen_BeamCalEfficiency.tcl
}

###################
# Photon efficiency
###################
module Efficiency PhotonEfficiency {
    set InputArray PhotonMerger/eflowPhotons
    set OutputArray photons

    source ILCgen/ILCgen_PhotonEfficiency.tcl
}

##################
# Photon isolation
##################
module Isolation PhotonIsolation {
    set CandidateInputArray PhotonEfficiency/photons
    set IsolationInputArray EFlowMerger/eflow
    set OutputArray photons

    source ILCgen/ILCgen_PhotonIsolation.tcl
}

#####################
# Electron efficiency
#####################
module Efficiency ElectronEfficiency {
    set InputArray ElectronFilter/electrons
    set OutputArray electrons

    source ILCgen/ILCgen_ElectronEfficiency.tcl
}

####################
# Electron isolation
####################
module Isolation ElectronIsolation {
    set CandidateInputArray ElectronEfficiency/electrons
    set IsolationInputArray EFlowMerger/eflow
    set OutputArray electrons
    
    source ILCgen/ILCgen_ElectronIsolation.tcl
}

#################
# Muon efficiency
#################
module Efficiency MuonEfficiency {
    set InputArray MuonMomentumSmearing/muons
    set OutputArray muons

    source ILCgen/ILCgen_MuonEfficiency.tcl
}

################
# Muon isolation
################
module Isolation MuonIsolation {
    set CandidateInputArray MuonEfficiency/muons
    set IsolationInputArray EFlowMerger/eflow
    set OutputArray muons
    
    source ILCgen/ILCgen_MuonIsolation.tcl
}

##################################
# EFlowFilter (UniqueObjectFinder)
##################################
module UniqueObjectFinder EFlowFilter {
    add InputArray PhotonIsolation/photons photons
    add InputArray ElectronIsolation/electrons electrons
    add InputArray MuonIsolation/muons muons
    add InputArray EFlowMerger/eflow eflow
}

###################
# Missing ET merger
###################
module Merger MissingET {
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

module Merger MissingET_MainCal {
  add InputArray EFlowMerger_MainCal/eflow
  set MomentumOutputArray momentum
}

##################
# Scalar HT merger
##################
module Merger ScalarHT {
  add InputArray EFlowMerger/eflow
  set EnergyOutputArray energy
}

module Merger ScalarHT_MainCal {
  add InputArray EFlowMerger_MainCal/eflow
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

  source ILCgen/ILCgen_GenJetFinder.tcl
}

#########################
# Gen Missing ET merger
########################
module Merger GenMissingET {
# add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}

############
# Jet finder
############
module FastJetFinder JetFinder {
    #set InputArray TowerMerger/towers
    set InputArray EFlowFilter/eflow
    set OutputArray jets
    
    source ILCgen/ILCgen_JetFinder.tcl
}

########################
# Jet Flavor Association
########################
module JetFlavorAssociation JetFlavorAssociation {

    set PartonInputArray Delphes/partons
    set ParticleInputArray Delphes/allParticles
    set ParticleLHEFInputArray Delphes/allParticlesLHEF
    set JetInputArray JetFinder/jets

    source ILCgen/ILCgen_JetFlavourAssoc.tcl
}

###########
# b-tagging
###########
module BTagging BTagging80 {
    set JetInputArray JetFinder/jets
    set BitNumber 0

    source ILCgen/ILCgen_BTagging_80.tcl
}

module BTagging BTagging70 {
    set JetInputArray JetFinder/jets
    set BitNumber 1

    source ILCgen/ILCgen_BTagging_70.tcl
}

module BTagging BTagging50 {
    set JetInputArray JetFinder/jets
    set BitNumber 2

    source ILCgen/ILCgen_BTagging_50.tcl
}

###########
# c-tagging
###########
module BTagging CTagging55 {
    set JetInputArray JetFinder/jets
    set BitNumber 4

    source ILCgen/ILCgen_CTagging_55.tcl
}

module BTagging CTagging30 {
    set JetInputArray JetFinder/jets
    set BitNumber 5

    source ILCgen/ILCgen_CTagging_30.tcl
}

module BTagging CTagging20 {
    set JetInputArray JetFinder/jets
    set BitNumber 6

    source ILCgen/ILCgen_CTagging_20.tcl
}

#############
# tau-tagging
#############
module TauTagging TauTagging {
    set ParticleInputArray Delphes/allParticles
    set PartonInputArray Delphes/partons
    set JetInputArray JetFinder/jets

    source ILCgen/ILCgen_TauTagging.tcl
}

##############################################
# Jet finder for inclusive clustering, N=2...6
##############################################
source ILCgen/ILCgen_JetFinder_N.tcl

##########################################################
# Jet Flavor Association for inclusive clustering, N=2...6
##########################################################
source ILCgen/ILCgen_JetFlavourAssoc_N.tcl

#############################################
# b-tagging for inclusive clustering, N=2...6
#############################################
source ILCgen/ILCgen_BTagging_N.tcl

#############################################
# c-tagging for inclusive clustering, N=2...6
#############################################
source ILCgen/ILCgen_CTagging_N.tcl

###############################################
# tau-tagging for inclusive clustering, N=2...6
###############################################
source ILCgen/ILCgen_TauTagging_N.tcl

####################################
# Photon efficiency central detector
####################################
module Efficiency PhotonEfficiency_MainCal {
    set InputArray ECal/eflowPhotons
    set OutputArray photons

    source ILCgen/ILCgen_PhotonEfficiency.tcl
}

####################################
# Photon isolation central detectors
####################################
module Isolation PhotonIsolation_MainCal {
    set CandidateInputArray PhotonEfficiency_MainCal/photons
    set IsolationInputArray EFlowMerger_MainCal/eflow
    set OutputArray photons

    source ILCgen/ILCgen_PhotonIsolation.tcl
}

######################################
# EFlowFilter for central calorimeters
######################################
module UniqueObjectFinder EFlowFilter_MainCal {
    add InputArray PhotonIsolation_MainCal/photons photons
    add InputArray ElectronIsolation/electrons electrons
    add InputArray MuonIsolation/muons muons
    add InputArray EFlowMerger_MainCal/eflow eflow
}

#######################################
# Jet finder for central detectors only
#######################################
module FastJetFinder JetFinder_MainCal {
    #set InputArray TowerMerger/towers
    set InputArray EFlowFilter_MainCal/eflow
    set OutputArray jets
    
    source ILCgen/ILCgen_JetFinder.tcl
}

#########################################
# Jet Flavor Association for central jets
#########################################
module JetFlavorAssociation JetFlavorAssociation_MainCal {

    set PartonInputArray Delphes/partons
    set ParticleInputArray Delphes/allParticles
    set ParticleLHEFInputArray Delphes/allParticlesLHEF
    set JetInputArray JetFinder_MainCal/jets

    source ILCgen/ILCgen_JetFlavourAssoc.tcl
}

############################
# b-tagging for central jets
############################
module BTagging BTagging80_MainCal {
    set JetInputArray JetFinder_MainCal/jets
    set BitNumber 0

    source ILCgen/ILCgen_BTagging_80.tcl
}

module BTagging BTagging70_MainCal {
    set JetInputArray JetFinder_MainCal/jets
    set BitNumber 1

    source ILCgen/ILCgen_BTagging_70.tcl
}

module BTagging BTagging50_MainCal {
    set JetInputArray JetFinder_MainCal/jets
    set BitNumber 2

    source ILCgen/ILCgen_BTagging_50.tcl
}

############################
# c-tagging for central jets
############################
module BTagging CTagging55_MainCal {
    set JetInputArray JetFinder_MainCal/jets
    set BitNumber 4

    source ILCgen/ILCgen_CTagging_55.tcl
}

module BTagging CTagging30_MainCal {
    set JetInputArray JetFinder_MainCal/jets
    set BitNumber 5

    source ILCgen/ILCgen_CTagging_30.tcl
}

module BTagging CTagging20_MainCal {
    set JetInputArray JetFinder_MainCal/jets
    set BitNumber 6

    source ILCgen/ILCgen_CTagging_20.tcl
}

##########################
# tau-tagging central jets
##########################
module TauTagging TauTagging_MainCal {
    set ParticleInputArray Delphes/allParticles
    set PartonInputArray Delphes/partons
    set JetInputArray JetFinder_MainCal/jets

    source ILCgen/ILCgen_TauTagging.tcl
}

##################
# ROOT tree writer
##################
module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass

#
# Generator level 
#
  add Branch Delphes/allParticles Particle GenParticle
  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenMissingET/momentum GenMissingET MissingET

#
# Raw detector response
#     uncomment for tests or for Dephes event display 
#
#  add Branch TrackMerger/tracks Track Track
#  add Branch Calorimeter/towers Tower Tower

#
# Additional raw data collections
#     for tests only
#
#  add Branch MainCalorimeter/towers Tower_MainCal Tower
#  add Branch BCalTowers/bcalTowers BCalTower Tower

#
# Particle flow objects 
#
  add Branch HCalMerger/eflowTracks EFlowTrack Track
  add Branch PhotonMerger/eflowPhotons EFlowPhoton Tower
  add Branch NeutralMerger/eflowNeutralHadrons EFlowNeutralHadron Tower

#
# Particle flow objects for main calorimeters only,
# exclusing LumiCal and LHCal - for tests only 
#
#  add Branch ECal/eflowPhotons EFlowPhoton_MainCal Tower
#  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron_MainCal Tower
#
  
#
# Final state reconstruction
#  
  add Branch EFlowFilter/electrons Electron Electron
  add Branch EFlowFilter/muons Muon Muon
  add Branch EFlowFilter/photons Photon Photon
  
  add Branch JetFinder/jets Jet Jet
  
  add Branch JetFinder_N2/jets Jet_N2 Jet
  add Branch JetFinder_N3/jets Jet_N3 Jet
  add Branch JetFinder_N4/jets Jet_N4 Jet
  add Branch JetFinder_N5/jets Jet_N5 Jet
  add Branch JetFinder_N6/jets Jet_N6 Jet

#
# Final state reconstruction with no LumiCal/LHCal
#   for special cases only, not in the default output stream  
#  
#  add Branch EFlowFilter_MainCal/photons Photon_MainCal Photon
#  add Branch JetFinder_MainCal/jets Jet_MainCal Jet

#
# Missing transverse momentum and transverse energy
# (vector and scalar sum of particle flow object momenta)
#  
  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT

#  
# Missing transverse momentum and transverse energy
#  without LumiCal and LHCal objects - for special cases only
#
#  add Branch MissingET_MainCal/momentum MissingET_MainCal MissingET
#  add Branch ScalarHT_MainCal/energy ScalarHT_MainCal ScalarHT

# BeamCal photons - not included in particle flow/clustering
#  nor in the transverse momentym/energy calculation
  
  add Branch BCalEfficiency/bcalPhotons BCalPhoton Photon

}

