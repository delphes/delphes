source delphes_card_CMS.tcl

set ExecutionPath [lreplace $ExecutionPath end end]
add ExecutionPath CaloJetFinder
add ExecutionPath PFJetFinder
add ExecutionPath CaloMissingET
add ExecutionPath PFMissingET
add ExecutionPath GenScalarHT
add ExecutionPath PionFilter
add ExecutionPath TreeWriter

module FastJetFinder GenJetFinder {
  set JetPTMin 1.0
}

module FastJetFinder FastJetFinder {
  set JetPTMin 1.0
}

########################
# Calorimeter jet finder
########################

module FastJetFinder CaloJetFinder {
  set InputArray Calorimeter/towers

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 1.0
}

########################
# PFJet finder
########################

module FastJetFinder PFJetFinder {
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 1.0
}

#########################
# Calo Missing ET merger
########################

module Merger CaloMissingET {
# add InputArray InputArray
  add InputArray Calorimeter/towers
  set MomentumOutputArray momentum
}

#########################
# Calo Missing ET merger
########################

module Merger PFMissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}


#################
# Gen Scalar HT
#################

module Merger GenScalarHT {
# add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set EnergyOutputArray energy
}

#################
# Pion filter
#################

module PdgCodeFilter PionFilter {
  set InputArray HCal/eflowTracks
  set OutputArray pions
  set Invert true
  add PdgCode {211}
  add PdgCode {-211}
}

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch CaloJetFinder/jets CaloJet Jet
  add Branch PFJetFinder/jets PFJet Jet
  add Branch CaloMissingET/momentum CaloMissingET MissingET
  add Branch PFMissingET/momentum PFMissingET MissingET
  add Branch GenScalarHT/energy GenScalarHT ScalarHT
  add Branch PionFilter/pions Pion Track 
  add Branch ElectronFilter/electrons ElectronPF Electron 
}
