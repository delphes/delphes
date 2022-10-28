set MaxEvents 1000
source delphes_card_IDEA.tcl

set ExecutionPath [lreplace $ExecutionPath end end]
add ExecutionPath CaloJetFinder
add ExecutionPath GenJetFlavorAssociation
add ExecutionPath GenJetTauTagging
add ExecutionPath CaloMissingET
add ExecutionPath PFMissingET
add ExecutionPath GenScalarHT
add ExecutionPath PionFilter
add ExecutionPath TreeWriter

module FastJetFinder GenJetFinder {
  set JetPTMin 10.0
}

module FastJetFinder FastJetFinder {
  set JetPTMin 10.0
}

##########################
# Jet Flavor Association
##########################

module JetFlavorAssociation GenJetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray GenJetFinder/jets

  set DeltaR 0.5
  set PartonPTMin 5.0
  set PartonEtaMax 6.0

}

#############
# tau-tagging
#############

module TauTagging GenJetTauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray GenJetFinder/jets

  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 4.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
  add EfficiencyFormula {0} {0.}
  add EfficiencyFormula {15} {1.0}

}


########################
# Calorimeter jet finder
########################

module FastJetFinder CaloJetFinder {
  set InputArray TowerMerger/towers

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 10
  set ParameterR 1.5
  set ParameterP -1.0
  set JetPTMin 10.0
}



#########################
# Calo Missing ET merger
########################

module Merger CaloMissingET {
# add InputArray InputArray
  add InputArray TowerMerger/towers
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
  set InputArray Calorimeter/eflowTracks
  set OutputArray pions
  set Invert true
  add PdgCode {211}
  add PdgCode {-211}
}

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch CaloJetFinder/jets CaloJet Jet
  add Branch FastJetFinder/jets PFJet Jet
  add Branch TrackMerger/tracks Track Track
  add Branch TowerMerger/towers Tower Tower
  add Branch CaloMissingET/momentum CaloMissingET MissingET
  add Branch PFMissingET/momentum PFMissingET MissingET
  add Branch GenScalarHT/energy GenScalarHT ScalarHT
  add Branch PionFilter/pions Pion Track
  add Branch ElectronFilter/electrons ElectronPF Electron
}
