source delphes_card_CMS.tcl

set ExecutionPath [lreplace $ExecutionPath end end]
add ExecutionPath CaloJetFinder
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
  set ParameterR 0.5

  set JetPTMin 1.0
}

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch CaloJetFinder/jets CaloJet Jet
}
