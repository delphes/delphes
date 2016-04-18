#
#  Phase II - No Pile-Up
#
#  Main authors: Michele Selvaggi (UCL)
#                                
#  Released on: 
#
#  Version: v01
#
#
#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
 
  PileUpMerger 
  
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  MuonMomentumSmearing

  TrackMerger
  
  ECal
  HCal

  ElectronFilter
  TrackPileUpSubtractor
  
  TowerMerger
  NeutralEFlowMerger
  EFlowMergerAllTracks
  EFlowMerger

  RunPUPPI

  PhotonEfficiency
  PhotonIsolation

  ElectronEfficiency
  ElectronIsolation

  MuonEfficiency
  MuonIsolation

  NeutrinoFilter

  MissingET
  GenMissingET
  GenPileUpMissingET

  GenJetFinder
  FastJetFinder

  ScalarHT
 
  JetEnergyScale

  JetFlavorAssociation

  BTaggingLoose
  BTaggingMedium
  BTaggingTight

  TauTagging
  
  GenParticleFilter
  
  TreeWriter
}



###############
# PileUp Merger
###############

module PileUpMerger PileUpMerger {
  set InputArray Delphes/stableParticles

  set ParticleOutputArray stableParticles
  set VertexOutputArray vertices

  # pre-generated minbias input file
  set PileUpFile MinBias.pileup

  # average expected pile up
  set MeanPileUp 200
  
  # maximum spread in the beam direction in m
  set ZVertexSpread 0.15

  # maximum spread in time in s
  set TVertexSpread 1.5E-09

  # vertex smearing formula f(z,t) (z,t need to be respectively given in m,s)
 
  #set VertexDistributionFormula {exp(-(t^2/(2*(0.05/2.99792458E8*exp(-(z^2/(2*(0.05)^2))))^2)))}
  
  
  set VertexDistributionFormula { (abs(t) <= 1.0e-09) * (abs(z) <= 0.15) * (1.00) +
                                  (abs(t) >  1.0e-09) * (abs(z) <= 0.15) * (0.00) +
   				  (abs(t) <= 1.0e-09) * (abs(z) > 0.15)  * (0.00) +
   				  (abs(t) >  1.0e-09) * (abs(z) > 0.15)  * (0.00)}

}


#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray PileUpMerger/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  # radius of the magnetic field coverage, in m
  set Radius 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.0

  # magnetic field
  set Bz 3.8
}


####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  ## particles after propagation
  set InputArray  ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  # tracking efficiency formula for charged hadrons
  set EfficiencyFormula {
      (pt <= 0.2) * (0.00) + \
	  (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.96) + \
	  (abs(eta) <= 1.2) * (pt > 1.0) * (0.97) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.85) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0) * (0.87) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.82) + \
	  (abs(eta) > 4.0) * (0.00)
  }
}


#####################################
# Electron tracking efficiency - ID
####################################

module Efficiency ElectronTrackingEfficiency {
  set InputArray  ParticlePropagator/electrons
  set OutputArray electrons
  # tracking efficiency formula for electrons
  set EfficiencyFormula { 
      (pt <= 0.2) * (0.00) + \
	  (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.96) + \
	  (abs(eta) <= 1.2) * (pt > 1.0) * (0.97) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.85) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0 && pt <= 10.0) * (0.82+pt*0.01) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 10.0) * (0.90) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.8) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0 && pt <= 10.0) * (0.8+pt*0.01) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 10.0) * (0.85) + \
	  (abs(eta) > 4.0) * (0.00)

  }
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons
  # tracking efficiency formula for muons
  set EfficiencyFormula {
      (pt <= 0.2) * (0.00) + \
	  (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 0.998) + \
	  (abs(eta) <= 1.2) * (pt > 1.0) * (0.998) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 0.2 && pt <= 1.0) * (pt*0.99) + \
	  (abs(eta) > 1.2 && abs(eta) <= 2.5) * (pt > 1.0) * (0.99) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.95) + \
	  (abs(eta) > 2.5 && abs(eta) <= 4.0) * (pt > 1.0) * (0.95) + \
	  (abs(eta) > 4.0) * (0.00)
	  
  }
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  ## hadrons after having applied the tracking efficiency
  set InputArray  ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons
  # resolution formula for charged hadrons
  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.013) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (0.017) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0   && pt <= 10.0)  * (0.03) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0  && pt <= 100.0) * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 1.0   && pt <= 10.0)  * (0.04) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 10.0  && pt <= 100.0) * (0.07) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.025) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 1.0   && pt <= 10.0)  * (0.05) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 10.0  && pt <= 100.0) * (0.20) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 100.0)                * (0.80)
  }
}

########################################
# Momentum resolution for electrons
########################################

module MomentumSmearing ElectronMomentumSmearing {
  ## hadrons after having applied the tracking efficiency
  set InputArray  ElectronTrackingEfficiency/electrons
  set OutputArray electrons
  # resolution formula for electrons
  set ResolutionFormula {                  (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.013) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.02) + \
                                           (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.015) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (0.017) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0   && pt <= 10.0)  * (0.03) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0  && pt <= 100.0) * (0.05) + \
                         (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 1.0   && pt <= 10.0)  * (0.04) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 10.0  && pt <= 100.0) * (0.07) + \
                         (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 100.0)                * (0.30) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.025) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 1.0   && pt <= 10.0)  * (0.05) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 10.0  && pt <= 100.0) * (0.20) + \
                         (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 100.0)                * (0.80)
  }
}


###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons
  # resolution formula for muons
  set ResolutionFormula {
    (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.012) + \
      (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e1) * (0.01) + \
      (abs(eta) <= 1.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.012) + \
      (abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.030) + \
      (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.011) + \
      (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e1) * (0.011) + \
      (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e1 && pt <= 2.0e2) * (0.022) + \
      (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 2.0e2)                * (0.030) + \
	  (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 0.1   && pt <= 1.0)   * (0.017) + \
	  (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 1.0   && pt <= 10.0)  * (0.03) + \
	  (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0  && pt <= 100.0) * (0.05) + \
	  (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 100.0)                * (0.30) + \
	  (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
	  (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 1.0   && pt <= 10.0)  * (0.04) + \
	  (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 10.0  && pt <= 100.0) * (0.07) + \
	  (abs(eta) > 3.0 && abs(eta) <= 3.5) * (pt > 100.0)                * (0.30) + \
	  (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 1.0)   * (0.025) + \
	  (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 1.0   && pt <= 10.0)  * (0.05) + \
	  (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 10.0  && pt <= 100.0) * (0.20) + \
	  (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 100.0)                * (0.80)
      
  }
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
 
  set EnergyMin 0.5
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # assume 0.02 x 0.02 resolution in eta,phi in the barrel |eta| < 1.5
  
  set PhiBins {}
  for {set i -157} {$i <= 157} {incr i} {
    add PhiBins [expr {$i * $pi/157.0}]
  }

  # 0.02 unit in eta up to eta = 1.5 (barrel)
  for {set i -74} {$i <= 75} {incr i} {
    set eta [expr {$i * 0.02}]
    add EtaPhiBins $eta $PhiBins
  }

  # assume 0.05 x 0.05 resolution in eta,phi in the endcaps 1.5 < |eta| < 3.0 (HGCAL)
  
  set PhiBins {}
  for {set i -65} {$i <= 65} {incr i} {
    add PhiBins [expr {$i * $pi/65.0}]
  }

  # 0.05 unit in eta up to eta = 3
  for {set i 1} {$i <= 30} {incr i} {
    set eta [expr { -3.00 + $i * 0.05}]
    add EtaPhiBins $eta $PhiBins
  }

  for {set i 1} {$i <= 30} {incr i} {
    set eta [expr { 1.5 + $i * 0.05}]
    add EtaPhiBins $eta $PhiBins
  }

  # take present CMS granularity for HF  
 
  # 0.175 x (0.175 - 0.35) resolution in eta,phi in the HF 3.0 < |eta| < 5.0
  set PhiBins {}
  for {set i -18} {$i <= 18} {incr i} {
    add PhiBins [expr {$i * $pi/18.0}]
  }
  
  foreach eta {-5 -4.7 -4.525 -4.35 -4.175 -4 -3.825 -3.65 -3.475 -3.3 -3.125 -3 3.125 3.3 3.475 3.65 3.825 4 4.175 4.35 4.525 4.7 5} {
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

  # set ResolutionFormula {resolution formula as a function of eta and energy}
  
  # for the ECAL barrel (|eta| < 1.5), see hep-ex/1306.2016.
  # for the endcaps (1.5 < |eta| < 3.0), we take HGCAL  see LHCC-P-008, Fig. 3.39, p.117

  set ResolutionFormula {                     (abs(eta) <= 1.50) * sqrt(energy^2*0.003^2 + energy*0.028^2 + 0.12^2) + \
                           (abs(eta) > 1.50 && abs(eta) <= 1.75) * sqrt(energy^2*0.006^2 + energy*0.20^2) + \
                           (abs(eta) > 1.75 && abs(eta) <= 2.15) * sqrt(energy^2*0.007^2 + energy*0.21^2) + \
                           (abs(eta) > 2.15 && abs(eta) <= 3.00) * sqrt(energy^2*0.008^2 + energy*0.24^2) + \
                           (abs(eta) > 3.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.08^2 + energy*1.97^2)}

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
 
  # assume 0.087 x 0.087 resolution in eta,phi in the barrel |eta| < 1.5

  set PhiBins {}
  for {set i -36} {$i <= 36} {incr i} {
    add PhiBins [expr {$i * $pi/36.0}]
  }
  foreach eta {-1.566 -1.479 -1.392 -1.305 -1.218 -1.131 -1.044 -0.957 -0.87 -0.783 -0.696 -0.609 -0.522 -0.435 -0.348 -0.261 -0.174 -0.087 0 0.087 0.174 0.261 0.348 0.435 0.522 0.609 0.696 0.783 0.87 0.957 1.044 1.131 1.218 1.305 1.392 1.479 1.566 1.65} {
    add EtaPhiBins $eta $PhiBins
  }

  # assume 0.05 x 0.05 resolution in eta,phi in the endcaps 1.5 < |eta| < 3.0 (HGCAL)

  set PhiBins {}
  for {set i -65} {$i <= 65} {incr i} {
    add PhiBins [expr {$i * $pi/65.0}]
  }

  # 0.05 unit in eta up to eta = 3
  for {set i 1} {$i <= 27} {incr i} {
    set eta [expr { -3.00 + $i * 0.05}]
    add EtaPhiBins $eta $PhiBins
  }

  for {set i 4} {$i <= 30} {incr i} {
    set eta [expr { 1.5 + $i * 0.05}]
    add EtaPhiBins $eta $PhiBins
  }  

  # take present CMS granularity for HF  
 
  # 0.175 x (0.175 - 0.35) resolution in eta,phi in the HF 3.0 < |eta| < 5.0
  set PhiBins {}
  for {set i -18} {$i <= 18} {incr i} {
    add PhiBins [expr {$i * $pi/18.0}]
  }
  
  foreach eta {-5 -4.7 -4.525 -4.35 -4.175 -4 -3.825 -3.65 -3.475 -3.3 -3.125 -3 3.125 3.3 3.475 3.65 3.825 4 4.175 4.35 4.525 4.7 5} {
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

# set ResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {                  (abs(eta) <= 1.7) * sqrt(energy^2*0.03^2 + energy*0.52^2 + 1.59^2) + \
						   (abs(eta) > 1.7 && abs(eta) <= 3.2) * sqrt(energy^2*0.05^2 + energy*0.706^2) + \
						   (abs(eta) > 3.0 && abs(eta) <= 4.9) * sqrt(energy^2*0.05^2 + energy*1.00^2)}

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

##########################
# Track pile-up subtractor
##########################

module TrackPileUpSubtractor TrackPileUpSubtractor {
# add InputArray InputArray OutputArray
  add InputArray HCal/eflowTracks eflowTracks
  add InputArray ElectronFilter/electrons electrons
  add InputArray MuonMomentumSmearing/muons muons

  set VertexInputArray PileUpMerger/vertices
  # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution
  # Z vertex resolution in m
  set ZVertexResolution 0.0001
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
# Neutral eflow erger
####################

module Merger NeutralEFlowMerger {
# add InputArray InputArray
  add InputArray ECal/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflowTowers
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

##################################
# Energy flow merger (all tracks)
##################################

module Merger EFlowMergerAllTracks {
# add InputArray InputArray
  add InputArray TrackMerger/tracks
  add InputArray ECal/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflow
}

#########################################
### Run the puppi code (to be tuned) ###
#########################################

module RunPUPPI RunPUPPI {
  ## input information
  set TrackInputArray   TrackMerger/tracks
  set NeutralInputArray NeutralEFlowMerger/eflowTowers
  set PVInputArray      PileUpMerger/vertices
  set MinPuppiWeight    0.05
  set UseExp            false
 
  ## define puppi algorithm parameters (more than one for the same eta region is possible)                                                                                      
  add EtaMinBin           0.    2.5    2.5    3.0   3.0
  add EtaMaxBin           2.5   3.0    3.0    10.0  10.0
  add PtMinBin            0.    0.5    0.5    0.5   0.5
  add ConeSizeBin         0.25  0.25   0.25   0.25  0.25
  add RMSPtMinBin         0.1   0.5    0.5    0.5   0.5
  add RMSScaleFactorBin   1.0   1.0    1.0    1.0   1.0
  add NeutralMinEBin      0.2   1.0    1.0    1.5   1.5
  add NeutralPtSlope      0.02  0.02   0.02   0.02  0.02
  add ApplyCHS            true  true   true   true  true
  add UseCharged          true  false  false  false false
  add ApplyLowPUCorr      true  true   true   true  true
  add MetricId            5     5      1      5     1

  ## output name
  set OutputArray         PuppiParticles
  set OutputArrayTracks   puppiTracks
  set OutputArrayNeutrals puppiNeutrals
} 



###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
#  add InputArray RunPUPPI/PuppiParticles
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

###################
# Ger PileUp Missing ET 
###################

module Merger GenPileUpMissingET {
# add InputArray InputArray
#  add InputArray RunPUPPI/PuppiParticles
  add InputArray ParticlePropagator/stableParticles
  set MomentumOutputArray momentum
}

##################
# Scalar HT merger
##################

module Merger ScalarHT {
# add InputArray InputArray
  add InputArray RunPUPPI/PuppiParticles
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
  set ParameterR 0.4

  set JetPTMin 15.0
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

module FastJetFinder FastJetFinder {
#  set InputArray TowerMerger/towers
  set InputArray RunPUPPI/PuppiParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 15.0
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


#####################
# Photon efficiency #
#####################

module Efficiency PhotonEfficiency {
  
  ## input particles
  set InputArray ECal/eflowPhotons
  ## output particles
  set OutputArray photons
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for photons
  set EfficiencyFormula {                      (pt <= 10.0) * (0.00) + \
	                   (abs(eta) <= 1.5) * (pt > 10.0)  * (0.9635) + \
	 (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 10.0)  * (0.9624) + \
	 (abs(eta) > 4.0)                                   * (0.00)}

}

####################
# Photon isolation #
####################

module Isolation PhotonIsolation {
  
  # particle for which calculate the isolation
  set CandidateInputArray PhotonEfficiency/photons 
  
  # isolation collection
  set IsolationInputArray RunPUPPI/PuppiParticles
 
  # output array
  set OutputArray photons
  
  # isolation cone
  set DeltaRMax 0.3
  
  # minimum pT  
  set PTMin     1.0
  
  # iso ratio to cut
  set PTRatioMax 9999.
}


#######################
# Electron efficiency #
#######################

module Efficiency ElectronEfficiency {
  
  set InputArray ElectronFilter/electrons
  set OutputArray electrons
  
  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for electrons
  set EfficiencyFormula {
                                      (pt <= 4.0)  * (0.00) + \
                         (abs(eta) <= 1.45 ) * (pt >  4.0 && pt <= 6.0)   * (0.50) + \
                         (abs(eta) <= 1.45 ) * (pt >  6.0 && pt <= 8.0)   * (0.70) + \
                         (abs(eta) <= 1.45 ) * (pt >  8.0 && pt <= 10.0)  * (0.85) + \
                         (abs(eta) <= 1.45 ) * (pt > 10.0 && pt <= 30.0)  * (0.94) + \                                                      
                         (abs(eta) <= 1.45 ) * (pt > 30.0 && pt <= 50.0)  * (0.97) + \                          
                         (abs(eta) <= 1.45 ) * (pt > 50.0 && pt <= 70.0)  * (0.98) + \          
                         (abs(eta) <= 1.45 ) * (pt > 70.0 )  * (1.0) + \                                                                                                 
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt >  4.0 && pt <= 10.0)   * (0.35) + \
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 10.0 && pt <= 30.0)   * (0.40) + \   
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 30.0 && pt <= 70.0)   * (0.45) + \                                 
                         (abs(eta) > 1.45  && abs(eta) <= 1.55) * (pt > 70.0 )  * (0.55) + \    
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt >  4.0 && pt <= 10.0)  * (0.75) + \
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 10.0 && pt <= 30.0)  * (0.85) + \                                                      
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 30.0 && pt <= 50.0)  * (0.95) + \                          
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 50.0 && pt <= 70.0)  * (0.95) + \          
                         (abs(eta) >= 1.55 && abs(eta) <= 2.0 ) * (pt > 70.0 )  * (1.0) + \   
                         (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt >  4.0 && pt <= 10.0)  * (0.65) + \
                         (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 10.0 && pt <= 30.0)  * (0.75) + \                                                      
                         (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 30.0 && pt <= 50.0)  * (0.90) + \                          
                         (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 50.0 && pt <= 70.0)  * (0.90) + \          
                         (abs(eta) >= 2.0 && abs(eta) <= 2.5 ) * (pt > 70.0 )  * (0.90) + \                                                                            
                         (abs(eta) > 2.5 && abs(eta) <= 4.0 ) * (pt > 4.0 && pt <= 10.0) * (0.65) + \
					  (abs(eta) > 2.5 && abs(eta) <= 4.0 ) * (pt > 10.0 && pt <= 30.0) * (0.75) + \
					  (abs(eta) > 2.5 && abs(eta) <= 4.0 ) * (pt > 30.0 && pt <= 50.0) * (0.90) + \
					  (abs(eta) > 2.5 && abs(eta) <= 4.0 ) * (pt > 50.0 && pt <= 70.0) * (0.90) + \
					  (abs(eta) > 2.5 && abs(eta) <= 4.0 ) * (pt > 70.0 ) * (0.90) + \
					  (abs(eta) > 4.0) * (0.00)
 
  }
}

######################
# Electron isolation #
######################

module Isolation ElectronIsolation {
  
  set CandidateInputArray ElectronEfficiency/electrons
  
  # isolation collection
  set IsolationInputArray RunPUPPI/PuppiParticles
  
  set OutputArray electrons
  
  set DeltaRMax 0.3
  set PTMin 1.0
  set PTRatioMax 9999.

}



###################
# Muon efficiency #
###################

module Efficiency MuonEfficiency {
  
  set InputArray MuonMomentumSmearing/muons
 
  set OutputArray muons
  # set EfficiencyFormula {efficiency as a function of eta and pt}
  # efficiency formula for muons
  set EfficiencyFormula { (pt <= 2.0)  * (0.00) + \
  			  (abs(eta) <= 4.00) * (pt >  2.0 && pt <= 3.0)  * (0.51) + \
                          (abs(eta) <= 4.00) * (pt >  3.0 && pt <= 4.0)  * (0.85) + \
                          (abs(eta) <= 4.00) * (pt >  4.0 && pt <= 11.0) * (0.93) + \
                          (abs(eta) <= 4.00) * (pt >  11. && pt <= 50.)  * (0.96) + \
                          (abs(eta) <= 4.00) * (pt >  50. && pt <= 70.)  * (0.98) + \
                          (abs(eta) <= 4.00) * (pt > 70.0 )  * (1.00) + \
			  (abs(eta) > 4.00)  * (0.00)
 }

}



##################
# Muon isolation #
##################

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
 
  # isolation collection
  set IsolationInputArray RunPUPPI/PuppiParticles
 
  set OutputArray muons
  
  set DeltaRMax 0.3
  set PTMin 1.0
  set PTRatioMax 9999.

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
  set PartonPTMin 10.0
  set PartonEtaMax 4.0

}


#############
# b-tagging #
#############
module BTagging BTaggingLoose {

  set JetInputArray JetEnergyScale/jets

  set BitNumber 0

  add EfficiencyFormula {0}      {(pt <= 20.0) * (0.000) + \
		                  (pt > 20.0) * (abs(eta) <= 3.4) * (0.1) + \
                                  (abs(eta) > 3.4) * (0.000) 
  }


  add EfficiencyFormula {4}          { (pt <= 20.0)                                    * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30)     * (0.341) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40)     * (0.389) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50)     * (0.416) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60)     * (0.432) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70)     * (0.447) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80)     * (0.451) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90)     * (0.451) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100)    * (0.447) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120)   * (0.447) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140)   * (0.444) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160)   * (0.430) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180)   * (0.424) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200)   * (0.419) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250)   * (0.407) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300)   * (0.384) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350)   * (0.364) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400)   * (0.353) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500)   * (0.303) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600)   * (0.273) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700)   * (0.249) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800)   * (0.233) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000)  * (0.209) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.194) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.227) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0)               * (0.000) + \

				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30)    * (0.236) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40)    * (0.276) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50)    * (0.302) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60)    * (0.330) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70)    * (0.338) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80)    * (0.355) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90)    * (0.353) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100)   * (0.358) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120)  * (0.350) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140)  * (0.327) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160)  * (0.313) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180)  * (0.307) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200)  * (0.322) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250)  * (0.299) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300)  * (0.274) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350)  * (0.240) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400)  * (0.203) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500)  * (0.218) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600)  * (0.203) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 3000) * (0.178) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 3000.0)              * (0.000) + \
				      
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 30)    * (0.179) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 30.0 && pt <= 40)    * (0.198) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 40.0 && pt <= 50)    * (0.219) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60)    * (0.225) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70)    * (0.232) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80)    * (0.251) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90)    * (0.267) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100)   * (0.252) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 120)  * (0.243) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 120.0 && pt <= 140)  * (0.244) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 140.0 && pt <= 160)  * (0.242) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 160.0 && pt <= 180)  * (0.246) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 180.0 && pt <= 200)  * (0.217) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 250)  * (0.223)+ \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 250.0 && pt <= 350)  * (0.194) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 350.0 && pt <= 3000) * (0.187) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 3000.0)              * (0.000) + \
				  
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 20.0 && pt <= 30)    * (0.09) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 30.0 && pt <= 40)    * (0.10) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 40.0 && pt <= 50)    * (0.11) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 50.0 && pt <= 60)    * (0.11) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 60.0 && pt <= 70)    * (0.12) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 70.0 && pt <= 80)    * (0.13) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 80.0 && pt <= 90)    * (0.13) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 90.0 && pt <= 100)   * (0.13) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 100.0 && pt <= 120)  * (0.12) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 120.0 && pt <= 140)  * (0.12) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 140.0 && pt <= 160)  * (0.12) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 160.0 && pt <= 180)  * (0.12) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 180.0 && pt <= 200)  * (0.11) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 200.0 && pt <= 250)  * (0.11) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 250.0 && pt <= 350)  * (0.10) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 350.0 && pt <= 3000) * (0.09) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 3000.0)              * (0.000) + \
				       (abs(eta) > 3.4)                                                 * (0.000)
  }

  add EfficiencyFormula {5}          { (pt <= 20.0)                                    * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30)     * (0.712) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40)     * (0.788) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50)     * (0.837) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60)     * (0.856) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70)     * (0.873) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80)     * (0.886) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90)     * (0.883) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100)    * (0.882) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120)   * (0.883) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140)   * (0.883) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160)   * (0.876) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180)   * (0.872) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200)   * (0.867) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250)   * (0.853) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300)   * (0.831) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350)   * (0.819) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400)   * (0.791) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500)   * (0.740) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600)   * (0.699) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700)   * (0.631) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800)   * (0.600) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000)  * (0.567) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.478) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.470) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0)               * (0.000) + \
				     
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30)    * (0.543) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40)    * (0.627) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50)    * (0.701) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60)    * (0.741) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70)    * (0.754) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80)    * (0.796) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90)    * (0.778) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100)   * (0.798) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120)  * (0.788) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140)  * (0.775) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160)  * (0.765) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180)  * (0.752) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200)  * (0.740) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250)  * (0.716) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300)  * (0.671) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350)  * (0.637) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400)  * (0.577) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500)  * (0.572) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600)  * (0.514) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 3000) * (0.516) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 3000.0)              * (0.000) + \
				      
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 30)    * (0.407) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 30.0 && pt <= 40)    * (0.479) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 40.0 && pt <= 50)    * (0.529) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60)    * (0.578) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70)    * (0.612) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80)    * (0.673) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90)    * (0.656) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100)   * (0.658) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 120)  * (0.647) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 120.0 && pt <= 140)  * (0.629) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 140.0 && pt <= 160)  * (0.630) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 160.0 && pt <= 180)  * (0.615) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 180.0 && pt <= 200)  * (0.604) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 250)  * (0.532)+ \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 250.0 && pt <= 350)  * (0.532) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 350.0 && pt <= 3000) * (0.503) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 3000.0)              * (0.000) + \
				  
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 20.0 && pt <= 30)    * (0.20) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 30.0 && pt <= 40)    * (0.235) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 40.0 && pt <= 50)    * (0.26) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 50.0 && pt <= 60)    * (0.28) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 60.0 && pt <= 70)    * (0.31) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 70.0 && pt <= 80)    * (0.33) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 80.0 && pt <= 90)    * (0.32) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 90.0 && pt <= 100)   * (0.32) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 100.0 && pt <= 120)  * (0.33) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 120.0 && pt <= 140)  * (0.31) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 140.0 && pt <= 160)  * (0.32) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 160.0 && pt <= 180)  * (0.31) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 180.0 && pt <= 200)  * (0.30) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 200.0 && pt <= 250)  * (0.26) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 250.0 && pt <= 350)  * (0.27) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 350.0 && pt <= 3000) * (0.25) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 3000.0)              * (0.000) + \
				       (abs(eta) > 3.4)                                                 * (0.000)
  }
}


#############
# b-tagging #
#############
module BTagging BTaggingMedium {

  set JetInputArray JetEnergyScale/jets

  set BitNumber 1
 
  add EfficiencyFormula {0}      {(pt <= 20.0) * (0.000) + \
		                          (pt > 20.0) * (abs(eta) <= 3.4) * (0.01) + \
                                  (abs(eta) > 3.4) * (0.000) 
  }

  add EfficiencyFormula {4}         { (pt <= 20.0)                                    * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30)     * (0.144) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40)     * (0.173) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50)     * (0.188) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60)     * (0.200) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70)     * (0.209) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80)     * (0.210) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90)     * (0.206) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100)    * (0.205) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120)   * (0.198) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140)   * (0.187) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160)   * (0.175) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180)   * (0.167) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200)   * (0.157) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250)   * (0.138) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300)   * (0.118) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350)   * (0.094) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400)   * (0.083) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500)   * (0.062) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600)   * (0.049) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700)   * (0.044) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800)   * (0.036) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000)  * (0.035) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.026) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.041) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0)               * (0.000) + \
				          
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30)    * (0.070) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40)    * (0.089) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50)    * (0.104) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60)    * (0.121) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70)    * (0.121) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80)    * (0.127) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90)    * (0.127) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100)   * (0.117) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120)  * (0.117) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140)  * (0.096) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160)  * (0.099) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180)  * (0.086) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200)  * (0.075) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250)  * (0.066) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300)  * (0.056) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350)  * (0.048) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400)  * (0.031) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500)  * (0.042) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600)  * (0.030) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 3000) * (0.030) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 3000.0)              * (0.000) + \
				      
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 30)    * (0.041) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 30.0 && pt <= 40)    * (0.047) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 40.0 && pt <= 50)    * (0.059) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60)    * (0.056) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70)    * (0.068) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80)    * (0.073) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90)    * (0.080) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100)   * (0.063) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 120)  * (0.065) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 120.0 && pt <= 140)  * (0.057) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 140.0 && pt <= 160)  * (0.049) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 160.0 && pt <= 180)  * (0.067) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 180.0 && pt <= 200)  * (0.065) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 250)  * (0.044)+ \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 250.0 && pt <= 350)  * (0.031) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 350.0 && pt <= 3000) * (0.026) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 3000.0)              * (0.000) + \
				  
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 20.0 && pt <= 30)    * (0.02) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 30.0 && pt <= 40)    * (0.02) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 40.0 && pt <= 50)    * (0.03) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 50.0 && pt <= 60)    * (0.03) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 60.0 && pt <= 70)    * (0.03) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 70.0 && pt <= 80)    * (0.035) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 80.0 && pt <= 90)    * (0.04) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 90.0 && pt <= 100)   * (0.03) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 100.0 && pt <= 120)  * (0.03) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 120.0 && pt <= 140)  * (0.025) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 140.0 && pt <= 160)  * (0.025) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 160.0 && pt <= 180)  * (0.035) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 180.0 && pt <= 200)  * (0.032) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 200.0 && pt <= 250)  * (0.022) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 250.0 && pt <= 350)  * (0.015) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 350.0 && pt <= 3000) * (0.013) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 3000.0)              * (0.000) + \
				       (abs(eta) > 3.4)                                                 * (0.000)
  }


  add EfficiencyFormula {5}        { (pt <= 20.0)                                    * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30)     * (0.528) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40)     * (0.611) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50)     * (0.666) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60)     * (0.701) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70)     * (0.734) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80)     * (0.741) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90)     * (0.746) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100)    * (0.750) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120)   * (0.742) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140)   * (0.729) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160)   * (0.711) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180)   * (0.694) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200)   * (0.684) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250)   * (0.642) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300)   * (0.602) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350)   * (0.554) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400)   * (0.488) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500)   * (0.417) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600)   * (0.343) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700)   * (0.284) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800)   * (0.250) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000)  * (0.221) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.164) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.173) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0)               * (0.000) + \

				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30)    * (0.322) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40)    * (0.414) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50)    * (0.488) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60)    * (0.521) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70)    * (0.554) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80)    * (0.590) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90)    * (0.571) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100)   * (0.583) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120)  * (0.575) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140)  * (0.557) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160)  * (0.536) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180)  * (0.516) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200)  * (0.483) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250)  * (0.432) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300)  * (0.366) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350)  * (0.323) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400)  * (0.241) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500)  * (0.252) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600)  * (0.237) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 3000) * (0.243) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 3000.0)              * (0.000) + \
				      
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 30)    * (0.200) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 30.0 && pt <= 40)    * (0.252) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 40.0 && pt <= 50)    * (0.310) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60)    * (0.353) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70)    * (0.377) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80)    * (0.425) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90)    * (0.406) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100)   * (0.414) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 120)  * (0.412) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 120.0 && pt <= 140)  * (0.366) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 140.0 && pt <= 160)  * (0.360) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 160.0 && pt <= 180)  * (0.347) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 180.0 && pt <= 200)  * (0.310) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 250)  * (0.307)+ \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 250.0 && pt <= 350)  * (0.224) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 350.0 && pt <= 3000) * (0.194) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 3000.0)              * (0.000) + \
				 
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 20.0 && pt <= 30)    * (0.100) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 30.0 && pt <= 40)    * (0.12) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 40.0 && pt <= 50)    * (0.15) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 50.0 && pt <= 60)    * (0.17) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 60.0 && pt <= 70)    * (0.19) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 70.0 && pt <= 80)    * (0.21) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 80.0 && pt <= 90)    * (0.20) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 90.0 && pt <= 100)   * (0.21) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 100.0 && pt <= 120)  * (0.20) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 120.0 && pt <= 140)  * (0.18) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 140.0 && pt <= 160)  * (0.18) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 160.0 && pt <= 180)  * (0.17) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 180.0 && pt <= 200)  * (0.15) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 200.0 && pt <= 250)  * (0.15) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 250.0 && pt <= 350)  * (0.11) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 350.0 && pt <= 3000) * (0.10) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 3000.0)              * (0.000) + \
				       (abs(eta) > 3.4)                                                 * (0.000)
  }
}

  

#############
# b-tagging #
#############
module BTagging BTaggingTight {

  set JetInputArray JetEnergyScale/jets

  set BitNumber 2
 
    add EfficiencyFormula {0}      {(pt <= 20.0) * (0.000) + \
		                    (pt > 20.0) * (abs(eta) <= 3.4) * (0.001) + \
                                    (abs(eta) > 3.4) * (0.000) 
  }


   add EfficiencyFormula {4}         { (pt <= 20.0)                                    * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30)     * (0.038) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40)     * (0.047) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50)     * (0.046) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60)     * (0.044) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70)     * (0.041) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80)     * (0.044) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90)     * (0.041) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100)    * (0.038) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120)   * (0.034) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140)   * (0.027) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160)   * (0.022) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180)   * (0.021) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200)   * (0.016) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250)   * (0.012) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300)   * (0.009) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350)   * (0.005) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400)   * (0.004) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500)   * (0.003) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600)   * (0.002) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700)   * (0.002) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800)   * (0.001) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000)  * (0.002) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.001) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.006) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0)               * (0.000) + \
				          
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30)    * (0.020) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40)    * (0.0242) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50)    * (0.0272) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60)    * (0.0349) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70)    * (0.0314) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80)    * (0.0220) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90)    * (0.028 ) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100)   * (0.0263) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120)  * (0.0245) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140)  * (0.0150) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160)  * (0.0088) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180)  * (0.0115) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200)  * (0.0063) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250)  * (0.0083) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300)  * (0.0027) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350)  * (0.0052) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400)  * (0.0007) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500)  * (0.0048) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600)  * (0.0013) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 3000) * (0.0013) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 3000.0)              * (0.000) + \
				      
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 30)    * (0.0092) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 30.0 && pt <= 40)    * (0.0094) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 40.0 && pt <= 50)    * (0.0141) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60)    * (0.0108) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70)    * (0.0161) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80)    * (0.0137) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90)    * (0.0135) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100)   * (0.0176) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 120)  * (0.0149) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 120.0 && pt <= 140)  * (0.0137) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 140.0 && pt <= 160)  * (0.0114) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 160.0 && pt <= 180)  * (0.0079) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 180.0 && pt <= 200)  * (0.0112) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 250)  * (0.0081)+ \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 250.0 && pt <= 350)  * (0.0033) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 350.0 && pt <= 3000) * (0.0016) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 3000.0)              * (0.000) + \
				 
				   (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 20.0 && pt <= 30)    * (0.005) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 30.0 && pt <= 40)    * (0.005) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 40.0 && pt <= 50)    * (0.005) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 50.0 && pt <= 60)    * (0.007) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 60.0 && pt <= 70)    * (0.005) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 70.0 && pt <= 80)    * (0.008) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 80.0 && pt <= 90)    * (0.006) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 90.0 && pt <= 100)   * (0.007) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 100.0 && pt <= 120)  * (0.008) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 120.0 && pt <= 140)  * (0.007) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 140.0 && pt <= 160)  * (0.006) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 160.0 && pt <= 180)  * (0.004) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 180.0 && pt <= 200)  * (0.005) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 200.0 && pt <= 250)  * (0.004) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 250.0 && pt <= 350)  * (0.001) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 350.0 && pt <= 3000) * (0.0005) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 3000.0)              * (0.000) + \
				       (abs(eta) > 3.4)                                                 * (0.000)
  }



  add EfficiencyFormula {5}        { (pt <= 20.0)                                    * (0.000) + \
				       (abs(eta) <= 1.8) * (pt > 20.0 && pt <= 30)     * (0.349) + \
				       (abs(eta) <= 1.8) * (pt > 30.0 && pt <= 40)     * (0.435) + \
				       (abs(eta) <= 1.8) * (pt > 40.0 && pt <= 50)     * (0.485) + \
				       (abs(eta) <= 1.8) * (pt > 50.0 && pt <= 60)     * (0.504) + \
				       (abs(eta) <= 1.8) * (pt > 60.0 && pt <= 70)     * (0.519) + \
				       (abs(eta) <= 1.8) * (pt > 70.0 && pt <= 80)     * (0.537) + \
				       (abs(eta) <= 1.8) * (pt > 80.0 && pt <= 90)     * (0.541) + \
				       (abs(eta) <= 1.8) * (pt > 90.0 && pt <= 100)    * (0.528) + \
				       (abs(eta) <= 1.8) * (pt > 100.0 && pt <= 120)   * (0.514) + \
				       (abs(eta) <= 1.8) * (pt > 120.0 && pt <= 140)   * (0.473) + \
				       (abs(eta) <= 1.8) * (pt > 140.0 && pt <= 160)   * (0.447) + \
				       (abs(eta) <= 1.8) * (pt > 160.0 && pt <= 180)   * (0.408) + \
				       (abs(eta) <= 1.8) * (pt > 180.0 && pt <= 200)   * (0.373) + \
				       (abs(eta) <= 1.8) * (pt > 200.0 && pt <= 250)   * (0.315) + \
				       (abs(eta) <= 1.8) * (pt > 250.0 && pt <= 300)   * (0.267) + \
				       (abs(eta) <= 1.8) * (pt > 300.0 && pt <= 350)   * (0.181) + \
				       (abs(eta) <= 1.8) * (pt > 350.0 && pt <= 400)   * (0.138) + \
				       (abs(eta) <= 1.8) * (pt > 400.0 && pt <= 500)   * (0.101) + \
				       (abs(eta) <= 1.8) * (pt > 500.0 && pt <= 600)   * (0.067) + \
				       (abs(eta) <= 1.8) * (pt > 600.0 && pt <= 700)   * (0.048) + \
				       (abs(eta) <= 1.8) * (pt > 700.0 && pt <= 800)   * (0.032) + \
				       (abs(eta) <= 1.8) * (pt > 800.0 && pt <= 1000)  * (0.029) + \
				       (abs(eta) <= 1.8) * (pt > 1000.0 && pt <= 1400) * (0.022) + \
				       (abs(eta) <= 1.8) * (pt > 1400.0 && pt <= 2000) * (0.038) + \
				       (abs(eta) <= 1.8) * (pt > 2000.0)               * (0.000) + \
				          
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 20.0 && pt <= 30)    * (0.186) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 30.0 && pt <= 40)    * (0.258) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 40.0 && pt <= 50)    * (0.306) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 50.0 && pt <= 60)    * (0.348) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 60.0 && pt <= 70)    * (0.359) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 70.0 && pt <= 80)    * (0.358) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 80.0 && pt <= 90)    * (0.367) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 90.0 && pt <= 100)   * (0.380) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 100.0 && pt <= 120)  * (0.359) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 120.0 && pt <= 140)  * (0.308) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 140.0 && pt <= 160)  * (0.238) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 160.0 && pt <= 180)  * (0.252) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 180.0 && pt <= 200)  * (0.229) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 200.0 && pt <= 250)  * (0.174) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 250.0 && pt <= 300)  * (0.125) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 300.0 && pt <= 350)  * (0.096) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 350.0 && pt <= 400)  * (0.066) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 400.0 && pt <= 500)  * (0.046) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 500.0 && pt <= 600)  * (0.076) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 600.0 && pt <= 3000) * (0.072) + \
				       (abs(eta) > 1.8 && abs(eta) <= 2.4) * (pt > 3000.0)              * (0.000) + \
				      
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 30)    * (0.096) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 30.0 && pt <= 40)    * (0.119) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 40.0 && pt <= 50)    * (0.156) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 60)    * (0.174) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 60.0 && pt <= 70)    * (0.192) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 70.0 && pt <= 80)    * (0.228) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 80.0 && pt <= 90)    * (0.218) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 90.0 && pt <= 100)   * (0.231) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 100.0 && pt <= 120)  * (0.224) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 120.0 && pt <= 140)  * (0.184) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 140.0 && pt <= 160)  * (0.164) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 160.0 && pt <= 180)  * (0.115) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 180.0 && pt <= 200)  * (0.123) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 200.0 && pt <= 250)  * (0.131)+ \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 250.0 && pt <= 350)  * (0.063) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 350.0 && pt <= 3000) * (0.039) + \
				       (abs(eta) > 2.4 && abs(eta) <= 3.0) * (pt > 3000.0)              * (0.000) + \
				  
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt <= 20.0)               * (0.000) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 20.0 && pt <= 30)    * (0.05) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 30.0 && pt <= 40)    * (0.06) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 40.0 && pt <= 50)    * (0.08) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 50.0 && pt <= 60)    * (0.09) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 60.0 && pt <= 70)    * (0.10) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 70.0 && pt <= 80)    * (0.11) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 80.0 && pt <= 90)    * (0.11) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 90.0 && pt <= 100)   * (0.11) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 100.0 && pt <= 120)  * (0.11) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 120.0 && pt <= 140)  * (0.09) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 140.0 && pt <= 160)  * (0.08) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 160.0 && pt <= 180)  * (0.06) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 180.0 && pt <= 200)  * (0.06) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 200.0 && pt <= 250)  * (0.06) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 250.0 && pt <= 350)  * (0.03) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 350.0 && pt <= 3000) * (0.02) + \
				       (abs(eta) > 3.0 && abs(eta) <= 3.4) * (pt > 3000.0)              * (0.000) + \
				       (abs(eta) > 3.4)                                                 * (0.000)
  }
}



#############
# tau-tagging
#############


module TauTagging TauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5

  set TauPTMin 20.0

  set TauEtaMax 2.3

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  add EfficiencyFormula {0}  { (abs(eta) < 2.3) * ((( -0.00621816+0.00130097*pt-2.19642e-5*pt^2+1.49393e-7*pt^3-4.58972e-10*pt^4+5.27983e-13*pt^5 )) * (pt<250) + 0.0032*(pt>250)) + \
                               (abs(eta) > 2.3) * (0.000)
                             }
  add EfficiencyFormula {15} { (abs(eta) < 2.3) * 0.97*0.77*( (0.32 + 0.01*pt - 0.000054*pt*pt )*(pt<100)+0.78*(pt>100) ) + \
                               (abs(eta) > 2.3) * (0.000)
                             } 
}


###############################################################################################################
# StatusPidFilter: this module removes all generated particles except electrons, muons, taus, and status == 3 #
###############################################################################################################

module StatusPidFilter GenParticleFilter {
  
    set InputArray  Delphes/allParticles
    set OutputArray filteredParticles
    set PTMin 5.0
    
}


##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch GenParticleFilter/filteredParticles Particle GenParticle
  add Branch PileUpMerger/vertices Vertex Vertex
 
  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenMissingET/momentum GenMissingET MissingET

#  add Branch HCal/eflowTracks EFlowTrack Track
#  add Branch ECal/eflowPhotons EFlowPhoton Tower
#  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower
  
  add Branch PhotonIsolation/photons Photon Photon
  add Branch ElectronIsolation/electrons Electron Electron
  add Branch MuonIsolation/muons Muon Muon
  
  add Branch JetEnergyScale/jets Jet Jet
  
  add Branch MissingET/momentum MissingET MissingET
  add Branch GenPileUpMissingET/momentum GenPileUpMissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}

