#######################################
# Uta Klein
# Order of execution of various modules
# 19.12.16 add paramterisations according to Peter Kostka's design
# muons not considered here - take CMS/FCC values
# 22.12.16 new extensions for FCC-eh based on new LHeC cards
# needs eflow and photons and JetFlavour Ass
# 3.3.17 add more Filter 16.3.17 smaller constant CAL terms c.f. Max Klein 
# 17.3.17 add track resolutions as given by Peter Kostka
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
  ImpactParameterSmearing
  Calorimeter
  EFlowMerger
  EFlowFilter

  PhotonEfficiency
  PhotonIsolation
 
  ElectronFilter
  ElectronEfficiency
  ElectronIsolation
  
 
  MuonEfficiency
  MuonIsolation

  NeutrinoFilter

  MissingET
  GenMissingET
  GenJetFinder
  FastJetFinder
  JetEnergyScale
 
  GenJetFlavorAssociation
  JetFlavorAssociation

  GenBTagging
  BTagging

  GenCTagging
  CTagging

  GenTauTagging
  TauTagging

  UniqueObjectFinder

  ScalarHT
  
  TreeWriter
}

# exclude since it adds also an entry to Jet->BTag
# TrackCountingBTagging

#################################
# Propagate particles in cylinder
#################################

module ParticlePropagator ParticlePropagator {
  set InputArray Delphes/stableParticles

  set OutputArray stableParticles
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

# FCC-eh
# radius of the magnetic field coverage, in m
  set Radius 1.36
  # half-length of the magnetic field coverage, in m
  set HalfLength 5.24

  # magnetic field
  set Bz 3.8

}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

# add EfficiencyFormula {efficiency formula as a function of eta and pt}

# tracking efficiency formula for charged hadrons
# set to 1 for full range (-0.2 eta already subtracted from edges)
  set EfficiencyFormula { (pt <= 0.5) * (0.00) + (pt > 0.5) *( eta <= 5.2 && eta >= -5.0)  * (1.0) + \
				(eta > 5.2 || eta < -5.0 ) * (0.00)}
}
##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for electrons
# set to 1
  set EfficiencyFormula {   (pt <= 0.5) * (0.00) + (pt > 0.5) *( (eta <=5.2 && eta >= -5.0) ) * (1.0) + \
			       ((eta > 5.2 || eta < -5.0)) * (0.00)}
 }
##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
 set InputArray ParticlePropagator/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for muons
# set to 1
  set EfficiencyFormula {    (pt <= 0.5) * (0.00) + (pt > 0.5) *(abs(eta) <= 4)  * (1.0) + \
				(abs(eta) > 4) * (0.00)}
}


########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons

  # set ResolutionFormula {resolution formula as a function of eta and pt}
### 17.3.17 Peter's values
  set ResolutionFormula {  (pt <= 0.1 ) * (0.00) + ( pt > 0.1   && eta < 5.2 && eta > -5.0)*sqrt( pt^2*5.e-5^2) + \ 
                            (eta > 5.2 || eta < -5.0) * (0.00)} }

  # resolution formula for charged hadrons
# taken from FCC but adjust to ep
#the abs(eta) <= 2.0 should be the "central best" and abs(eta) >
#2.0  and  <= 4.9 the fwd    /    <= 4.3 the bwd    (LHeC)
#        and abs(eta) > 2.0 and <= 5.2 fwd   /   <= 5.0   the bwd    (FCC)
#
#  set ResolutionFormula {                  (abs(eta) <= 2.0) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
#                                           (abs(eta) <= 2.0) * (pt > 1.0   && pt <= 1.0e1) * (0.01) + \
#                                           (abs(eta) <= 2.0) * (pt > 1.0e1 && pt <= 2.0e2) * (0.03) + \
#                                           (abs(eta) <= 2.0) * (pt > 2.0e2)                * (0.05) + \
#                         (abs(eta) > 2.0 && ( eta <= 5.2 && eta >= -5.0 ) ) * (pt > 0.1   && pt <= 1.0)   * (0.03) + \
#                         (abs(eta) > 2.0 && ( eta <= 5.2 && eta >= -5.0 ) ) * (pt > 1.0   && pt <= 1.0e1) * (0.02) + \
#                         (abs(eta) > 2.0 && ( eta <= 5.2 && eta >= -5.0 ) ) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
#                        (abs(eta) > 2.0 && ( eta <= 5.2 && eta >= -5.0 ) ) * (pt > 2.0e2)                * (0.05)}
#}


#################################
# Energy resolution for electrons
#################################

module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons


### 17.3.17 Peter's values
  set ResolutionFormula {  (pt <= 0.1 ) * (0.00) + ( pt > 0.1   && eta < 5.2 && eta > -5.0)*sqrt( pt^2*5.e-5^2) + \ 
                            (eta > 5.2 || eta < -5.0) * (0.00)} }

# 22.12.16 Uta
## take value from ECAL WITHOUT some energy ranges
## usually there were energy ranges given

# set ResolutionFormula {                  ( eta <= 5.5 && eta > 3.0  && energy > 0.1 ) * sqrt(energy*0.1^2 + energy^2*0.01^2)  + \
#                             ( eta <= 3.0 && eta >= -2.7 && energy > 0.1 ) * sqrt(energy*0.09^2 + energy^2*0.02^2)         +
#     \                            ( eta < -2.7 && eta >= -5.2 && energy > 0.1  ) * sqrt(energy*0.1^2 + energy^2*0.01^2) }
#
#}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons

  # set ResolutionFormula {resolution formula as a function of eta and pt}

### 17.3.17 Peter's values
  set ResolutionFormula {  (pt <= 0.1 ) * (0.00) + ( pt > 0.1   && eta < 5.2 && eta > -5.0)*sqrt( pt^2*5.e-5^2) + \ 
                            (eta > 5.2 || eta < -5.0) * (0.00)} }


  # resolution formula for muons
# taken from FCC
#   set ResolutionFormula {                  (abs(eta) <= 0.5) * (pt > 0.1   && pt <= 5.0)   * (0.02) + \
#                                           (abs(eta) <= 0.5) * (pt > 5.0   && pt <= 1.0e2) * (0.015) + \
#                                           (abs(eta) <= 0.5) * (pt > 1.0e2 && pt <= 2.0e2) * (0.03) + \
#                                           (abs(eta) <= 0.5) * (pt > 2.0e2)                * (0.05 + pt*1.e-4) + \
#                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 0.1   && pt <= 5.0)   * (0.03) + \
#                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 5.0   && pt <= 1.0e2) * (0.02) + \
#                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 1.0e2 && pt <= 2.0e2) * (0.04) + \
#                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05 + pt*1.e-4) + \
#                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 5.0)   * (0.04) + \
#                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 5.0   && pt <= 1.0e2) * (0.035) + \
#                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 1.0e2 && pt <= 2.0e2) * (0.05) + \
#                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 2.0e2)                * (0.05 + pt*1.e-4)}
#}


##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
#  add InputArray ElectronEnergySmearing/electrons
# change here to the momentum smearing
  add InputArray ElectronMomentumSmearing/electrons
  add InputArray MuonMomentumSmearing/muons
  set OutputArray tracks
}
###################
# Impact Parameters
###################

module ImpactParameterSmearing ImpactParameterSmearing {
  set TrackInputArray TrackMerger/tracks
  set OutputArray tracks

# set ResolutionFormula {resolution formula as a function of eta and pt}
# taken from FCC example - and about achieved at ATLAS Run-2
# absolute impact parameter smearing formula (in mm) as a function of pt and eta
#set ResolutionFormula {(pt > 0.1  && pt <= 5.0)   * (0.010) + \
#                       (pt > 5.0)                 * (0.005)}
#
#}
# value from Peter 5mu m x 20 mu m /( p sin theta ^3/2)
set ResolutionFormula { (pt > 0.5  )  * ( sqrt( (0.0001/pt)^2/(sin(2*atan(exp(-eta))))^3 ) )  }

##

##################
# LHeC Calorimeter
##################

module Calorimeter Calorimeter {
  set ParticleInputArray ParticlePropagator/stableParticles
## 17.3.17 use also TrackMerger/tracks
  set TrackInputArray TrackMerger/tracks
#  set TrackInputArray ImpactParameterSmearing/tracks

  set TowerOutputArray towers
  set PhotonOutputArray photons

  set EFlowTrackOutputArray eflowTracks
  set EFlowPhotonOutputArray eflowPhotons
  set EFlowNeutralHadronOutputArray eflowNeutralHadrons

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # use a delta_phi=0.025 for 2*phi range == 252 cells
  # use a delta_eta-0.025 for eta = -5 to 5 == 400 cells
  # Uta 22.12.16 use eta range from Peter Kostka -4.7 - 5.2
  set PhiBins {}
  for {set i -126} {$i <= 126} {incr i} {
      add PhiBins [expr {$i * $pi/126.0}]
  }
  foreach eta { -5.600 -5.575 -5.550 -5.525 -5.500 -5.475 -5.450 -5.425 -5.400 -5.375 -5.350 -5.325 -5.300 -5.275 -5.250 -5.225 -5.200 -5.175 -5.150 -5.125 -5.100 -5.075 -5.050 -5.025 -5.000  -4.975 -4.950 -4.925 -4.900 -4.875 -4.850 -4.825 -4.800 -4.775 -4.750 -4.725 -4.700 -4.675 -4.650 -4.625 -4.600 -4.575 -4.550 -4.525 -4.500 -4.475 -4.450 -4.425 -4.400 -4.375 -4.350 -4.325 -4.300 -4.275 -4.250 -4.225 -4.200 -4.175 -4.150 -4.125 -4.100 -4.075 -4.050 -4.025 -4.000 -3.975 -3.950 -3.925 -3.900 -3.875 -3.850 -3.825 -3.800 -3.775 -3.750 -3.725 -3.700 -3.675 -3.650 -3.625 -3.600 -3.575 -3.550 -3.525 -3.500 -3.475 -3.450 -3.425 -3.400 -3.375 -3.350 -3.325 -3.300 -3.275 -3.250 -3.225 -3.200 -3.175 -3.150 -3.125 -3.100 -3.075 -3.050 -3.025 -3.000 -2.975 -2.950 -2.925 -2.900 -2.875 -2.850 -2.825 -2.800 -2.775 -2.750 -2.725 -2.700 -2.675 -2.650 -2.625 -2.600 -2.575 -2.550 -2.525 -2.500 -2.475 -2.450 -2.425 -2.400 -2.375 -2.350 -2.325 -2.300 -2.275 -2.250 -2.225 -2.200 -2.175 -2.150 -2.125 -2.100 -2.075 -2.050 -2.025 -2.000 -1.975 -1.950 -1.925 -1.900 -1.875 -1.850 -1.825 -1.800 -1.775 -1.750 -1.725 -1.700 -1.675 -1.650 -1.625 -1.600 -1.575 -1.550 -1.525 -1.500 -1.475 -1.450 -1.425 -1.400 -1.375 -1.350 -1.325 -1.300 -1.275 -1.250 -1.225 -1.200 -1.175 -1.150 -1.125 -1.100 -1.075 -1.050 -1.025 -1.000 -0.975 -0.950 -0.925 -0.900 -0.875 -0.850 -0.825 -0.800 -0.775 -0.750 -0.725 -0.700 -0.675 -0.650 -0.625 -0.600 -0.575 -0.550 -0.525 -0.500 -0.475 -0.450 -0.425 -0.400 -0.375 -0.350 -0.325 -0.300 -0.275 -0.250 -0.225 -0.200 -0.175 -0.150 -0.125 -0.100 -0.075 -0.050 -0.025  0.000  0.025  0.050  0.075  0.100  0.125  0.150  0.175  0.200  0.225  0.250  0.275  0.300  0.325  0.350  0.375  0.400  0.425  0.450  0.475  0.500  0.525  0.550  0.575  0.600  0.625  0.650  0.675  0.700  0.725  0.750  0.775  0.800  0.825  0.850  0.875  0.900  0.925  0.950  0.975  1.000  1.025  1.050  1.075  1.100  1.125  1.150  1.175  1.200  1.225  1.250  1.275  1.300  1.325  1.350  1.375  1.400  1.425  1.450  1.475  1.500  1.525  1.550  1.575  1.600  1.625  1.650  1.675  1.700  1.725  1.750  1.775  1.800  1.825  1.850  1.875  1.900  1.925  1.950  1.975  2.000  2.025  2.050  2.075  2.100  2.125  2.150  2.175  2.200  2.225  2.250  2.275  2.300  2.325  2.350  2.375  2.400  2.425  2.450  2.475  2.500  2.525  2.550  2.575  2.600  2.625  2.650  2.675  2.700  2.725  2.750  2.775  2.800  2.825  2.850  2.875  2.900  2.925  2.950  2.975  3.000  3.025  3.050  3.075  3.100  3.125  3.150  3.175  3.200  3.225  3.250  3.275  3.300  3.325  3.350  3.375  3.400  3.425  3.450  3.475  3.500  3.525  3.550  3.575  3.600  3.625  3.650  3.675  3.700  3.725  3.750  3.775  3.800  3.825  3.850  3.875  3.900  3.925  3.950  3.975  4.000  4.025  4.050  4.075  4.100  4.125  4.150  4.175  4.200  4.225  4.250  4.275  4.300  4.325  4.350  4.375  4.400  4.425  4.450  4.475  4.500  4.525  4.550  4.575  4.600  4.625  4.650  4.675  4.700  4.725  4.750  4.775  4.800  4.825  4.850  4.875  4.900  4.925  4.950  4.975  5.00  5.025  5.050  5.075  5.100  5.125  5.150  5.175  5.200 5.225  5.250  5.275  5.300  5.325  5.350  5.375  5.400  5.425  5.450  5.475  5.500  5.525  5.550  5.575  5.600  5.625  5.650  5.675  5.700  5.725  5.750  5.775  5.800  5.825  5.850  5.875  5.900  5.925  5.950  5.975  6.000
  } {
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
  add EnergyFraction {3122} {0.3 0.7}

## no minimum p range given (Uta 22.12.16) - update with max 16.3.17 improved constant terms
## sigma_E = E* (a/sqrt(E) + b) = sqrt ( E*a^2 E +E^2*b^2 )
## a = sampling term, b = constant term 
# ECAL eta : 5.0 - 2.7 - -2.1 - -4.4 (subtract -0.1 in eta as safety distance to edges) 
    set ECalResolutionFormula {                  ( eta <= 5.5 && eta > 3.0  && energy > 0.1 ) * sqrt(energy*0.1^2 + energy^2*0.01^2)  + \
                             ( eta <= 3.0 && eta >= -2.7 && energy > 0.1 ) * sqrt(energy*0.09^2 + energy^2*0.01^2)         +
     \                            ( eta < -2.7 && eta >= -5.2 && energy > 0.1  ) * sqrt(energy*0.1^2 + energy^2*0.01^2) }
    
#### HCAL eta : 5.1 - 2.1 - -1.7 - -4.6 (subtract -0.1 in eta as safety distance to edges) 
#set HCalResolutionFormula {resolution formula as a function of eta and energy}

    set HCalResolutionFormula {  ( eta >= 5.9 && eta > 2.5 ) * sqrt(energy*0.3^2 + energy^2*0.01^2)  + \
				     ( eta <=2.5 && eta >= -2.2 ) * sqrt(energy*0.4^2 + energy^2*0.01^2)         + \   
	( eta < -2.2 && eta >= -5.5 ) * sqrt(energy*0.4^2 + energy^2*0.02^2) }

### HCal from Tanaka-san LHEC Hbb study
#set ECalResolutionFormula {                  (abs(eta) <= 3.0) * sqrt(energy^2*0.007^2 + energy*0.07^2 + 0.35^2)  + \
#                             (abs(eta) > 3.0 && abs(eta) <= 4.0) * sqrt(energy^2*0.02^2 + energy*0.2^2)           +
# \                             (abs(eta) > 4.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.10^2 + energy*0.4^2)}

# uta gues inspired by ILC particle flow resolution
  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
#  set HCalResolutionFormula {                  (abs(eta) <= 2.0) * sqrt(energy^2*0.030^2 + energy*0.35^2) + \
#                             (abs(eta) > 2.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.050^2 + energy*0.45^2)}


    
}
####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
# exchange UK:  25.7.16 change back to all eflowtracks!
#  add InputArray ImpactParameterSmearing/tracks
  add InputArray Calorimeter/eflowTracks
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
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
#  set InputArray Calorimeter/eflowPhotons
# 21.7.16 use photons
  set InputArray Calorimeter/photons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for photons - same as for  electrons
# expand and set flat
 # set to 1
  set EfficiencyFormula {  (pt <= 0.5) * (0.00) + (pt > 0.5) *( (eta <=5.4 && eta >= -5.1) ) * (1.0) + \
			       ((eta > 5.4 || eta < -5.1)) * (0.00)} 

} 

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
#  set IsolationInputArray EFlowMerger/eflow
 	set IsolationInputArray EFlowFilter/eflow

  set OutputArray photons

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.1
}
#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray Calorimeter/eflowTracks
  set OutputArray electrons
  set Invert true
    add PdgCode {11}
    add PdgCode {-11}
}



#################
# Muon filter
#################

module PdgCodeFilter MuonFilter {
  set InputArray Calorimeter/eflowTracks
  set OutputArray muons
  set Invert true
    add PdgCode {13}
    add PdgCode {-13}
}


#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
#  set InputArray ElectronEnergySmearing/electrons
  set InputArray ElectronFilter/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
# set to 1 for eta<4.7 1 degree 
   set EfficiencyFormula {   (pt <= 0.5) * (0.00) + (pt > 0.5) *( (eta <=5.4 && eta >= -5.1) ) * (1.0) + \
			       ((eta > 5.4 || eta < -5.1)) * (0.00)} 
# set EfficiencyFormula {                 (abs(eta) <= 4.7)  * (1.0) + \
#                         (abs(eta) > 4.7)    * (0.00)}
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
#  set IsolationInputArray EFlowMerger/eflow
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray electrons

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.1
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
  set InputArray MuonMomentumSmearing/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency as a function of eta and pt}

  # efficiency formula for muons
# set to 1 for eta<4
  set EfficiencyFormula {     (pt <= 0.5) * (0.00) + (pt > 0.5) *(abs(eta) <= 4)  * (1.0) + \
                         (abs(eta) > 4)                                   * (0.00)}
}

################
# Muon isolation
################

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
#  set IsolationInputArray EFlowMerger/eflow
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray muons

  set DeltaRMax 0.3

  set PTMin 0.5

  set PTRatioMax 0.1
}

###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}
#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {

# add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}


##################
# Scalar HT merger
##################

module Merger ScalarHT {
# add InputArray InputArray
  add InputArray UniqueObjectFinder/jets
  add InputArray UniqueObjectFinder/electrons
  add InputArray UniqueObjectFinder/photons
  add InputArray UniqueObjectFinder/muons
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
  set InputArray Delphes/stableParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 3.0
}

############
# Jet finder
############

module FastJetFinder FastJetFinder {
# ATLAS uses towers
#  set InputArray Calorimeter/towers
# CMS and we  use eflow
 set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.4

  set JetPTMin 3.0
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
########################
# GenJet Flavor Association
########################

module JetFlavorAssociation GenJetFlavorAssociation {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray GenJetFinder/jets
# 18.4.17 enhance pt cut and eta range -> eta cut in tagging modules
# change deltaR of 0.4 back to FCC-hh delta_R cut of 0.5
  set DeltaR 0.5
#  set PartonPTMin 0.5
  set PartonPTMin 1.0
 # set PartonEtaMax 4.0
 set PartonEtaMax 6.0

}
########################
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale/jets

# 18.4.17 enhance pt cut and eta range -> eta cut in tagging modules
  set DeltaR 0.5
#  set PartonPTMin 0.5
  set PartonPTMin 1.0
 # set PartonEtaMax 4.0
 set PartonEtaMax 6.0

}

#################
# Gen b-tagging
#################
module BTagging GenBTagging {
#  set PartonInputArray Delphes/partons
  set JetInputArray GenJetFinder/jets

# 18.4.17 taken from FCC hh but lower pt cut from 10 to 4
  set BitNumber 0	
 
  add EfficiencyFormula {0} {

  (pt <= 4.0)                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 4.0 && pt < 500)      * (0.001) + \
  (abs(eta) < 2.5)                   * (pt > 500.0 && pt < 20000.0) * (0.001)*(1.0 - pt/20000.) + \
  (abs(eta) < 2.5)                   * (pt > 20000.0)               * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 4.0 && pt < 500)      * (0.00075) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 20000.0) * (0.00075)*(1.0 - pt/20000.) + \
  (abs(eta) < 2.5 && abs(eta) < 4.0) * (pt > 20000.0)               * (0.000) + \
  (abs(eta) > 4.0) * (0.00)}

  add EfficiencyFormula {4} {

  (pt <= 4.0)                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 4.0 && pt < 500)      * (0.04) + \
  (abs(eta) < 2.5)                   * (pt > 500.0 && pt < 20000.0) * (0.04)*(1.0 - pt/20000.) + \
  (abs(eta) < 2.5)                   * (pt > 20000.0)               * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 4.0 && pt < 500)      * (0.03) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 20000.0) * (0.03)*(1.0 - pt/20000.) + \
  (abs(eta) < 2.5 && abs(eta) < 4.0) * (pt > 20000.0)               * (0.000) + \
  (abs(eta) > 4.0) * (0.00)}

  add EfficiencyFormula {5} {

  (pt <= 4.0)                                                       * (0.00) +
  (abs(eta) < 2.5)                    * (pt > 4.0 && pt < 500)      * (0.85) + 
  (abs(eta) < 2.5)                    * (pt > 500.0 && pt < 20000.0) * (0.85)*(1.0 - pt/20000.) + 
  (abs(eta) < 2.5)                    * (pt > 20000.0)               * (0.000) + 
  (abs(eta) >= 2.5 && abs(eta) < 4.0) * (pt > 4.0 && pt < 500)      * (0.64) + 
  (abs(eta) >= 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 20000.0) * (0.64)*(1.0 - pt/20000.) + 
  (abs(eta) <= 2.5 && abs(eta) < 4.0) * (pt > 20000.0)               * (0.000) + 
  (abs(eta) >= 4.0) * (0.00)}

}

#################
# Rec b-tagging
#################
module BTagging BTagging {
#  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

# 18.4.17 taken from FCC hh but lower pt cut from 10 to 4
  set BitNumber 0	
 
  add EfficiencyFormula {0} {

  (pt <= 4.0)                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 4.0 && pt < 500)      * (0.001) + \
  (abs(eta) < 2.5)                   * (pt > 500.0 && pt < 20000.0) * (0.001)*(1.0 - pt/20000.) + \
  (abs(eta) < 2.5)                   * (pt > 20000.0)               * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 4.0 && pt < 500)      * (0.00075) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 20000.0) * (0.00075)*(1.0 - pt/20000.) + \
  (abs(eta) < 2.5 && abs(eta) < 4.0) * (pt > 20000.0)               * (0.000) + \
  (abs(eta) > 4.0) * (0.00)}

  add EfficiencyFormula {4} {

  (pt <= 4.0)                       * (0.00) +
  (abs(eta) < 2.5)                   * (pt > 4.0 && pt < 500)      * (0.04) + \
  (abs(eta) < 2.5)                   * (pt > 500.0 && pt < 20000.0) * (0.04)*(1.0 - pt/20000.) + \
  (abs(eta) < 2.5)                   * (pt > 20000.0)               * (0.000) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 4.0 && pt < 500)      * (0.03) + \
  (abs(eta) > 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 20000.0) * (0.03)*(1.0 - pt/20000.) + \
  (abs(eta) < 2.5 && abs(eta) < 4.0) * (pt > 20000.0)               * (0.000) + \
  (abs(eta) > 4.0) * (0.00)}

  add EfficiencyFormula {5} {

  (pt <= 4.0)                                                       * (0.00) +
  (abs(eta) < 2.5)                    * (pt > 4.0 && pt < 500)      * (0.85) + 
  (abs(eta) < 2.5)                    * (pt > 500.0 && pt < 20000.0) * (0.85)*(1.0 - pt/20000.) + 
  (abs(eta) < 2.5)                    * (pt > 20000.0)               * (0.000) + 
  (abs(eta) >= 2.5 && abs(eta) < 4.0) * (pt > 4.0 && pt < 500)      * (0.64) + 
  (abs(eta) >= 2.5 && abs(eta) < 4.0) * (pt > 500.0 && pt < 20000.0) * (0.64)*(1.0 - pt/20000.) + 
  (abs(eta) <= 2.5 && abs(eta) < 4.0) * (pt > 20000.0)               * (0.000) + 
  (abs(eta) >= 4.0) * (0.00)}

}
#################
# Gen  c-tagging
#################

module BTagging GenCTagging {
#  set PartonInputArray Delphes/partons
  set JetInputArray GenJetFinder/jets


  set BitNumber 1

  add EfficiencyFormula {0} {

  (pt <= 4.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 4.0) * (0.01) + \
  (abs(eta) > 4.0) * (pt > 4.0) * (0.00)}

  add EfficiencyFormula {4} {

  (pt <= 4.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 4.0) * (0.10) + \
  (abs(eta) > 4.0) * (pt > 4.0) * (0.00)}

  add EfficiencyFormula {5} {

  (pt <= 4.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 4.0) * (0.25) + \
  (abs(eta) > 4.0) * (pt > 4.0) * (0.00)}

}

#################
# Rec  c-tagging
#################

module BTagging CTagging {
#  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets	


  set BitNumber 1

  add EfficiencyFormula {0} {

  (pt <= 4.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 4.0) * (0.01) + \
  (abs(eta) > 4.0) * (pt > 4.0) * (0.00)}

  add EfficiencyFormula {4} {

  (pt <= 4.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 4.0) * (0.10) + \
  (abs(eta) > 4.0) * (pt > 4.0) * (0.00)}

  add EfficiencyFormula {5} {

  (pt <= 4.0)     * (0.00) +
  (abs(eta) < 4.0) * (pt > 4.0) * (0.25) + \
  (abs(eta) > 4.0) * (pt > 4.0) * (0.00)}

}
module TauTagging TauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5

  set TauPTMin 0.5

  set TauEtaMax 6.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} { (pt <= 0.5) * (0.00) + (pt > 0.5) * (0.4) }
}

module TauTagging GenTauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray GenJetFinder/jets

  set DeltaR 0.4

  set TauPTMin 0.5

  set TauEtaMax 6.0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} { (pt <= 0.5) * (0.00) + (pt > 0.5) * (0.4) }
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

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle

  add Branch TrackMerger/tracks Track Track
  #add Branch Calorimeter/towers Tower Tower

  #add Branch Calorimeter/eflowTracks EFlowTrack Track
  #add Branch Calorimeter/eflowPhotons EFlowPhoton Tower
  #add Branch Calorimeter/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch GenJetFinder/jets GenJet Jet
  add Branch UniqueObjectFinder/jets Jet Jet
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/photons Photon Photon
  add Branch UniqueObjectFinder/muons Muon Muon
  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}
