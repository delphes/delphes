######################################################################################################################
# EIC detector model
# based on parameters from EIC detector matrix from EIC yellow report https://physdiv.jlab.org/DetectorMatrix/
# as well as on assumptions on calorimeter granularity and tracking efficiency (not specified in handbook).
#Berkeley all-silicon tracker 3.0 T. Taken from slides from Rey Cruz-Torres https://indico.bnl.gov/event/7913/
# email: miguel.arratia@ucr.edu
#######################################################################################################################

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
  TrackSmearing

  ECal
  HCal

  Calorimeter
  EFlowMerger
  EFlowFilter

  PhotonEfficiency
  PhotonIsolation

  ElectronFilter
  ElectronEfficiency
  ElectronIsolation

  ChargedHadronFilter
  MissingET

  NeutrinoFilter
  GenJetFinder
  GenMissingET

  FastJetFinder

  JetEnergyScale

  JetFlavorAssociation
  GenJetFlavorAssociation

  UniqueObjectFinder

  ScalarHT

  TrackCountingBTagging

  mRICH
  barrelDIRC
  dualRICH_aerogel
  dualRICH_c2f6

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
    set Radius 0.8
    # half-length of the magnetic field coverage, in m
    set HalfLength 1.00
    # magnetic field
    set Bz 3.0
}


####################################
# Common Tracking Efficiency Model
####################################
#Dummy efficiency (100%). Leaving structure to show how tracking dependent on pt and eta can be incorporated)
#

#From EIC YR detector matrix:
#Minimum pT for B = 3 T:
#150 MeV/c for -3.0 < eta < -2.5
#220 MeV/c for -2.5 < eta < -2.0
#160 MeV/c for -2.0 < eta < -1.5
#300 MeV/c for -1.5 < eta < -1.0
#(For B = 3T: minimum pT = 400 MeV/c with 90% acceptance (similar for pi and K))

set CommonTrackingEfficiency {

    (abs(eta) <= 1.0) * (pt > 0.400)                     * (1.0) +
    (abs(eta) > 1.0 && abs(eta) <= 1.5) * (pt > 0.300)   * (1.0) +
    (abs(eta) > 1.5 && abs(eta) <= 2.0) * (pt > 0.160)   * (1.0) +
    (abs(eta) > 2.0 && abs(eta) <= 2.5) * (pt > 0.220)   * (1.0) +
    (abs(eta) > 2.5 && abs(eta) <= 3.5) * (pt > 0.150)   * (1.0) +
    (abs(eta) > 3.5 && abs(eta) <= 4.0) * (pt > 0.150)   * (1.0) +

    (abs(eta) > 4.0)                                                  * (0.00)+
    0.0
}

#Berkeley all-silicon tracker 3.0 T. Taken from slides from Rey Cruz-Torres https://indico.bnl.gov/event/7913/
set CommonTrackingResolution {
    (abs(eta)<=0.5)                  * (sqrt( (3.69e-3)^2 + (pt*cosh(eta)*1.8e-4)^2  ) )  +
    (abs(eta)<=1.0 && abs(eta)>0.5) * (sqrt( (4.28e-3)^2 + (pt*cosh(eta)*1.6e-4)^2   ) )  +
    (abs(eta)<=1.5 && abs(eta)>1.0) * (sqrt( (4.27e-3)^2 + (pt*cosh(eta)*1.6e-4)^2   ) )  +
    (abs(eta)<=2.0 && abs(eta)>1.5) * (sqrt( (4.62e-3)^2 + (pt*cosh(eta)*1.2e-4)^2   ) )  +
    (abs(eta)<=2.5 && abs(eta)>2.0) * (sqrt( (7.19e-3)^2 + (pt*cosh(eta)*1.8e-4)^2   ) )  +
    (abs(eta)<=3.0 && abs(eta)>2.5) * (sqrt( (1.34e-2)^2 + (pt*cosh(eta)*3.9e-4)^2   ) )  +
    (abs(eta)<=3.5 && abs(eta)>3.0) * (sqrt( (2.43e-2)^2 + (pt*cosh(eta)*1.03e-3)^2  ) )  +
    (abs(eta)<=4.0 && abs(eta)>3.5) * (sqrt( (4.56e-2)^2 + (pt*cosh(eta)*2.95e-3)^2  ) )
}


####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  set EfficiencyFormula $CommonTrackingEfficiency
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons
  set EfficiencyFormula $CommonTrackingEfficiency

}

##############################
# Muon tracking efficiency
##############################
module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons
  set EfficiencyFormula $CommonTrackingEfficiency

}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons
  set ResolutionFormula  $CommonTrackingResolution
}

###################################
# Momentum resolution for muons
###################################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons
  set ResolutionFormula $CommonTrackingResolution
}



###################################
# Momentum resolution for electrons
###################################
module MomentumSmearing ElectronMomentumSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons
  set ResolutionFormula $CommonTrackingResolution
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

################################
# Track impact parameter smearing
################################


module TrackSmearing TrackSmearing {
  set InputArray TrackMerger/tracks
#  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
  set OutputArray tracks
#  set ApplyToPileUp true
  # magnetic field
  set Bz 3.0
  set PResolutionFormula { 0.0 }
  set CtgThetaResolutionFormula { 0.0 }
  set PhiResolutionFormula { 0.0 }
#Berkeley all-silicon tracker 3.0 T. Taken from slides from Rey Cruz-Torres https://indico.bnl.gov/event/7913/
  set D0ResolutionFormula "
    (abs(eta)<=0.5)                   * (sqrt( (0.0048)^2 +   (0.025/pt)^2   ) )  +
    (abs(eta)<=1.0 && abs(eta)>0.5)   * (sqrt( (0.0045)^2 +   (0.029/pt)^2   ) )  +
    (abs(eta)<=1.5 && abs(eta)>1.0)   * (sqrt( (0.0055)^2 +   (0.033/pt)^2   ) )  +
    (abs(eta)<=2.0 && abs(eta)>1.5)   * (sqrt( (0.0055)^2 +   (0.039/pt)^2   ) )  +
    (abs(eta)<=2.5 && abs(eta)>2.0)   * (sqrt( (0.0095)^2 +   (0.045/pt)^2   ) )
  "


  set DZResolutionFormula "
    (abs(eta)<=0.5)                   * (sqrt( (0.0032)^2 +   (0.027/pt)^2   ) )  +
    (abs(eta)<=1.0 && abs(eta)>0.5)   * (sqrt( (0.0038)^2 +   (0.037/pt)^2   ) )  +
    (abs(eta)<=1.5 && abs(eta)>1.0)   * (sqrt( (0.0059)^2 +   (0.056/pt)^2   ) )  +
    (abs(eta)<=2.0 && abs(eta)>1.5)   * (sqrt( (0.0087)^2 +   (0.108/pt)^2   ) )  +
    (abs(eta)<=2.5 && abs(eta)>2.0)   * (sqrt( (0.0198)^2 +   (0.207/pt)^2   ) )
  "

}


#############
#   ECAL
#############

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackSmearing/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true
  set EnergyMin 0.050
  #does not seem possible to set minimum dependent on eta as spec in the YR.


  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

  # Granularity is not discussed in EIC detector handbook.
  ##BARREL
  #assume 0.1 x 0.1 (real cell size will be smaller, so this is to represent some cluster)

    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    for {set i -10} {$i <=10} {incr i} {
	set eta [expr {$i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }

    ## Coverage is -3.5, -1.0 , and +1.0 to 3.5.
   ## assume 0.1 x 0.1 (real cell size will be smaller, so this is to represent some cluster)
    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }

    for {set i 1} {$i <=26} {incr i} {
	set eta [expr {-3.6 + $i*0.1}]
	add EtaPhiBins $eta $PhiBins
    }
    for {set i 1} {$i <=26} {incr i} {
	set eta [expr {0.9 + $i*0.1 }]
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
 # add EnergyFraction {310} {0.3}
 # add EnergyFraction {3122} {0.3}

  set ResolutionFormula {          (eta <= -2.0 && eta>-3.5)                          * sqrt(energy^2*0.01^2 + energy*0.025^2 + 0.01^2)+ \
                 		   (eta <= -1.0 && eta>-2.0 )                         * sqrt(energy^2*0.02^2 + energy*0.08^2 + 0.02^2 )+ \
				   (eta <= 1.0  && eta> -1.0 )                        * sqrt(energy^2*0.03^2 + energy*0.14^2 + 0.02^2 )+ \
				   (eta <= 3.5  &&  eta>1.0 )                         * sqrt(energy^2*0.02^2 + energy*0.12^2 + 0.02^2)}

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

  ##Assumes noise 100 MeV per tower.
  set EnergyMin 0.5
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    for {set i -10} {$i <=10} {incr i} {
	set eta [expr {$i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }

    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }

    for {set i 1} {$i <=26} {incr i} {
	set eta [expr {-3.6 + $i*0.1 }]
	add EtaPhiBins $eta $PhiBins
    }
    for {set i 1} {$i <=26} {incr i} {
	set eta [expr {0.9 + $i*0.1 }]
	add EtaPhiBins $eta $PhiBins
    }


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
  #add EnergyFraction {310} {0.7}
  #add EnergyFraction {3122} {0.7}

  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {    (eta <= -1.0 && eta>-3.5)                       * sqrt(energy^2*0.10^2 + energy*0.50^2)+
                             (eta <= 1.0 && eta>-1.0 )                       * sqrt(energy^2*0.10^2 + energy*1.00^2)+
                             (eta <= 3.5 && eta>1.0 )                       * sqrt(energy^2*0.10^2 + energy*0.50^2)
  }

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
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray ECal/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for photons
    set EfficiencyFormula { 1}
}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray photons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
}


#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray ElectronFilter/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
    set EfficiencyFormula {1}
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray electrons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
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
  add InputArray UniqueObjectFinder/jets
  add InputArray UniqueObjectFinder/electrons
  add InputArray UniqueObjectFinder/photons

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
  set ParameterR 1.0

  set JetPTMin 3.0
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
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.0

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeTrimming 1
  set RTrim 0.4
  set PtFracTrim 0.20
  #set PtFracTrim 0.05

  set ComputePruning 1
  set ZcutPrun 0.1
  set RcutPrun 0.5
  set RPrun 0.8

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

  set JetPTMin 3.0}





##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder/jets
  set OutputArray jets

  # scale formula for jets (do not apply it)
  set ScaleFormula {1.0}
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
  set PartonPTMin 4.0
  set PartonEtaMax 4.0

}

module JetFlavorAssociation GenJetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray GenJetFinder/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 4.0

}



#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

module UniqueObjectFinder UniqueObjectFinder {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray JetEnergyScale/jets jets
}

############################
# b-tagging (track counting)
############################

module TrackCountingBTagging TrackCountingBTagging {
    set JetInputArray JetEnergyScale/jets
    set TrackInputArray HCal/eflowTracks
    set BitNumber 0
    # maximum distance between jet and track
    set DeltaR 0.5
    # minimum pt of tracks
    set TrackPtMin 1.0
    # maximum transverse impact parameter (in mm)
    set TrackIPMax 3
    # minimum ip significance for the track to be counted
    set SigMin 2.0
    set Use3D true
    # alternate setting for 2D IP (default)
    #  set SigMin 1.3
    #  set Use3D false
    # minimum number of tracks (high efficiency n=2, high purity n=3)
    #set Ntracks 3

}

##################
# Particle ID Systems
##################

# Efficiency maps determined from EIC PID Code, https://gitlab.com/preghenella/pid/


#My name is "mRICH Pixel Size =3x3 (mm^{2})" and I am described as follows:
#    Momentum coverage for k/pi separation=  [3 to 10 (GeV/c) ]
#    Since mRICH detector system consists of array of identical shoebox-like modules, in principle, 
#    the mRICH pid performance is invariant relative to the tracking polar angle but the momentum of the 
#    charge particle and the entrance point location on the front face of a given mRICH module.
#    Assumed time precision = 1 ns
#    Assumed track resolution = 0.00175 mrad

module IdentificationMap mRICH { 
  set InputArray TrackSmearing/tracks
  set OutputArray tracks 

   #--- kaons ---

    add EfficiencyFormula {321} {321} {
      (eta<-3.50 || eta>=-1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
      (-3.50 <= eta && eta < -1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 3.10) * (0.000000) + 
      (-3.50 <= eta && eta < -1.00) * (3.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.998780) + 
      (-3.50 <= eta && eta < -1.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.979055) + 
      (-3.50 <= eta && eta < -1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.957718) + 
      (-3.50 <= eta && eta < -1.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.936379) + 
      (-3.50 <= eta && eta < -1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.915859) + 
      (-3.50 <= eta && eta < -1.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.894273) + 
      (-3.50 <= eta && eta < -1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.872402) + 
      (-3.50 <= eta && eta < -1.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.850872) + 
      (-3.50 <= eta && eta < -1.00) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.826813) + 
      (-3.50 <= eta && eta < -1.00) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.801204) + 
      (-3.50 <= eta && eta < -1.00) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.777961) + 
      (-3.50 <= eta && eta < -1.00) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.754789) + 
      (-3.50 <= eta && eta < -1.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000)} 
  
    add EfficiencyFormula {321} {-211} {
      (eta<-3.50 || eta>=-1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
      (-3.50 <= eta && eta < -1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.000713) + 
      (-3.50 <= eta && eta < -1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.021562) + 
      (-3.50 <= eta && eta < -1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.043027) + 
      (-3.50 <= eta && eta < -1.00) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.065594) + 
      (-3.50 <= eta && eta < -1.00) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.087964) + 
      (-3.50 <= eta && eta < -1.00) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.108606) + 
      (-3.50 <= eta && eta < -1.00) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.129309) + 
      (-3.50 <= eta && eta < -1.00) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.152113) + 
      (-3.50 <= eta && eta < -1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.169280) + 
      (-3.50 <= eta && eta < -1.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000)} 
  
  
    add EfficiencyFormula {321} {2212} {
      (eta<-3.50 || eta>=-1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
      (-3.50 <= eta && eta < -1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.000713) + 
      (-3.50 <= eta && eta < -1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.021562) + 
      (-3.50 <= eta && eta < -1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.043027) + 
      (-3.50 <= eta && eta < -1.00) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.065594) + 
      (-3.50 <= eta && eta < -1.00) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.087964) + 
      (-3.50 <= eta && eta < -1.00) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.108606) + 
      (-3.50 <= eta && eta < -1.00) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.129309) + 
      (-3.50 <= eta && eta < -1.00) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.152113) + 
      (-3.50 <= eta && eta < -1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.169280) + 
      (-3.50 <= eta && eta < -1.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
      (-1.00 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 
  
    # --- pions ---
  
    add EfficiencyFormula {-211} {321} {
      (eta<-3.50 || eta>=-1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
      (-3.50 <= eta && eta < -1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.000713) + 
      (-3.50 <= eta && eta < -1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.021562) + 
      (-3.50 <= eta && eta < -1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.043027) + 
      (-3.50 <= eta && eta < -1.00) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.065594) + 
      (-3.50 <= eta && eta < -1.00) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.087964) + 
      (-3.50 <= eta && eta < -1.00) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.108606) + 
      (-3.50 <= eta && eta < -1.00) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.129309) + 
      (-3.50 <= eta && eta < -1.00) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.152113) + 
      (-3.50 <= eta && eta < -1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.169280) + 
      (-3.50 <= eta && eta < -1.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 
  
    add EfficiencyFormula {-211} {-211} {
      (eta<-3.50 || eta>=-1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
      (-3.50 <= eta && eta < -1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 3.10) * (0.000000) + 
      (-3.50 <= eta && eta < -1.00) * (3.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.998780) + 
      (-3.50 <= eta && eta < -1.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.979055) + 
      (-3.50 <= eta && eta < -1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.957718) + 
      (-3.50 <= eta && eta < -1.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.936379) + 
      (-3.50 <= eta && eta < -1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.915859) + 
      (-3.50 <= eta && eta < -1.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.894273) + 
      (-3.50 <= eta && eta < -1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.872402) + 
      (-3.50 <= eta && eta < -1.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.850872) + 
      (-3.50 <= eta && eta < -1.00) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.826813) + 
      (-3.50 <= eta && eta < -1.00) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.801204) + 
      (-3.50 <= eta && eta < -1.00) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.777961) + 
      (-3.50 <= eta && eta < -1.00) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.754789) + 
      (-3.50 <= eta && eta < -1.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 
  
    # --- protons ---
  
    add EfficiencyFormula {2212} {321} {
      (eta<-3.50 || eta>=-1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
      (-3.50 <= eta && eta < -1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.000713) + 
      (-3.50 <= eta && eta < -1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.021562) + 
      (-3.50 <= eta && eta < -1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.043027) + 
      (-3.50 <= eta && eta < -1.00) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.065594) + 
      (-3.50 <= eta && eta < -1.00) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.087964) + 
      (-3.50 <= eta && eta < -1.00) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.108606) + 
      (-3.50 <= eta && eta < -1.00) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.129309) + 
      (-3.50 <= eta && eta < -1.00) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.152113) + 
      (-3.50 <= eta && eta < -1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.169280) + 
      (-3.50 <= eta && eta < -1.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 
  
    add EfficiencyFormula {2212} {2212} {
      (eta<-3.50 || eta>=-1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
      (-3.50 <= eta && eta < -1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 3.10) * (0.000000) + 
      (-3.50 <= eta && eta < -1.00) * (3.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.998780) + 
      (-3.50 <= eta && eta < -1.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.979055) + 
      (-3.50 <= eta && eta < -1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.957718) + 
      (-3.50 <= eta && eta < -1.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.936379) + 
      (-3.50 <= eta && eta < -1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.915859) + 
      (-3.50 <= eta && eta < -1.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.894273) + 
      (-3.50 <= eta && eta < -1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.872402) + 
      (-3.50 <= eta && eta < -1.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.850872) + 
      (-3.50 <= eta && eta < -1.00) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.826813) + 
      (-3.50 <= eta && eta < -1.00) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.801204) + 
      (-3.50 <= eta && eta < -1.00) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.777961) + 
      (-3.50 <= eta && eta < -1.00) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.754789) + 
      (-3.50 <= eta && eta < -1.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 


# Everything else is not ID'd at all.
  add EfficiencyFormula {0} {0} { 0.00 }

}


#My name is "Barrel DIRC TR=0.5 [mrad] dT=0.1 ns QE = 27 %" and I am described as follows:
#    Eta coverage =  [-1,1]
#    Assumed time precision = 0.1 ns
#    Assumed track resolution = 0.5 mrad
#    Assumed quantum efficiency of the MCP-PMT = 27%

module IdentificationMap barrelDIRC { 
  set InputArray TrackSmearing/tracks
  set OutputArray tracks 

  # --- kaons ---

  add EfficiencyFormula {321} {321} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-1.00 <= eta && eta < -0.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.70) * (0.999363) + 
    (-1.00 <= eta && eta < -0.90) * (5.70 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.977469) + 
    (-1.00 <= eta && eta < -0.90) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.953824) + 
    (-1.00 <= eta && eta < -0.90) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.930396) + 
    (-1.00 <= eta && eta < -0.90) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.908335) + 
    (-1.00 <= eta && eta < -0.90) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.885172) + 
    (-1.00 <= eta && eta < -0.90) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.861520) + 
    (-1.00 <= eta && eta < -0.90) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.838567) + 
    (-1.00 <= eta && eta < -0.90) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.817634) + 
    (-1.00 <= eta && eta < -0.90) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.793978) + 
    (-1.00 <= eta && eta < -0.90) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.769187) + 
    (-1.00 <= eta && eta < -0.90) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.745592) + 
    (-1.00 <= eta && eta < -0.90) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.721481) + 
    (-1.00 <= eta && eta < -0.90) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.698201) + 
    (-1.00 <= eta && eta < -0.90) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.676260) + 
    (-1.00 <= eta && eta < -0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.654166) + 
    (-1.00 <= eta && eta < -0.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.629819) + 
    (-1.00 <= eta && eta < -0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604030) + 
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.90 <= eta && eta < -0.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.60) * (0.999415) + 
    (-0.90 <= eta && eta < -0.80) * (5.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.978320) + 
    (-0.90 <= eta && eta < -0.80) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.954946) + 
    (-0.90 <= eta && eta < -0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.931253) + 
    (-0.90 <= eta && eta < -0.80) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.908742) + 
    (-0.90 <= eta && eta < -0.80) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.884654) + 
    (-0.90 <= eta && eta < -0.80) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.860913) + 
    (-0.90 <= eta && eta < -0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.837747) + 
    (-0.90 <= eta && eta < -0.80) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.815568) + 
    (-0.90 <= eta && eta < -0.80) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.791569) + 
    (-0.90 <= eta && eta < -0.80) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.767848) + 
    (-0.90 <= eta && eta < -0.80) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.747406) + 
    (-0.90 <= eta && eta < -0.80) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.723516) + 
    (-0.90 <= eta && eta < -0.80) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.699697) + 
    (-0.90 <= eta && eta < -0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.677096) + 
    (-0.90 <= eta && eta < -0.80) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.654458) + 
    (-0.90 <= eta && eta < -0.80) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.630047) + 
    (-0.90 <= eta && eta < -0.80) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603999) + 
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.80 <= eta && eta < -0.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999438) + 
    (-0.80 <= eta && eta < -0.70) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.978782) + 
    (-0.80 <= eta && eta < -0.70) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.953835) + 
    (-0.80 <= eta && eta < -0.70) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.929657) + 
    (-0.80 <= eta && eta < -0.70) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.907134) + 
    (-0.80 <= eta && eta < -0.70) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.882151) + 
    (-0.80 <= eta && eta < -0.70) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.858132) + 
    (-0.80 <= eta && eta < -0.70) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.835697) + 
    (-0.80 <= eta && eta < -0.70) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.813827) + 
    (-0.80 <= eta && eta < -0.70) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.790188) + 
    (-0.80 <= eta && eta < -0.70) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.765614) + 
    (-0.80 <= eta && eta < -0.70) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.742033) + 
    (-0.80 <= eta && eta < -0.70) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.720207) + 
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.698161) + 
    (-0.80 <= eta && eta < -0.70) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.675683) + 
    (-0.80 <= eta && eta < -0.70) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.653272) + 
    (-0.80 <= eta && eta < -0.70) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.628628) + 
    (-0.80 <= eta && eta < -0.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603794) + 
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.70 <= eta && eta < -0.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999386) + 
    (-0.70 <= eta && eta < -0.60) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.977686) + 
    (-0.70 <= eta && eta < -0.60) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.952568) + 
    (-0.70 <= eta && eta < -0.60) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.927697) + 
    (-0.70 <= eta && eta < -0.60) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.904394) + 
    (-0.70 <= eta && eta < -0.60) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.880318) + 
    (-0.70 <= eta && eta < -0.60) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.856544) + 
    (-0.70 <= eta && eta < -0.60) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.833257) + 
    (-0.70 <= eta && eta < -0.60) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.810933) + 
    (-0.70 <= eta && eta < -0.60) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.786986) + 
    (-0.70 <= eta && eta < -0.60) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.763520) + 
    (-0.70 <= eta && eta < -0.60) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.740261) + 
    (-0.70 <= eta && eta < -0.60) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.718074) + 
    (-0.70 <= eta && eta < -0.60) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.696228) + 
    (-0.70 <= eta && eta < -0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.674142) + 
    (-0.70 <= eta && eta < -0.60) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.651446) + 
    (-0.70 <= eta && eta < -0.60) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.625653) + 
    (-0.70 <= eta && eta < -0.60) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603457) + 
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.60 <= eta && eta < -0.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999412) + 
    (-0.60 <= eta && eta < -0.50) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.978142) + 
    (-0.60 <= eta && eta < -0.50) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.953041) + 
    (-0.60 <= eta && eta < -0.50) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.928653) + 
    (-0.60 <= eta && eta < -0.50) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.906049) + 
    (-0.60 <= eta && eta < -0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.882734) + 
    (-0.60 <= eta && eta < -0.50) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.858632) + 
    (-0.60 <= eta && eta < -0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.835440) + 
    (-0.60 <= eta && eta < -0.50) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.813533) + 
    (-0.60 <= eta && eta < -0.50) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.789624) + 
    (-0.60 <= eta && eta < -0.50) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.765618) + 
    (-0.60 <= eta && eta < -0.50) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.744122) + 
    (-0.60 <= eta && eta < -0.50) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.722552) + 
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.700157) + 
    (-0.60 <= eta && eta < -0.50) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.677085) + 
    (-0.60 <= eta && eta < -0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.654123) + 
    (-0.60 <= eta && eta < -0.50) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.629276) + 
    (-0.60 <= eta && eta < -0.50) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603840) + 
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.50 <= eta && eta < -0.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999451) + 
    (-0.50 <= eta && eta < -0.40) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979067) + 
    (-0.50 <= eta && eta < -0.40) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.958187) + 
    (-0.50 <= eta && eta < -0.40) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.937021) + 
    (-0.50 <= eta && eta < -0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.913728) + 
    (-0.50 <= eta && eta < -0.40) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.889767) + 
    (-0.50 <= eta && eta < -0.40) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.865592) + 
    (-0.50 <= eta && eta < -0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.842135) + 
    (-0.50 <= eta && eta < -0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.819614) + 
    (-0.50 <= eta && eta < -0.40) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.795484) + 
    (-0.50 <= eta && eta < -0.40) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.770150) + 
    (-0.50 <= eta && eta < -0.40) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.747397) + 
    (-0.50 <= eta && eta < -0.40) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.725362) + 
    (-0.50 <= eta && eta < -0.40) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.703515) + 
    (-0.50 <= eta && eta < -0.40) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.681967) + 
    (-0.50 <= eta && eta < -0.40) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.659657) + 
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.636152) + 
    (-0.50 <= eta && eta < -0.40) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604442) + 
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.40 <= eta && eta < -0.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999472) + 
    (-0.40 <= eta && eta < -0.30) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979617) + 
    (-0.40 <= eta && eta < -0.30) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.955630) + 
    (-0.40 <= eta && eta < -0.30) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.932152) + 
    (-0.40 <= eta && eta < -0.30) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.909173) + 
    (-0.40 <= eta && eta < -0.30) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.885363) + 
    (-0.40 <= eta && eta < -0.30) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.860582) + 
    (-0.40 <= eta && eta < -0.30) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.836574) + 
    (-0.40 <= eta && eta < -0.30) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.813641) + 
    (-0.40 <= eta && eta < -0.30) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.789231) + 
    (-0.40 <= eta && eta < -0.30) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.764564) + 
    (-0.40 <= eta && eta < -0.30) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.740178) + 
    (-0.40 <= eta && eta < -0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715912) + 
    (-0.40 <= eta && eta < -0.30) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.693289) + 
    (-0.40 <= eta && eta < -0.30) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.670721) + 
    (-0.40 <= eta && eta < -0.30) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647669) + 
    (-0.40 <= eta && eta < -0.30) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605441) + 
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.30 <= eta && eta < -0.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999491) + 
    (-0.30 <= eta && eta < -0.20) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979938) + 
    (-0.30 <= eta && eta < -0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.959867) + 
    (-0.30 <= eta && eta < -0.20) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.938816) + 
    (-0.30 <= eta && eta < -0.20) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.916434) + 
    (-0.30 <= eta && eta < -0.20) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.892350) + 
    (-0.30 <= eta && eta < -0.20) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.868247) + 
    (-0.30 <= eta && eta < -0.20) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.845108) + 
    (-0.30 <= eta && eta < -0.20) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.821653) + 
    (-0.30 <= eta && eta < -0.20) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.799669) + 
    (-0.30 <= eta && eta < -0.20) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.776530) + 
    (-0.30 <= eta && eta < -0.20) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.754283) + 
    (-0.30 <= eta && eta < -0.20) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.732042) + 
    (-0.30 <= eta && eta < -0.20) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.708971) + 
    (-0.30 <= eta && eta < -0.20) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.685968) + 
    (-0.30 <= eta && eta < -0.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.663406) + 
    (-0.30 <= eta && eta < -0.20) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.640318) + 
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604825) + 
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.20 <= eta && eta < -0.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999335) + 
    (-0.20 <= eta && eta < -0.10) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.976774) + 
    (-0.20 <= eta && eta < -0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.951314) + 
    (-0.20 <= eta && eta < -0.10) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.926371) + 
    (-0.20 <= eta && eta < -0.10) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.902907) + 
    (-0.20 <= eta && eta < -0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.878307) + 
    (-0.20 <= eta && eta < -0.10) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.855267) + 
    (-0.20 <= eta && eta < -0.10) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.831969) + 
    (-0.20 <= eta && eta < -0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.806045) + 
    (-0.20 <= eta && eta < -0.10) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.780629) + 
    (-0.20 <= eta && eta < -0.10) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.757736) + 
    (-0.20 <= eta && eta < -0.10) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.734043) + 
    (-0.20 <= eta && eta < -0.10) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.711369) + 
    (-0.20 <= eta && eta < -0.10) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.688232) + 
    (-0.20 <= eta && eta < -0.10) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.665130) + 
    (-0.20 <= eta && eta < -0.10) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.641972) + 
    (-0.20 <= eta && eta < -0.10) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605010) + 
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.10 <= eta && eta < -0.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.60) * (0.999436) + 
    (-0.10 <= eta && eta < -0.00) * (5.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.978964) + 
    (-0.10 <= eta && eta < -0.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.958369) + 
    (-0.10 <= eta && eta < -0.00) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.938389) + 
    (-0.10 <= eta && eta < -0.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.916194) + 
    (-0.10 <= eta && eta < -0.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.893059) + 
    (-0.10 <= eta && eta < -0.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.869718) + 
    (-0.10 <= eta && eta < -0.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.846143) + 
    (-0.10 <= eta && eta < -0.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.825154) + 
    (-0.10 <= eta && eta < -0.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.800818) + 
    (-0.10 <= eta && eta < -0.00) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.778480) + 
    (-0.10 <= eta && eta < -0.00) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.755699) + 
    (-0.10 <= eta && eta < -0.00) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.733208) + 
    (-0.10 <= eta && eta < -0.00) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.710395) + 
    (-0.10 <= eta && eta < -0.00) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.687940) + 
    (-0.10 <= eta && eta < -0.00) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.665437) + 
    (-0.10 <= eta && eta < -0.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.642619) + 
    (-0.10 <= eta && eta < -0.00) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605218) + 
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.00 <= eta && eta < 0.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999487) + 
    (-0.00 <= eta && eta < 0.10) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.980086) + 
    (-0.00 <= eta && eta < 0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.955878) + 
    (-0.00 <= eta && eta < 0.10) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.931602) + 
    (-0.00 <= eta && eta < 0.10) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.908797) + 
    (-0.00 <= eta && eta < 0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.885676) + 
    (-0.00 <= eta && eta < 0.10) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.862280) + 
    (-0.00 <= eta && eta < 0.10) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.839159) + 
    (-0.00 <= eta && eta < 0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.815846) + 
    (-0.00 <= eta && eta < 0.10) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.795129) + 
    (-0.00 <= eta && eta < 0.10) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.773860) + 
    (-0.00 <= eta && eta < 0.10) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.751980) + 
    (-0.00 <= eta && eta < 0.10) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.728013) + 
    (-0.00 <= eta && eta < 0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.705518) + 
    (-0.00 <= eta && eta < 0.10) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.683919) + 
    (-0.00 <= eta && eta < 0.10) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.661467) + 
    (-0.00 <= eta && eta < 0.10) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.638022) + 
    (-0.00 <= eta && eta < 0.10) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604731) + 
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.10 <= eta && eta < 0.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999390) + 
    (0.10 <= eta && eta < 0.20) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978436) + 
    (0.10 <= eta && eta < 0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.954026) + 
    (0.10 <= eta && eta < 0.20) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.928980) + 
    (0.10 <= eta && eta < 0.20) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.905261) + 
    (0.10 <= eta && eta < 0.20) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.880827) + 
    (0.10 <= eta && eta < 0.20) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.855798) + 
    (0.10 <= eta && eta < 0.20) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.831816) + 
    (0.10 <= eta && eta < 0.20) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.811017) + 
    (0.10 <= eta && eta < 0.20) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.790324) + 
    (0.10 <= eta && eta < 0.20) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.766845) + 
    (0.10 <= eta && eta < 0.20) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.742732) + 
    (0.10 <= eta && eta < 0.20) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.721120) + 
    (0.10 <= eta && eta < 0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.698097) + 
    (0.10 <= eta && eta < 0.20) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.675216) + 
    (0.10 <= eta && eta < 0.20) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.652610) + 
    (0.10 <= eta && eta < 0.20) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.627707) + 
    (0.10 <= eta && eta < 0.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603627) + 
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.20 <= eta && eta < 0.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999486) + 
    (0.20 <= eta && eta < 0.30) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.979519) + 
    (0.20 <= eta && eta < 0.30) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.955406) + 
    (0.20 <= eta && eta < 0.30) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.930947) + 
    (0.20 <= eta && eta < 0.30) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.906986) + 
    (0.20 <= eta && eta < 0.30) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.882189) + 
    (0.20 <= eta && eta < 0.30) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.857189) + 
    (0.20 <= eta && eta < 0.30) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.833268) + 
    (0.20 <= eta && eta < 0.30) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.811348) + 
    (0.20 <= eta && eta < 0.30) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.787352) + 
    (0.20 <= eta && eta < 0.30) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.762475) + 
    (0.20 <= eta && eta < 0.30) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.737585) + 
    (0.20 <= eta && eta < 0.30) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.716144) + 
    (0.20 <= eta && eta < 0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.693021) + 
    (0.20 <= eta && eta < 0.30) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.670061) + 
    (0.20 <= eta && eta < 0.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.646891) + 
    (0.20 <= eta && eta < 0.30) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605262) + 
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.30 <= eta && eta < 0.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999459) + 
    (0.30 <= eta && eta < 0.40) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.979305) + 
    (0.30 <= eta && eta < 0.40) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.955111) + 
    (0.30 <= eta && eta < 0.40) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.929350) + 
    (0.30 <= eta && eta < 0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.904229) + 
    (0.30 <= eta && eta < 0.40) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.880114) + 
    (0.30 <= eta && eta < 0.40) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.855471) + 
    (0.30 <= eta && eta < 0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.831596) + 
    (0.30 <= eta && eta < 0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.808837) + 
    (0.30 <= eta && eta < 0.40) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.785228) + 
    (0.30 <= eta && eta < 0.40) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.756982) + 
    (0.30 <= eta && eta < 0.40) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.732645) + 
    (0.30 <= eta && eta < 0.40) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.709880) + 
    (0.30 <= eta && eta < 0.40) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.687128) + 
    (0.30 <= eta && eta < 0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.663925) + 
    (0.30 <= eta && eta < 0.40) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.640354) + 
    (0.30 <= eta && eta < 0.40) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604740) + 
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.40 <= eta && eta < 0.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999416) + 
    (0.40 <= eta && eta < 0.50) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.977866) + 
    (0.40 <= eta && eta < 0.50) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.956023) + 
    (0.40 <= eta && eta < 0.50) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.934804) + 
    (0.40 <= eta && eta < 0.50) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.910945) + 
    (0.40 <= eta && eta < 0.50) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.886157) + 
    (0.40 <= eta && eta < 0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.861028) + 
    (0.40 <= eta && eta < 0.50) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.836587) + 
    (0.40 <= eta && eta < 0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.813125) + 
    (0.40 <= eta && eta < 0.50) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.788762) + 
    (0.40 <= eta && eta < 0.50) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.764916) + 
    (0.40 <= eta && eta < 0.50) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.739594) + 
    (0.40 <= eta && eta < 0.50) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.715349) + 
    (0.40 <= eta && eta < 0.50) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.691540) + 
    (0.40 <= eta && eta < 0.50) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.668971) + 
    (0.40 <= eta && eta < 0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.645517) + 
    (0.40 <= eta && eta < 0.50) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605109) + 
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.50 <= eta && eta < 0.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999374) + 
    (0.50 <= eta && eta < 0.60) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.976949) + 
    (0.50 <= eta && eta < 0.60) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.951356) + 
    (0.50 <= eta && eta < 0.60) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.926120) + 
    (0.50 <= eta && eta < 0.60) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.902566) + 
    (0.50 <= eta && eta < 0.60) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.878135) + 
    (0.50 <= eta && eta < 0.60) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.853567) + 
    (0.50 <= eta && eta < 0.60) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.830007) + 
    (0.50 <= eta && eta < 0.60) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.807715) + 
    (0.50 <= eta && eta < 0.60) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.784361) + 
    (0.50 <= eta && eta < 0.60) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.760491) + 
    (0.50 <= eta && eta < 0.60) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.736736) + 
    (0.50 <= eta && eta < 0.60) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.712414) + 
    (0.50 <= eta && eta < 0.60) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.690695) + 
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.668691) + 
    (0.50 <= eta && eta < 0.60) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.645647) + 
    (0.50 <= eta && eta < 0.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605221) + 
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.60 <= eta && eta < 0.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999428) + 
    (0.60 <= eta && eta < 0.70) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978093) + 
    (0.60 <= eta && eta < 0.70) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.952896) + 
    (0.60 <= eta && eta < 0.70) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.928743) + 
    (0.60 <= eta && eta < 0.70) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.904835) + 
    (0.60 <= eta && eta < 0.70) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.880239) + 
    (0.60 <= eta && eta < 0.70) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.855800) + 
    (0.60 <= eta && eta < 0.70) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.832129) + 
    (0.60 <= eta && eta < 0.70) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.809971) + 
    (0.60 <= eta && eta < 0.70) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.785533) + 
    (0.60 <= eta && eta < 0.70) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.762384) + 
    (0.60 <= eta && eta < 0.70) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.739140) + 
    (0.60 <= eta && eta < 0.70) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715270) + 
    (0.60 <= eta && eta < 0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.692367) + 
    (0.60 <= eta && eta < 0.70) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.669998) + 
    (0.60 <= eta && eta < 0.70) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647168) + 
    (0.60 <= eta && eta < 0.70) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605392) + 
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.70 <= eta && eta < 0.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999431) + 
    (0.70 <= eta && eta < 0.80) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978854) + 
    (0.70 <= eta && eta < 0.80) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.957614) + 
    (0.70 <= eta && eta < 0.80) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.936746) + 
    (0.70 <= eta && eta < 0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.914572) + 
    (0.70 <= eta && eta < 0.80) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.889412) + 
    (0.70 <= eta && eta < 0.80) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.863632) + 
    (0.70 <= eta && eta < 0.80) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.839744) + 
    (0.70 <= eta && eta < 0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.814966) + 
    (0.70 <= eta && eta < 0.80) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.787803) + 
    (0.70 <= eta && eta < 0.80) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.762097) + 
    (0.70 <= eta && eta < 0.80) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.738550) + 
    (0.70 <= eta && eta < 0.80) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715206) + 
    (0.70 <= eta && eta < 0.80) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.692151) + 
    (0.70 <= eta && eta < 0.80) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.669829) + 
    (0.70 <= eta && eta < 0.80) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647051) + 
    (0.70 <= eta && eta < 0.80) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605381) + 
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.80 <= eta && eta < 0.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999463) + 
    (0.80 <= eta && eta < 0.90) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.979480) + 
    (0.80 <= eta && eta < 0.90) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.955658) + 
    (0.80 <= eta && eta < 0.90) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.931795) + 
    (0.80 <= eta && eta < 0.90) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.909021) + 
    (0.80 <= eta && eta < 0.90) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.884511) + 
    (0.80 <= eta && eta < 0.90) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.859876) + 
    (0.80 <= eta && eta < 0.90) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.837206) + 
    (0.80 <= eta && eta < 0.90) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.815556) + 
    (0.80 <= eta && eta < 0.90) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.792098) + 
    (0.80 <= eta && eta < 0.90) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.767572) + 
    (0.80 <= eta && eta < 0.90) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.742973) + 
    (0.80 <= eta && eta < 0.90) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.719838) + 
    (0.80 <= eta && eta < 0.90) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.696957) + 
    (0.80 <= eta && eta < 0.90) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.674906) + 
    (0.80 <= eta && eta < 0.90) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.652868) + 
    (0.80 <= eta && eta < 0.90) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.628233) + 
    (0.80 <= eta && eta < 0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603774) + 
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.90 <= eta && eta < 1.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999490) + 
    (0.90 <= eta && eta < 1.00) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.980505) + 
    (0.90 <= eta && eta < 1.00) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.957300) + 
    (0.90 <= eta && eta < 1.00) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.933891) + 
    (0.90 <= eta && eta < 1.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.910870) + 
    (0.90 <= eta && eta < 1.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.887137) + 
    (0.90 <= eta && eta < 1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.863393) + 
    (0.90 <= eta && eta < 1.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.840634) + 
    (0.90 <= eta && eta < 1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.818273) + 
    (0.90 <= eta && eta < 1.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.794422) + 
    (0.90 <= eta && eta < 1.00) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.769639) + 
    (0.90 <= eta && eta < 1.00) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.748157) + 
    (0.90 <= eta && eta < 1.00) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.724553) + 
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.699984) + 
    (0.90 <= eta && eta < 1.00) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.677138) + 
    (0.90 <= eta && eta < 1.00) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.654326) + 
    (0.90 <= eta && eta < 1.00) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.629520) + 
    (0.90 <= eta && eta < 1.00) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603898) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 

  add EfficiencyFormula {321} {-211} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.000629) + 
    (-1.00 <= eta && eta < -0.90) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.020818) + 
    (-1.00 <= eta && eta < -0.90) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.043238) + 
    (-1.00 <= eta && eta < -0.90) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.064987) + 
    (-1.00 <= eta && eta < -0.90) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.086643) + 
    (-1.00 <= eta && eta < -0.90) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.108246) + 
    (-1.00 <= eta && eta < -0.90) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.131038) + 
    (-1.00 <= eta && eta < -0.90) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.153265) + 
    (-1.00 <= eta && eta < -0.90) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.173948) + 
    (-1.00 <= eta && eta < -0.90) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.196026) + 
    (-1.00 <= eta && eta < -0.90) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.219035) + 
    (-1.00 <= eta && eta < -0.90) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.242086) + 
    (-1.00 <= eta && eta < -0.90) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.264737) + 
    (-1.00 <= eta && eta < -0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.286490) + 
    (-1.00 <= eta && eta < -0.90) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.308418) + 
    (-1.00 <= eta && eta < -0.90) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.329730) + 
    (-1.00 <= eta && eta < -0.90) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.350771) + 
    (-1.00 <= eta && eta < -0.90) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.372200) + 
    (-1.00 <= eta && eta < -0.90) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.393662) + 
    (-1.00 <= eta && eta < -0.90) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.415409) + 
    (-1.00 <= eta && eta < -0.90) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.437641) + 
    (-1.00 <= eta && eta < -0.90) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.461287) + 
    (-1.00 <= eta && eta < -0.90) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.480405) + 
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000601) + 
    (-0.90 <= eta && eta < -0.80) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.020257) + 
    (-0.90 <= eta && eta < -0.80) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.043042) + 
    (-0.90 <= eta && eta < -0.80) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.065328) + 
    (-0.90 <= eta && eta < -0.80) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.087115) + 
    (-0.90 <= eta && eta < -0.80) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.109928) + 
    (-0.90 <= eta && eta < -0.80) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.132996) + 
    (-0.90 <= eta && eta < -0.80) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.156649) + 
    (-0.90 <= eta && eta < -0.80) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.180600) + 
    (-0.90 <= eta && eta < -0.80) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.203148) + 
    (-0.90 <= eta && eta < -0.80) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.225705) + 
    (-0.90 <= eta && eta < -0.80) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.248439) + 
    (-0.90 <= eta && eta < -0.80) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.270714) + 
    (-0.90 <= eta && eta < -0.80) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.292054) + 
    (-0.90 <= eta && eta < -0.80) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.313521) + 
    (-0.90 <= eta && eta < -0.80) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.335443) + 
    (-0.90 <= eta && eta < -0.80) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.357624) + 
    (-0.90 <= eta && eta < -0.80) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.379405) + 
    (-0.90 <= eta && eta < -0.80) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.401051) + 
    (-0.90 <= eta && eta < -0.80) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.422825) + 
    (-0.90 <= eta && eta < -0.80) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.445280) + 
    (-0.90 <= eta && eta < -0.80) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 41.00) * (0.470610) + 
    (-0.90 <= eta && eta < -0.80) * (41.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484124) + 
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000757) + 
    (-0.80 <= eta && eta < -0.70) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.022320) + 
    (-0.80 <= eta && eta < -0.70) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.044193) + 
    (-0.80 <= eta && eta < -0.70) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.067162) + 
    (-0.80 <= eta && eta < -0.70) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.088655) + 
    (-0.80 <= eta && eta < -0.70) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.111401) + 
    (-0.80 <= eta && eta < -0.70) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.134165) + 
    (-0.80 <= eta && eta < -0.70) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.156920) + 
    (-0.80 <= eta && eta < -0.70) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.177718) + 
    (-0.80 <= eta && eta < -0.70) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.199864) + 
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.223340) + 
    (-0.80 <= eta && eta < -0.70) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.246613) + 
    (-0.80 <= eta && eta < -0.70) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.269390) + 
    (-0.80 <= eta && eta < -0.70) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.291176) + 
    (-0.80 <= eta && eta < -0.70) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.313049) + 
    (-0.80 <= eta && eta < -0.70) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.335333) + 
    (-0.80 <= eta && eta < -0.70) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.357820) + 
    (-0.80 <= eta && eta < -0.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.379834) + 
    (-0.80 <= eta && eta < -0.70) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.401637) + 
    (-0.80 <= eta && eta < -0.70) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.423486) + 
    (-0.80 <= eta && eta < -0.70) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.445919) + 
    (-0.80 <= eta && eta < -0.70) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 41.30) * (0.471428) + 
    (-0.80 <= eta && eta < -0.70) * (41.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484680) + 
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000627) + 
    (-0.70 <= eta && eta < -0.60) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.019743) + 
    (-0.70 <= eta && eta < -0.60) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.041036) + 
    (-0.70 <= eta && eta < -0.60) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.063090) + 
    (-0.70 <= eta && eta < -0.60) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.084632) + 
    (-0.70 <= eta && eta < -0.60) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.107944) + 
    (-0.70 <= eta && eta < -0.60) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.131352) + 
    (-0.70 <= eta && eta < -0.60) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.153268) + 
    (-0.70 <= eta && eta < -0.60) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.174594) + 
    (-0.70 <= eta && eta < -0.60) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.197491) + 
    (-0.70 <= eta && eta < -0.60) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.221492) + 
    (-0.70 <= eta && eta < -0.60) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.243056) + 
    (-0.70 <= eta && eta < -0.60) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.264487) + 
    (-0.70 <= eta && eta < -0.60) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.286988) + 
    (-0.70 <= eta && eta < -0.60) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.309562) + 
    (-0.70 <= eta && eta < -0.60) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.331387) + 
    (-0.70 <= eta && eta < -0.60) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.352810) + 
    (-0.70 <= eta && eta < -0.60) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374492) + 
    (-0.70 <= eta && eta < -0.60) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.396058) + 
    (-0.70 <= eta && eta < -0.60) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417745) + 
    (-0.70 <= eta && eta < -0.60) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.439965) + 
    (-0.70 <= eta && eta < -0.60) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.463862) + 
    (-0.70 <= eta && eta < -0.60) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481927) + 
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000606) + 
    (-0.60 <= eta && eta < -0.50) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.019510) + 
    (-0.60 <= eta && eta < -0.50) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.040196) + 
    (-0.60 <= eta && eta < -0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.061442) + 
    (-0.60 <= eta && eta < -0.50) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.083297) + 
    (-0.60 <= eta && eta < -0.50) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.105846) + 
    (-0.60 <= eta && eta < -0.50) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.129010) + 
    (-0.60 <= eta && eta < -0.50) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.151227) + 
    (-0.60 <= eta && eta < -0.50) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.173018) + 
    (-0.60 <= eta && eta < -0.50) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.197095) + 
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.220621) + 
    (-0.60 <= eta && eta < -0.50) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.242223) + 
    (-0.60 <= eta && eta < -0.50) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.263700) + 
    (-0.60 <= eta && eta < -0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.286256) + 
    (-0.60 <= eta && eta < -0.50) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.308892) + 
    (-0.60 <= eta && eta < -0.50) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.330782) + 
    (-0.60 <= eta && eta < -0.50) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.352274) + 
    (-0.60 <= eta && eta < -0.50) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374028) + 
    (-0.60 <= eta && eta < -0.50) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.395670) + 
    (-0.60 <= eta && eta < -0.50) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417435) + 
    (-0.60 <= eta && eta < -0.50) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.439737) + 
    (-0.60 <= eta && eta < -0.50) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.463725) + 
    (-0.60 <= eta && eta < -0.50) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481858) + 
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000754) + 
    (-0.50 <= eta && eta < -0.40) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.022057) + 
    (-0.50 <= eta && eta < -0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.044626) + 
    (-0.50 <= eta && eta < -0.40) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.067546) + 
    (-0.50 <= eta && eta < -0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.089482) + 
    (-0.50 <= eta && eta < -0.40) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.112586) + 
    (-0.50 <= eta && eta < -0.40) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.135083) + 
    (-0.50 <= eta && eta < -0.40) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.158756) + 
    (-0.50 <= eta && eta < -0.40) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.180504) + 
    (-0.50 <= eta && eta < -0.40) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.204102) + 
    (-0.50 <= eta && eta < -0.40) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.227320) + 
    (-0.50 <= eta && eta < -0.40) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.250611) + 
    (-0.50 <= eta && eta < -0.40) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.273341) + 
    (-0.50 <= eta && eta < -0.40) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.295022) + 
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.316735) + 
    (-0.50 <= eta && eta < -0.40) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.338798) + 
    (-0.50 <= eta && eta < -0.40) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.361005) + 
    (-0.50 <= eta && eta < -0.40) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.382690) + 
    (-0.50 <= eta && eta < -0.40) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.404114) + 
    (-0.50 <= eta && eta < -0.40) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425866) + 
    (-0.50 <= eta && eta < -0.40) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448463) + 
    (-0.50 <= eta && eta < -0.40) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.475401) + 
    (-0.50 <= eta && eta < -0.40) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486661) + 
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000728) + 
    (-0.40 <= eta && eta < -0.30) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.021373) + 
    (-0.40 <= eta && eta < -0.30) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.042787) + 
    (-0.40 <= eta && eta < -0.30) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.065124) + 
    (-0.40 <= eta && eta < -0.30) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.087923) + 
    (-0.40 <= eta && eta < -0.30) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.111671) + 
    (-0.40 <= eta && eta < -0.30) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.134974) + 
    (-0.40 <= eta && eta < -0.30) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.157718) + 
    (-0.40 <= eta && eta < -0.30) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.179824) + 
    (-0.40 <= eta && eta < -0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.202629) + 
    (-0.40 <= eta && eta < -0.30) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.225577) + 
    (-0.40 <= eta && eta < -0.30) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.248957) + 
    (-0.40 <= eta && eta < -0.30) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.271791) + 
    (-0.40 <= eta && eta < -0.30) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.293585) + 
    (-0.40 <= eta && eta < -0.30) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.315422) + 
    (-0.40 <= eta && eta < -0.30) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.337621) + 
    (-0.40 <= eta && eta < -0.30) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.359975) + 
    (-0.40 <= eta && eta < -0.30) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.381809) + 
    (-0.40 <= eta && eta < -0.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.403386) + 
    (-0.40 <= eta && eta < -0.30) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425299) + 
    (-0.40 <= eta && eta < -0.30) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448066) + 
    (-0.40 <= eta && eta < -0.30) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.474824) + 
    (-0.40 <= eta && eta < -0.30) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486349) + 
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000713) + 
    (-0.30 <= eta && eta < -0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.021188) + 
    (-0.30 <= eta && eta < -0.20) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.043291) + 
    (-0.30 <= eta && eta < -0.20) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.065636) + 
    (-0.30 <= eta && eta < -0.20) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.087346) + 
    (-0.30 <= eta && eta < -0.20) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.111189) + 
    (-0.30 <= eta && eta < -0.20) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.135124) + 
    (-0.30 <= eta && eta < -0.20) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.157289) + 
    (-0.30 <= eta && eta < -0.20) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.178343) + 
    (-0.30 <= eta && eta < -0.20) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.201490) + 
    (-0.30 <= eta && eta < -0.20) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.225445) + 
    (-0.30 <= eta && eta < -0.20) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.248831) + 
    (-0.30 <= eta && eta < -0.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.271673) + 
    (-0.30 <= eta && eta < -0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.293475) + 
    (-0.30 <= eta && eta < -0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.315321) + 
    (-0.30 <= eta && eta < -0.20) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.337531) + 
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.359896) + 
    (-0.30 <= eta && eta < -0.20) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.381742) + 
    (-0.30 <= eta && eta < -0.20) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.403330) + 
    (-0.30 <= eta && eta < -0.20) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425256) + 
    (-0.30 <= eta && eta < -0.20) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448036) + 
    (-0.30 <= eta && eta < -0.20) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.474810) + 
    (-0.30 <= eta && eta < -0.20) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486341) + 
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000670) + 
    (-0.20 <= eta && eta < -0.10) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.020397) + 
    (-0.20 <= eta && eta < -0.10) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.041717) + 
    (-0.20 <= eta && eta < -0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.064356) + 
    (-0.20 <= eta && eta < -0.10) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.085635) + 
    (-0.20 <= eta && eta < -0.10) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.109110) + 
    (-0.20 <= eta && eta < -0.10) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.131602) + 
    (-0.20 <= eta && eta < -0.10) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.153486) + 
    (-0.20 <= eta && eta < -0.10) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.175586) + 
    (-0.20 <= eta && eta < -0.10) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.198903) + 
    (-0.20 <= eta && eta < -0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.222061) + 
    (-0.20 <= eta && eta < -0.10) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.243601) + 
    (-0.20 <= eta && eta < -0.10) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.265002) + 
    (-0.20 <= eta && eta < -0.10) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.287466) + 
    (-0.20 <= eta && eta < -0.10) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.310000) + 
    (-0.20 <= eta && eta < -0.10) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.331781) + 
    (-0.20 <= eta && eta < -0.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.353161) + 
    (-0.20 <= eta && eta < -0.10) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374795) + 
    (-0.20 <= eta && eta < -0.10) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.396312) + 
    (-0.20 <= eta && eta < -0.10) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417947) + 
    (-0.20 <= eta && eta < -0.10) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.440351) + 
    (-0.20 <= eta && eta < -0.10) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 35.00) * (0.464506) + 
    (-0.20 <= eta && eta < -0.10) * (35.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482178) + 
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.000742) + 
    (-0.10 <= eta && eta < -0.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.021676) + 
    (-0.10 <= eta && eta < -0.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.043000) + 
    (-0.10 <= eta && eta < -0.00) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.064909) + 
    (-0.10 <= eta && eta < -0.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.086553) + 
    (-0.10 <= eta && eta < -0.00) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.108414) + 
    (-0.10 <= eta && eta < -0.00) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.131561) + 
    (-0.10 <= eta && eta < -0.00) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.155629) + 
    (-0.10 <= eta && eta < -0.00) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.177181) + 
    (-0.10 <= eta && eta < -0.00) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.200013) + 
    (-0.10 <= eta && eta < -0.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.223039) + 
    (-0.10 <= eta && eta < -0.00) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.245902) + 
    (-0.10 <= eta && eta < -0.00) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.268329) + 
    (-0.10 <= eta && eta < -0.00) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.289835) + 
    (-0.10 <= eta && eta < -0.00) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.311487) + 
    (-0.10 <= eta && eta < -0.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.333613) + 
    (-0.10 <= eta && eta < -0.00) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.356016) + 
    (-0.10 <= eta && eta < -0.00) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.378024) + 
    (-0.10 <= eta && eta < -0.00) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.399906) + 
    (-0.10 <= eta && eta < -0.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.421924) + 
    (-0.10 <= eta && eta < -0.00) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.90) * (0.444432) + 
    (-0.10 <= eta && eta < -0.00) * (26.90 <= pt * cosh(eta) && pt * cosh(eta) < 40.00) * (0.469406) + 
    (-0.10 <= eta && eta < -0.00) * (40.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.483535) + 
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000700) + 
    (-0.00 <= eta && eta < 0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.021257) + 
    (-0.00 <= eta && eta < 0.10) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.043040) + 
    (-0.00 <= eta && eta < 0.10) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.064435) + 
    (-0.00 <= eta && eta < 0.10) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.086042) + 
    (-0.00 <= eta && eta < 0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.109607) + 
    (-0.00 <= eta && eta < 0.10) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.132930) + 
    (-0.00 <= eta && eta < 0.10) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.153931) + 
    (-0.00 <= eta && eta < 0.10) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.175705) + 
    (-0.00 <= eta && eta < 0.10) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.197877) + 
    (-0.00 <= eta && eta < 0.10) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220292) + 
    (-0.00 <= eta && eta < 0.10) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.241706) + 
    (-0.00 <= eta && eta < 0.10) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263021) + 
    (-0.00 <= eta && eta < 0.10) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285436) + 
    (-0.00 <= eta && eta < 0.10) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.307965) + 
    (-0.00 <= eta && eta < 0.10) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.329785) + 
    (-0.00 <= eta && eta < 0.10) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351244) + 
    (-0.00 <= eta && eta < 0.10) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373004) + 
    (-0.00 <= eta && eta < 0.10) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394692) + 
    (-0.00 <= eta && eta < 0.10) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416547) + 
    (-0.00 <= eta && eta < 0.10) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.438994) + 
    (-0.00 <= eta && eta < 0.10) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.462884) + 
    (-0.00 <= eta && eta < 0.10) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481363) + 
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000617) + 
    (0.10 <= eta && eta < 0.20) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.019051) + 
    (0.10 <= eta && eta < 0.20) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.040330) + 
    (0.10 <= eta && eta < 0.20) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.062862) + 
    (0.10 <= eta && eta < 0.20) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085532) + 
    (0.10 <= eta && eta < 0.20) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.107916) + 
    (0.10 <= eta && eta < 0.20) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.131693) + 
    (0.10 <= eta && eta < 0.20) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.155745) + 
    (0.10 <= eta && eta < 0.20) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.177142) + 
    (0.10 <= eta && eta < 0.20) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.196263) + 
    (0.10 <= eta && eta < 0.20) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.218643) + 
    (0.10 <= eta && eta < 0.20) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.240921) + 
    (0.10 <= eta && eta < 0.20) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.262860) + 
    (0.10 <= eta && eta < 0.20) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.285861) + 
    (0.10 <= eta && eta < 0.20) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.308892) + 
    (0.10 <= eta && eta < 0.20) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.331107) + 
    (0.10 <= eta && eta < 0.20) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.352856) + 
    (0.10 <= eta && eta < 0.20) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.374803) + 
    (0.10 <= eta && eta < 0.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.396560) + 
    (0.10 <= eta && eta < 0.20) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.418357) + 
    (0.10 <= eta && eta < 0.20) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.440831) + 
    (0.10 <= eta && eta < 0.20) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.465032) + 
    (0.10 <= eta && eta < 0.20) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482582) + 
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000736) + 
    (0.20 <= eta && eta < 0.30) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.021548) + 
    (0.20 <= eta && eta < 0.30) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.043947) + 
    (0.20 <= eta && eta < 0.30) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.067742) + 
    (0.20 <= eta && eta < 0.30) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.090454) + 
    (0.20 <= eta && eta < 0.30) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.113386) + 
    (0.20 <= eta && eta < 0.30) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.136843) + 
    (0.20 <= eta && eta < 0.30) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.160311) + 
    (0.20 <= eta && eta < 0.30) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.185362) + 
    (0.20 <= eta && eta < 0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.207279) + 
    (0.20 <= eta && eta < 0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.230256) + 
    (0.20 <= eta && eta < 0.30) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.253794) + 
    (0.20 <= eta && eta < 0.30) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.276522) + 
    (0.20 <= eta && eta < 0.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.298150) + 
    (0.20 <= eta && eta < 0.30) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.319758) + 
    (0.20 <= eta && eta < 0.30) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.341663) + 
    (0.20 <= eta && eta < 0.30) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.363659) + 
    (0.20 <= eta && eta < 0.30) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.385086) + 
    (0.20 <= eta && eta < 0.30) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.406698) + 
    (0.20 <= eta && eta < 0.30) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.428629) + 
    (0.20 <= eta && eta < 0.30) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.60) * (0.451291) + 
    (0.20 <= eta && eta < 0.30) * (27.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.478658) + 
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000753) + 
    (0.30 <= eta && eta < 0.40) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.021899) + 
    (0.30 <= eta && eta < 0.40) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.045558) + 
    (0.30 <= eta && eta < 0.40) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.068436) + 
    (0.30 <= eta && eta < 0.40) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.091463) + 
    (0.30 <= eta && eta < 0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.115435) + 
    (0.30 <= eta && eta < 0.40) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.139105) + 
    (0.30 <= eta && eta < 0.40) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.162598) + 
    (0.30 <= eta && eta < 0.40) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.182795) + 
    (0.30 <= eta && eta < 0.40) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.205458) + 
    (0.30 <= eta && eta < 0.40) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.228047) + 
    (0.30 <= eta && eta < 0.40) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.249350) + 
    (0.30 <= eta && eta < 0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.270616) + 
    (0.30 <= eta && eta < 0.40) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.292865) + 
    (0.30 <= eta && eta < 0.40) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.315110) + 
    (0.30 <= eta && eta < 0.40) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.336546) + 
    (0.30 <= eta && eta < 0.40) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.357525) + 
    (0.30 <= eta && eta < 0.40) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.378695) + 
    (0.30 <= eta && eta < 0.40) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.400236) + 
    (0.30 <= eta && eta < 0.40) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.422294) + 
    (0.30 <= eta && eta < 0.40) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.444843) + 
    (0.30 <= eta && eta < 0.40) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.469965) + 
    (0.30 <= eta && eta < 0.40) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484558) + 
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.70) * (0.000629) + 
    (0.40 <= eta && eta < 0.50) * (5.70 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.019727) + 
    (0.40 <= eta && eta < 0.50) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.041772) + 
    (0.40 <= eta && eta < 0.50) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.064756) + 
    (0.40 <= eta && eta < 0.50) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.087848) + 
    (0.40 <= eta && eta < 0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.112141) + 
    (0.40 <= eta && eta < 0.50) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.135408) + 
    (0.40 <= eta && eta < 0.50) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.157790) + 
    (0.40 <= eta && eta < 0.50) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.180619) + 
    (0.40 <= eta && eta < 0.50) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.203357) + 
    (0.40 <= eta && eta < 0.50) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.227838) + 
    (0.40 <= eta && eta < 0.50) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.249584) + 
    (0.40 <= eta && eta < 0.50) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.271029) + 
    (0.40 <= eta && eta < 0.50) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.293437) + 
    (0.40 <= eta && eta < 0.50) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.315808) + 
    (0.40 <= eta && eta < 0.50) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.337334) + 
    (0.40 <= eta && eta < 0.50) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.358365) + 
    (0.40 <= eta && eta < 0.50) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.379551) + 
    (0.40 <= eta && eta < 0.50) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.401066) + 
    (0.40 <= eta && eta < 0.50) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.423052) + 
    (0.40 <= eta && eta < 0.50) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.445684) + 
    (0.40 <= eta && eta < 0.50) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.471163) + 
    (0.40 <= eta && eta < 0.50) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.485174) + 
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000656) + 
    (0.50 <= eta && eta < 0.60) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.020312) + 
    (0.50 <= eta && eta < 0.60) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.041999) + 
    (0.50 <= eta && eta < 0.60) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.064530) + 
    (0.50 <= eta && eta < 0.60) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.087046) + 
    (0.50 <= eta && eta < 0.60) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.110549) + 
    (0.50 <= eta && eta < 0.60) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.133710) + 
    (0.50 <= eta && eta < 0.60) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.156249) + 
    (0.50 <= eta && eta < 0.60) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.178337) + 
    (0.50 <= eta && eta < 0.60) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.199630) + 
    (0.50 <= eta && eta < 0.60) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.220540) + 
    (0.50 <= eta && eta < 0.60) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.242222) + 
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.264089) + 
    (0.50 <= eta && eta < 0.60) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.287002) + 
    (0.50 <= eta && eta < 0.60) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.309935) + 
    (0.50 <= eta && eta < 0.60) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.332047) + 
    (0.50 <= eta && eta < 0.60) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.353688) + 
    (0.50 <= eta && eta < 0.60) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.375520) + 
    (0.50 <= eta && eta < 0.60) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.397160) + 
    (0.50 <= eta && eta < 0.60) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.418834) + 
    (0.50 <= eta && eta < 0.60) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.441179) + 
    (0.50 <= eta && eta < 0.60) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.465437) + 
    (0.50 <= eta && eta < 0.60) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482784) + 
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000610) + 
    (0.60 <= eta && eta < 0.70) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.019574) + 
    (0.60 <= eta && eta < 0.70) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.040602) + 
    (0.60 <= eta && eta < 0.70) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.063109) + 
    (0.60 <= eta && eta < 0.70) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085463) + 
    (0.60 <= eta && eta < 0.70) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.108771) + 
    (0.60 <= eta && eta < 0.70) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.132650) + 
    (0.60 <= eta && eta < 0.70) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.154383) + 
    (0.60 <= eta && eta < 0.70) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.175861) + 
    (0.60 <= eta && eta < 0.70) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.198600) + 
    (0.60 <= eta && eta < 0.70) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.222323) + 
    (0.60 <= eta && eta < 0.70) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.244096) + 
    (0.60 <= eta && eta < 0.70) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.265662) + 
    (0.60 <= eta && eta < 0.70) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.288269) + 
    (0.60 <= eta && eta < 0.70) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.310912) + 
    (0.60 <= eta && eta < 0.70) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.332765) + 
    (0.60 <= eta && eta < 0.70) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.354177) + 
    (0.60 <= eta && eta < 0.70) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.375807) + 
    (0.60 <= eta && eta < 0.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.397278) + 
    (0.60 <= eta && eta < 0.70) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.418823) + 
    (0.60 <= eta && eta < 0.70) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.441082) + 
    (0.60 <= eta && eta < 0.70) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.465310) + 
    (0.60 <= eta && eta < 0.70) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482651) + 
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000768) + 
    (0.70 <= eta && eta < 0.80) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.022252) + 
    (0.70 <= eta && eta < 0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.042122) + 
    (0.70 <= eta && eta < 0.80) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.062694) + 
    (0.70 <= eta && eta < 0.80) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085610) + 
    (0.70 <= eta && eta < 0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.107607) + 
    (0.70 <= eta && eta < 0.80) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.130626) + 
    (0.70 <= eta && eta < 0.80) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.154667) + 
    (0.70 <= eta && eta < 0.80) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.176433) + 
    (0.70 <= eta && eta < 0.80) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.198617) + 
    (0.70 <= eta && eta < 0.80) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.222620) + 
    (0.70 <= eta && eta < 0.80) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.244391) + 
    (0.70 <= eta && eta < 0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.265940) + 
    (0.70 <= eta && eta < 0.80) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.288527) + 
    (0.70 <= eta && eta < 0.80) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.311148) + 
    (0.70 <= eta && eta < 0.80) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.332977) + 
    (0.70 <= eta && eta < 0.80) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.354366) + 
    (0.70 <= eta && eta < 0.80) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.375969) + 
    (0.70 <= eta && eta < 0.80) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.397969) + 
    (0.70 <= eta && eta < 0.80) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.420136) + 
    (0.70 <= eta && eta < 0.80) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.442635) + 
    (0.70 <= eta && eta < 0.80) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.60) * (0.467206) + 
    (0.70 <= eta && eta < 0.80) * (36.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.483337) + 
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000728) + 
    (0.80 <= eta && eta < 0.90) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.021330) + 
    (0.80 <= eta && eta < 0.90) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.042978) + 
    (0.80 <= eta && eta < 0.90) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.065797) + 
    (0.80 <= eta && eta < 0.90) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.087589) + 
    (0.80 <= eta && eta < 0.90) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.109904) + 
    (0.80 <= eta && eta < 0.90) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.132432) + 
    (0.80 <= eta && eta < 0.90) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.155030) + 
    (0.80 <= eta && eta < 0.90) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.177015) + 
    (0.80 <= eta && eta < 0.90) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.198109) + 
    (0.80 <= eta && eta < 0.90) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220572) + 
    (0.80 <= eta && eta < 0.90) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.241975) + 
    (0.80 <= eta && eta < 0.90) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263275) + 
    (0.80 <= eta && eta < 0.90) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285672) + 
    (0.80 <= eta && eta < 0.90) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.308182) + 
    (0.80 <= eta && eta < 0.90) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.329981) + 
    (0.80 <= eta && eta < 0.90) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351418) + 
    (0.80 <= eta && eta < 0.90) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373154) + 
    (0.80 <= eta && eta < 0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394818) + 
    (0.80 <= eta && eta < 0.90) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416648) + 
    (0.80 <= eta && eta < 0.90) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.439069) + 
    (0.80 <= eta && eta < 0.90) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.463037) + 
    (0.80 <= eta && eta < 0.90) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481439) + 
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000687) + 
    (0.90 <= eta && eta < 1.00) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.020501) + 
    (0.90 <= eta && eta < 1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.041727) + 
    (0.90 <= eta && eta < 1.00) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.063800) + 
    (0.90 <= eta && eta < 1.00) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.084994) + 
    (0.90 <= eta && eta < 1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.107765) + 
    (0.90 <= eta && eta < 1.00) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.130564) + 
    (0.90 <= eta && eta < 1.00) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.152720) + 
    (0.90 <= eta && eta < 1.00) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.174241) + 
    (0.90 <= eta && eta < 1.00) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.197247) + 
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220892) + 
    (0.90 <= eta && eta < 1.00) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.242281) + 
    (0.90 <= eta && eta < 1.00) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263564) + 
    (0.90 <= eta && eta < 1.00) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285942) + 
    (0.90 <= eta && eta < 1.00) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.308428) + 
    (0.90 <= eta && eta < 1.00) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.330204) + 
    (0.90 <= eta && eta < 1.00) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351616) + 
    (0.90 <= eta && eta < 1.00) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373325) + 
    (0.90 <= eta && eta < 1.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394961) + 
    (0.90 <= eta && eta < 1.00) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416763) + 
    (0.90 <= eta && eta < 1.00) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.439153) + 
    (0.90 <= eta && eta < 1.00) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.463089) + 
    (0.90 <= eta && eta < 1.00) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481465) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 

  add EfficiencyFormula {321} {2212} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.000705) + 
    (-1.00 <= eta && eta < -0.90) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.021025) + 
    (-1.00 <= eta && eta < -0.90) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.042481) + 
    (-1.00 <= eta && eta < -0.90) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.063884) + 
    (-1.00 <= eta && eta < -0.90) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.085969) + 
    (-1.00 <= eta && eta < -0.90) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.108965) + 
    (-1.00 <= eta && eta < -0.90) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.132030) + 
    (-1.00 <= eta && eta < -0.90) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.154564) + 
    (-1.00 <= eta && eta < -0.90) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.176174) + 
    (-1.00 <= eta && eta < -0.90) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.198024) + 
    (-1.00 <= eta && eta < -0.90) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.219729) + 
    (-1.00 <= eta && eta < -0.90) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.240906) + 
    (-1.00 <= eta && eta < -0.90) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.262363) + 
    (-1.00 <= eta && eta < -0.90) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.283549) + 
    (-1.00 <= eta && eta < -0.90) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.304817) + 
    (-1.00 <= eta && eta < -0.90) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.326209) + 
    (-1.00 <= eta && eta < -0.90) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.347442) + 
    (-1.00 <= eta && eta < -0.90) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.368549) + 
    (-1.00 <= eta && eta < -0.90) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 30.70) * (0.389801) + 
    (-1.00 <= eta && eta < -0.90) * (30.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.411336) + 
    (-1.00 <= eta && eta < -0.90) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.90) * (0.433234) + 
    (-1.00 <= eta && eta < -0.90) * (40.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.453570) + 
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.000649) + 
    (-0.90 <= eta && eta < -0.80) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.020231) + 
    (-0.90 <= eta && eta < -0.80) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.041709) + 
    (-0.90 <= eta && eta < -0.80) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.063327) + 
    (-0.90 <= eta && eta < -0.80) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.085708) + 
    (-0.90 <= eta && eta < -0.80) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.107361) + 
    (-0.90 <= eta && eta < -0.80) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.129120) + 
    (-0.90 <= eta && eta < -0.80) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.152075) + 
    (-0.90 <= eta && eta < -0.80) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.174114) + 
    (-0.90 <= eta && eta < -0.80) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.196407) + 
    (-0.90 <= eta && eta < -0.80) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.218543) + 
    (-0.90 <= eta && eta < -0.80) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.240121) + 
    (-0.90 <= eta && eta < -0.80) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.261958) + 
    (-0.90 <= eta && eta < -0.80) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.283483) + 
    (-0.90 <= eta && eta < -0.80) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.305052) + 
    (-0.90 <= eta && eta < -0.80) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.326699) + 
    (-0.90 <= eta && eta < -0.80) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.348135) + 
    (-0.90 <= eta && eta < -0.80) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.369386) + 
    (-0.90 <= eta && eta < -0.80) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.390721) + 
    (-0.90 <= eta && eta < -0.80) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.412270) + 
    (-0.90 <= eta && eta < -0.80) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.434265) + 
    (-0.90 <= eta && eta < -0.80) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.454684) + 
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000698) + 
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.020156) + 
    (-0.80 <= eta && eta < -0.70) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.040719) + 
    (-0.80 <= eta && eta < -0.70) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.062503) + 
    (-0.80 <= eta && eta < -0.70) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.083481) + 
    (-0.80 <= eta && eta < -0.70) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.103705) + 
    (-0.80 <= eta && eta < -0.70) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.124110) + 
    (-0.80 <= eta && eta < -0.70) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.145883) + 
    (-0.80 <= eta && eta < -0.70) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.168517) + 
    (-0.80 <= eta && eta < -0.70) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.189997) + 
    (-0.80 <= eta && eta < -0.70) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.211530) + 
    (-0.80 <= eta && eta < -0.70) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.234017) + 
    (-0.80 <= eta && eta < -0.70) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.255655) + 
    (-0.80 <= eta && eta < -0.70) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.277170) + 
    (-0.80 <= eta && eta < -0.70) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.298995) + 
    (-0.80 <= eta && eta < -0.70) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.320333) + 
    (-0.80 <= eta && eta < -0.70) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.341878) + 
    (-0.80 <= eta && eta < -0.70) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.50) * (0.363611) + 
    (-0.80 <= eta && eta < -0.70) * (26.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.384994) + 
    (-0.80 <= eta && eta < -0.70) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.406381) + 
    (-0.80 <= eta && eta < -0.70) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.10) * (0.428249) + 
    (-0.80 <= eta && eta < -0.70) * (38.10 <= pt * cosh(eta) && pt * cosh(eta) < 47.90) * (0.450901) + 
    (-0.80 <= eta && eta < -0.70) * (47.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462545) + 
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000644) + 
    (-0.70 <= eta && eta < -0.60) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019416) + 
    (-0.70 <= eta && eta < -0.60) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.039800) + 
    (-0.70 <= eta && eta < -0.60) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.061565) + 
    (-0.70 <= eta && eta < -0.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.082602) + 
    (-0.70 <= eta && eta < -0.60) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.102923) + 
    (-0.70 <= eta && eta < -0.60) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.123448) + 
    (-0.70 <= eta && eta < -0.60) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.145363) + 
    (-0.70 <= eta && eta < -0.60) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.168148) + 
    (-0.70 <= eta && eta < -0.60) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.189770) + 
    (-0.70 <= eta && eta < -0.60) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.211440) + 
    (-0.70 <= eta && eta < -0.60) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.232802) + 
    (-0.70 <= eta && eta < -0.60) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.253505) + 
    (-0.70 <= eta && eta < -0.60) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.274359) + 
    (-0.70 <= eta && eta < -0.60) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.295736) + 
    (-0.70 <= eta && eta < -0.60) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.317652) + 
    (-0.70 <= eta && eta < -0.60) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.339116) + 
    (-0.70 <= eta && eta < -0.60) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.360427) + 
    (-0.70 <= eta && eta < -0.60) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.381782) + 
    (-0.70 <= eta && eta < -0.60) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.403198) + 
    (-0.70 <= eta && eta < -0.60) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.424964) + 
    (-0.70 <= eta && eta < -0.60) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.447490) + 
    (-0.70 <= eta && eta < -0.60) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461199) + 
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000623) + 
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019050) + 
    (-0.60 <= eta && eta < -0.50) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.039230) + 
    (-0.60 <= eta && eta < -0.50) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.060853) + 
    (-0.60 <= eta && eta < -0.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.081798) + 
    (-0.60 <= eta && eta < -0.50) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.102059) + 
    (-0.60 <= eta && eta < -0.50) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.122547) + 
    (-0.60 <= eta && eta < -0.50) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.144443) + 
    (-0.60 <= eta && eta < -0.50) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.167226) + 
    (-0.60 <= eta && eta < -0.50) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.188860) + 
    (-0.60 <= eta && eta < -0.50) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.210554) + 
    (-0.60 <= eta && eta < -0.50) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.231950) + 
    (-0.60 <= eta && eta < -0.50) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.252693) + 
    (-0.60 <= eta && eta < -0.50) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.273595) + 
    (-0.60 <= eta && eta < -0.50) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.295028) + 
    (-0.60 <= eta && eta < -0.50) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.317006) + 
    (-0.60 <= eta && eta < -0.50) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.338536) + 
    (-0.60 <= eta && eta < -0.50) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.359915) + 
    (-0.60 <= eta && eta < -0.50) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.381344) + 
    (-0.60 <= eta && eta < -0.50) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.402835) + 
    (-0.60 <= eta && eta < -0.50) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.424681) + 
    (-0.60 <= eta && eta < -0.50) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.60) * (0.447175) + 
    (-0.60 <= eta && eta < -0.50) * (45.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460966) + 
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000665) + 
    (-0.50 <= eta && eta < -0.40) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.020045) + 
    (-0.50 <= eta && eta < -0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040989) + 
    (-0.50 <= eta && eta < -0.40) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063254) + 
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084688) + 
    (-0.50 <= eta && eta < -0.40) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105322) + 
    (-0.50 <= eta && eta < -0.40) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.126101) + 
    (-0.50 <= eta && eta < -0.40) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146555) + 
    (-0.50 <= eta && eta < -0.40) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167958) + 
    (-0.50 <= eta && eta < -0.40) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189862) + 
    (-0.50 <= eta && eta < -0.40) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211794) + 
    (-0.50 <= eta && eta < -0.40) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233391) + 
    (-0.50 <= eta && eta < -0.40) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254294) + 
    (-0.50 <= eta && eta < -0.40) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275318) + 
    (-0.50 <= eta && eta < -0.40) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296835) + 
    (-0.50 <= eta && eta < -0.40) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318854) + 
    (-0.50 <= eta && eta < -0.40) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340377) + 
    (-0.50 <= eta && eta < -0.40) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361703) + 
    (-0.50 <= eta && eta < -0.40) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.383028) + 
    (-0.50 <= eta && eta < -0.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404363) + 
    (-0.50 <= eta && eta < -0.40) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426192) + 
    (-0.50 <= eta && eta < -0.40) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448789) + 
    (-0.50 <= eta && eta < -0.40) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462263) + 
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000645) + 
    (-0.40 <= eta && eta < -0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019306) + 
    (-0.40 <= eta && eta < -0.30) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.039824) + 
    (-0.40 <= eta && eta < -0.30) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.061801) + 
    (-0.40 <= eta && eta < -0.30) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.083051) + 
    (-0.40 <= eta && eta < -0.30) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.103567) + 
    (-0.40 <= eta && eta < -0.30) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.124275) + 
    (-0.40 <= eta && eta < -0.30) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.144695) + 
    (-0.40 <= eta && eta < -0.30) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.166094) + 
    (-0.40 <= eta && eta < -0.30) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.188023) + 
    (-0.40 <= eta && eta < -0.30) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.210005) + 
    (-0.40 <= eta && eta < -0.30) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.231670) + 
    (-0.40 <= eta && eta < -0.30) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.252656) + 
    (-0.40 <= eta && eta < -0.30) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.273778) + 
    (-0.40 <= eta && eta < -0.30) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.295407) + 
    (-0.40 <= eta && eta < -0.30) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.316782) + 
    (-0.40 <= eta && eta < -0.30) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.337879) + 
    (-0.40 <= eta && eta < -0.30) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.359062) + 
    (-0.40 <= eta && eta < -0.30) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.380434) + 
    (-0.40 <= eta && eta < -0.30) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.402053) + 
    (-0.40 <= eta && eta < -0.30) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.423693) + 
    (-0.40 <= eta && eta < -0.30) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.446017) + 
    (-0.40 <= eta && eta < -0.30) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460865) + 
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000624) + 
    (-0.30 <= eta && eta < -0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019209) + 
    (-0.30 <= eta && eta < -0.20) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.039736) + 
    (-0.30 <= eta && eta < -0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.061691) + 
    (-0.30 <= eta && eta < -0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.082927) + 
    (-0.30 <= eta && eta < -0.20) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.103434) + 
    (-0.30 <= eta && eta < -0.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.124136) + 
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.144554) + 
    (-0.30 <= eta && eta < -0.20) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.165953) + 
    (-0.30 <= eta && eta < -0.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.187883) + 
    (-0.30 <= eta && eta < -0.20) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.209869) + 
    (-0.30 <= eta && eta < -0.20) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.231539) + 
    (-0.30 <= eta && eta < -0.20) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.252531) + 
    (-0.30 <= eta && eta < -0.20) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.273661) + 
    (-0.30 <= eta && eta < -0.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.295298) + 
    (-0.30 <= eta && eta < -0.20) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.316682) + 
    (-0.30 <= eta && eta < -0.20) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.337790) + 
    (-0.30 <= eta && eta < -0.20) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.358983) + 
    (-0.30 <= eta && eta < -0.20) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.380366) + 
    (-0.30 <= eta && eta < -0.20) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.401997) + 
    (-0.30 <= eta && eta < -0.20) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.423648) + 
    (-0.30 <= eta && eta < -0.20) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.445986) + 
    (-0.30 <= eta && eta < -0.20) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460842) + 
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000667) + 
    (-0.20 <= eta && eta < -0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019658) + 
    (-0.20 <= eta && eta < -0.10) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.040176) + 
    (-0.20 <= eta && eta < -0.10) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.062033) + 
    (-0.20 <= eta && eta < -0.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.083131) + 
    (-0.20 <= eta && eta < -0.10) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.103490) + 
    (-0.20 <= eta && eta < -0.10) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.124039) + 
    (-0.20 <= eta && eta < -0.10) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.145966) + 
    (-0.20 <= eta && eta < -0.10) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.168751) + 
    (-0.20 <= eta && eta < -0.10) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.190365) + 
    (-0.20 <= eta && eta < -0.10) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.212019) + 
    (-0.20 <= eta && eta < -0.10) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.233359) + 
    (-0.20 <= eta && eta < -0.10) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.254035) + 
    (-0.20 <= eta && eta < -0.10) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.274858) + 
    (-0.20 <= eta && eta < -0.10) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.296199) + 
    (-0.20 <= eta && eta < -0.10) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.318075) + 
    (-0.20 <= eta && eta < -0.10) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.339496) + 
    (-0.20 <= eta && eta < -0.10) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.360760) + 
    (-0.20 <= eta && eta < -0.10) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.382069) + 
    (-0.20 <= eta && eta < -0.10) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.403435) + 
    (-0.20 <= eta && eta < -0.10) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.425149) + 
    (-0.20 <= eta && eta < -0.10) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.447621) + 
    (-0.20 <= eta && eta < -0.10) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461296) + 
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.000679) + 
    (-0.10 <= eta && eta < -0.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.019921) + 
    (-0.10 <= eta && eta < -0.00) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.039923) + 
    (-0.10 <= eta && eta < -0.00) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.061117) + 
    (-0.10 <= eta && eta < -0.00) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.083207) + 
    (-0.10 <= eta && eta < -0.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.104676) + 
    (-0.10 <= eta && eta < -0.00) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.126328) + 
    (-0.10 <= eta && eta < -0.00) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.149235) + 
    (-0.10 <= eta && eta < -0.00) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.171280) + 
    (-0.10 <= eta && eta < -0.00) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.192190) + 
    (-0.10 <= eta && eta < -0.00) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.213160) + 
    (-0.10 <= eta && eta < -0.00) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.235081) + 
    (-0.10 <= eta && eta < -0.00) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.257304) + 
    (-0.10 <= eta && eta < -0.00) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.279240) + 
    (-0.10 <= eta && eta < -0.00) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.300402) + 
    (-0.10 <= eta && eta < -0.00) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.321151) + 
    (-0.10 <= eta && eta < -0.00) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.342173) + 
    (-0.10 <= eta && eta < -0.00) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.363463) + 
    (-0.10 <= eta && eta < -0.00) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.384883) + 
    (-0.10 <= eta && eta < -0.00) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.406518) + 
    (-0.10 <= eta && eta < -0.00) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 38.90) * (0.428380) + 
    (-0.10 <= eta && eta < -0.00) * (38.90 <= pt * cosh(eta) && pt * cosh(eta) < 49.00) * (0.451065) + 
    (-0.10 <= eta && eta < -0.00) * (49.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461887) + 
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000652) + 
    (-0.00 <= eta && eta < 0.10) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.019853) + 
    (-0.00 <= eta && eta < 0.10) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041350) + 
    (-0.00 <= eta && eta < 0.10) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063083) + 
    (-0.00 <= eta && eta < 0.10) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.085621) + 
    (-0.00 <= eta && eta < 0.10) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.107436) + 
    (-0.00 <= eta && eta < 0.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.127685) + 
    (-0.00 <= eta && eta < 0.10) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.149221) + 
    (-0.00 <= eta && eta < 0.10) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.171557) + 
    (-0.00 <= eta && eta < 0.10) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.192724) + 
    (-0.00 <= eta && eta < 0.10) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.213928) + 
    (-0.00 <= eta && eta < 0.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.236064) + 
    (-0.00 <= eta && eta < 0.10) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.258471) + 
    (-0.00 <= eta && eta < 0.10) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.280550) + 
    (-0.00 <= eta && eta < 0.10) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.301813) + 
    (-0.00 <= eta && eta < 0.10) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.322623) + 
    (-0.00 <= eta && eta < 0.10) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.343665) + 
    (-0.00 <= eta && eta < 0.10) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.364932) + 
    (-0.00 <= eta && eta < 0.10) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.386281) + 
    (-0.00 <= eta && eta < 0.10) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.407794) + 
    (-0.00 <= eta && eta < 0.10) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 38.80) * (0.429658) + 
    (-0.00 <= eta && eta < 0.10) * (38.80 <= pt * cosh(eta) && pt * cosh(eta) < 49.20) * (0.452424) + 
    (-0.00 <= eta && eta < 0.10) * (49.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463042) + 
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000688) + 
    (0.10 <= eta && eta < 0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.020283) + 
    (0.10 <= eta && eta < 0.20) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.041390) + 
    (0.10 <= eta && eta < 0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063752) + 
    (0.10 <= eta && eta < 0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.085249) + 
    (0.10 <= eta && eta < 0.20) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105921) + 
    (0.10 <= eta && eta < 0.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.126724) + 
    (0.10 <= eta && eta < 0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.147189) + 
    (0.10 <= eta && eta < 0.20) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.168592) + 
    (0.10 <= eta && eta < 0.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.190487) + 
    (0.10 <= eta && eta < 0.20) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.212403) + 
    (0.10 <= eta && eta < 0.20) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233976) + 
    (0.10 <= eta && eta < 0.20) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254850) + 
    (0.10 <= eta && eta < 0.20) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275841) + 
    (0.10 <= eta && eta < 0.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.297319) + 
    (0.10 <= eta && eta < 0.20) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.319295) + 
    (0.10 <= eta && eta < 0.20) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340773) + 
    (0.10 <= eta && eta < 0.20) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.362051) + 
    (0.10 <= eta && eta < 0.20) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.383326) + 
    (0.10 <= eta && eta < 0.20) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404609) + 
    (0.10 <= eta && eta < 0.20) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426384) + 
    (0.10 <= eta && eta < 0.20) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.90) * (0.449033) + 
    (0.10 <= eta && eta < 0.20) * (45.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462443) + 
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000702) + 
    (0.20 <= eta && eta < 0.30) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.021415) + 
    (0.20 <= eta && eta < 0.30) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.041984) + 
    (0.20 <= eta && eta < 0.30) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.063242) + 
    (0.20 <= eta && eta < 0.30) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.085050) + 
    (0.20 <= eta && eta < 0.30) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.106045) + 
    (0.20 <= eta && eta < 0.30) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.127172) + 
    (0.20 <= eta && eta < 0.30) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.147941) + 
    (0.20 <= eta && eta < 0.30) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.169640) + 
    (0.20 <= eta && eta < 0.30) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.191806) + 
    (0.20 <= eta && eta < 0.30) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.213956) + 
    (0.20 <= eta && eta < 0.30) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.235721) + 
    (0.20 <= eta && eta < 0.30) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.256741) + 
    (0.20 <= eta && eta < 0.30) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.277837) + 
    (0.20 <= eta && eta < 0.30) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.299379) + 
    (0.20 <= eta && eta < 0.30) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.321372) + 
    (0.20 <= eta && eta < 0.30) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.342817) + 
    (0.20 <= eta && eta < 0.30) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.364014) + 
    (0.20 <= eta && eta < 0.30) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.385555) + 
    (0.20 <= eta && eta < 0.30) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.407162) + 
    (0.20 <= eta && eta < 0.30) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.429002) + 
    (0.20 <= eta && eta < 0.30) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.451798) + 
    (0.20 <= eta && eta < 0.30) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.464291) + 
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000726) + 
    (0.30 <= eta && eta < 0.40) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.020369) + 
    (0.30 <= eta && eta < 0.40) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.041656) + 
    (0.30 <= eta && eta < 0.40) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.064296) + 
    (0.30 <= eta && eta < 0.40) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.086047) + 
    (0.30 <= eta && eta < 0.40) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.106941) + 
    (0.30 <= eta && eta < 0.40) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.127942) + 
    (0.30 <= eta && eta < 0.40) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.148575) + 
    (0.30 <= eta && eta < 0.40) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.170125) + 
    (0.30 <= eta && eta < 0.40) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.192139) + 
    (0.30 <= eta && eta < 0.40) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.214143) + 
    (0.30 <= eta && eta < 0.40) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.235774) + 
    (0.30 <= eta && eta < 0.40) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.256674) + 
    (0.30 <= eta && eta < 0.40) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.277665) + 
    (0.30 <= eta && eta < 0.40) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.299113) + 
    (0.30 <= eta && eta < 0.40) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.321029) + 
    (0.30 <= eta && eta < 0.40) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.342418) + 
    (0.30 <= eta && eta < 0.40) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.363580) + 
    (0.30 <= eta && eta < 0.40) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.385108) + 
    (0.30 <= eta && eta < 0.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.406727) + 
    (0.30 <= eta && eta < 0.40) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.428608) + 
    (0.30 <= eta && eta < 0.40) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.451376) + 
    (0.30 <= eta && eta < 0.40) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463881) + 
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.000671) + 
    (0.40 <= eta && eta < 0.50) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.020318) + 
    (0.40 <= eta && eta < 0.50) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.040535) + 
    (0.40 <= eta && eta < 0.50) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.061663) + 
    (0.40 <= eta && eta < 0.50) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.083458) + 
    (0.40 <= eta && eta < 0.50) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.104511) + 
    (0.40 <= eta && eta < 0.50) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.125739) + 
    (0.40 <= eta && eta < 0.50) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.146635) + 
    (0.40 <= eta && eta < 0.50) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.168484) + 
    (0.40 <= eta && eta < 0.50) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.190813) + 
    (0.40 <= eta && eta < 0.50) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.213130) + 
    (0.40 <= eta && eta < 0.50) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.235057) + 
    (0.40 <= eta && eta < 0.50) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.256227) + 
    (0.40 <= eta && eta < 0.50) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.277468) + 
    (0.40 <= eta && eta < 0.50) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.299144) + 
    (0.40 <= eta && eta < 0.50) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.321260) + 
    (0.40 <= eta && eta < 0.50) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.342810) + 
    (0.40 <= eta && eta < 0.50) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.364092) + 
    (0.40 <= eta && eta < 0.50) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.385699) + 
    (0.40 <= eta && eta < 0.50) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.407349) + 
    (0.40 <= eta && eta < 0.50) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.429205) + 
    (0.40 <= eta && eta < 0.50) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.70) * (0.451986) + 
    (0.40 <= eta && eta < 0.50) * (46.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.464527) + 
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000629) + 
    (0.50 <= eta && eta < 0.60) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.019027) + 
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.039568) + 
    (0.50 <= eta && eta < 0.60) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.061693) + 
    (0.50 <= eta && eta < 0.60) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.083113) + 
    (0.50 <= eta && eta < 0.60) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.103800) + 
    (0.50 <= eta && eta < 0.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.124675) + 
    (0.50 <= eta && eta < 0.60) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.145250) + 
    (0.50 <= eta && eta < 0.60) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.166798) + 
    (0.50 <= eta && eta < 0.60) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.188860) + 
    (0.50 <= eta && eta < 0.60) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.210955) + 
    (0.50 <= eta && eta < 0.60) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.232711) + 
    (0.50 <= eta && eta < 0.60) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.253762) + 
    (0.50 <= eta && eta < 0.60) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.274929) + 
    (0.50 <= eta && eta < 0.60) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.296580) + 
    (0.50 <= eta && eta < 0.60) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.317952) + 
    (0.50 <= eta && eta < 0.60) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.339024) + 
    (0.50 <= eta && eta < 0.60) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.360155) + 
    (0.50 <= eta && eta < 0.60) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.381450) + 
    (0.50 <= eta && eta < 0.60) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.402963) + 
    (0.50 <= eta && eta < 0.60) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.424674) + 
    (0.50 <= eta && eta < 0.60) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.70) * (0.447121) + 
    (0.50 <= eta && eta < 0.60) * (44.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461667) + 
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000675) + 
    (0.60 <= eta && eta < 0.70) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019837) + 
    (0.60 <= eta && eta < 0.70) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040690) + 
    (0.60 <= eta && eta < 0.70) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.062882) + 
    (0.60 <= eta && eta < 0.70) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084270) + 
    (0.60 <= eta && eta < 0.70) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.104874) + 
    (0.60 <= eta && eta < 0.70) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.125635) + 
    (0.60 <= eta && eta < 0.70) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146081) + 
    (0.60 <= eta && eta < 0.70) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167483) + 
    (0.60 <= eta && eta < 0.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189393) + 
    (0.60 <= eta && eta < 0.70) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211339) + 
    (0.60 <= eta && eta < 0.70) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.232954) + 
    (0.60 <= eta && eta < 0.70) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.253877) + 
    (0.60 <= eta && eta < 0.70) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.274927) + 
    (0.60 <= eta && eta < 0.70) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296472) + 
    (0.60 <= eta && eta < 0.70) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318523) + 
    (0.60 <= eta && eta < 0.70) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340081) + 
    (0.60 <= eta && eta < 0.70) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361442) + 
    (0.60 <= eta && eta < 0.70) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.382805) + 
    (0.60 <= eta && eta < 0.70) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404178) + 
    (0.60 <= eta && eta < 0.70) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426049) + 
    (0.60 <= eta && eta < 0.70) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448688) + 
    (0.60 <= eta && eta < 0.70) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462189) + 
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000676) + 
    (0.70 <= eta && eta < 0.80) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019969) + 
    (0.70 <= eta && eta < 0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040896) + 
    (0.70 <= eta && eta < 0.80) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063138) + 
    (0.70 <= eta && eta < 0.80) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084558) + 
    (0.70 <= eta && eta < 0.80) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105182) + 
    (0.70 <= eta && eta < 0.80) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.125956) + 
    (0.70 <= eta && eta < 0.80) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146407) + 
    (0.70 <= eta && eta < 0.80) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167810) + 
    (0.70 <= eta && eta < 0.80) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189716) + 
    (0.70 <= eta && eta < 0.80) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211653) + 
    (0.70 <= eta && eta < 0.80) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233255) + 
    (0.70 <= eta && eta < 0.80) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254164) + 
    (0.70 <= eta && eta < 0.80) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275196) + 
    (0.70 <= eta && eta < 0.80) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296722) + 
    (0.70 <= eta && eta < 0.80) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318751) + 
    (0.70 <= eta && eta < 0.80) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340285) + 
    (0.70 <= eta && eta < 0.80) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361622) + 
    (0.70 <= eta && eta < 0.80) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.382958) + 
    (0.70 <= eta && eta < 0.80) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404306) + 
    (0.70 <= eta && eta < 0.80) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426148) + 
    (0.70 <= eta && eta < 0.80) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448757) + 
    (0.70 <= eta && eta < 0.80) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462240) + 
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000666) + 
    (0.80 <= eta && eta < 0.90) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.019973) + 
    (0.80 <= eta && eta < 0.90) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041538) + 
    (0.80 <= eta && eta < 0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063316) + 
    (0.80 <= eta && eta < 0.90) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.085883) + 
    (0.80 <= eta && eta < 0.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.107717) + 
    (0.80 <= eta && eta < 0.90) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.127976) + 
    (0.80 <= eta && eta < 0.90) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.149517) + 
    (0.80 <= eta && eta < 0.90) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.171852) + 
    (0.80 <= eta && eta < 0.90) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.193015) + 
    (0.80 <= eta && eta < 0.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.214211) + 
    (0.80 <= eta && eta < 0.90) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.236336) + 
    (0.80 <= eta && eta < 0.90) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.258728) + 
    (0.80 <= eta && eta < 0.90) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.280790) + 
    (0.80 <= eta && eta < 0.90) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.302035) + 
    (0.80 <= eta && eta < 0.90) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.322825) + 
    (0.80 <= eta && eta < 0.90) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.343847) + 
    (0.80 <= eta && eta < 0.90) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.365091) + 
    (0.80 <= eta && eta < 0.90) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.386417) + 
    (0.80 <= eta && eta < 0.90) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.407905) + 
    (0.80 <= eta && eta < 0.90) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 38.80) * (0.429743) + 
    (0.80 <= eta && eta < 0.90) * (38.80 <= pt * cosh(eta) && pt * cosh(eta) < 49.20) * (0.452482) + 
    (0.80 <= eta && eta < 0.90) * (49.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463087) + 
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000631) + 
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.020111) + 
    (0.90 <= eta && eta < 1.00) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041754) + 
    (0.90 <= eta && eta < 1.00) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063581) + 
    (0.90 <= eta && eta < 1.00) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.084516) + 
    (0.90 <= eta && eta < 1.00) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.104653) + 
    (0.90 <= eta && eta < 1.00) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.124945) + 
    (0.90 <= eta && eta < 1.00) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.146583) + 
    (0.90 <= eta && eta < 1.00) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.169067) + 
    (0.90 <= eta && eta < 1.00) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.190404) + 
    (0.90 <= eta && eta < 1.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.211798) + 
    (0.90 <= eta && eta < 1.00) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.234147) + 
    (0.90 <= eta && eta < 1.00) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.255662) + 
    (0.90 <= eta && eta < 1.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.277068) + 
    (0.90 <= eta && eta < 1.00) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.298797) + 
    (0.90 <= eta && eta < 1.00) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.320057) + 
    (0.90 <= eta && eta < 1.00) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.341541) + 
    (0.90 <= eta && eta < 1.00) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.363233) + 
    (0.90 <= eta && eta < 1.00) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.384597) + 
    (0.90 <= eta && eta < 1.00) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.80) * (0.405989) + 
    (0.90 <= eta && eta < 1.00) * (32.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.427888) + 
    (0.90 <= eta && eta < 1.00) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 48.00) * (0.450605) + 
    (0.90 <= eta && eta < 1.00) * (48.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462220) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 

  # --- pions ---

  add EfficiencyFormula {-211} {321} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.000629) + 
    (-1.00 <= eta && eta < -0.90) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.020818) + 
    (-1.00 <= eta && eta < -0.90) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.043238) + 
    (-1.00 <= eta && eta < -0.90) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.064987) + 
    (-1.00 <= eta && eta < -0.90) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.086643) + 
    (-1.00 <= eta && eta < -0.90) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.108246) + 
    (-1.00 <= eta && eta < -0.90) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.131038) + 
    (-1.00 <= eta && eta < -0.90) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.153265) + 
    (-1.00 <= eta && eta < -0.90) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.173948) + 
    (-1.00 <= eta && eta < -0.90) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.196026) + 
    (-1.00 <= eta && eta < -0.90) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.219035) + 
    (-1.00 <= eta && eta < -0.90) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.242086) + 
    (-1.00 <= eta && eta < -0.90) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.264737) + 
    (-1.00 <= eta && eta < -0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.286490) + 
    (-1.00 <= eta && eta < -0.90) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.308418) + 
    (-1.00 <= eta && eta < -0.90) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.329730) + 
    (-1.00 <= eta && eta < -0.90) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.350771) + 
    (-1.00 <= eta && eta < -0.90) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.372200) + 
    (-1.00 <= eta && eta < -0.90) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.393662) + 
    (-1.00 <= eta && eta < -0.90) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.415409) + 
    (-1.00 <= eta && eta < -0.90) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.437641) + 
    (-1.00 <= eta && eta < -0.90) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.461287) + 
    (-1.00 <= eta && eta < -0.90) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.480405) + 
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000601) + 
    (-0.90 <= eta && eta < -0.80) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.020257) + 
    (-0.90 <= eta && eta < -0.80) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.043042) + 
    (-0.90 <= eta && eta < -0.80) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.065328) + 
    (-0.90 <= eta && eta < -0.80) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.087115) + 
    (-0.90 <= eta && eta < -0.80) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.109928) + 
    (-0.90 <= eta && eta < -0.80) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.132996) + 
    (-0.90 <= eta && eta < -0.80) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.156649) + 
    (-0.90 <= eta && eta < -0.80) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.180600) + 
    (-0.90 <= eta && eta < -0.80) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.203148) + 
    (-0.90 <= eta && eta < -0.80) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.225705) + 
    (-0.90 <= eta && eta < -0.80) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.248439) + 
    (-0.90 <= eta && eta < -0.80) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.270714) + 
    (-0.90 <= eta && eta < -0.80) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.292054) + 
    (-0.90 <= eta && eta < -0.80) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.313521) + 
    (-0.90 <= eta && eta < -0.80) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.335443) + 
    (-0.90 <= eta && eta < -0.80) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.357624) + 
    (-0.90 <= eta && eta < -0.80) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.379405) + 
    (-0.90 <= eta && eta < -0.80) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.401051) + 
    (-0.90 <= eta && eta < -0.80) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.422825) + 
    (-0.90 <= eta && eta < -0.80) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.445280) + 
    (-0.90 <= eta && eta < -0.80) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 41.00) * (0.470610) + 
    (-0.90 <= eta && eta < -0.80) * (41.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484124) + 
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000757) + 
    (-0.80 <= eta && eta < -0.70) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.022320) + 
    (-0.80 <= eta && eta < -0.70) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.044193) + 
    (-0.80 <= eta && eta < -0.70) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.067162) + 
    (-0.80 <= eta && eta < -0.70) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.088655) + 
    (-0.80 <= eta && eta < -0.70) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.111401) + 
    (-0.80 <= eta && eta < -0.70) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.134165) + 
    (-0.80 <= eta && eta < -0.70) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.156920) + 
    (-0.80 <= eta && eta < -0.70) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.177718) + 
    (-0.80 <= eta && eta < -0.70) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.199864) + 
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.223340) + 
    (-0.80 <= eta && eta < -0.70) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.246613) + 
    (-0.80 <= eta && eta < -0.70) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.269390) + 
    (-0.80 <= eta && eta < -0.70) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.291176) + 
    (-0.80 <= eta && eta < -0.70) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.313049) + 
    (-0.80 <= eta && eta < -0.70) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.335333) + 
    (-0.80 <= eta && eta < -0.70) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.357820) + 
    (-0.80 <= eta && eta < -0.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.379834) + 
    (-0.80 <= eta && eta < -0.70) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.401637) + 
    (-0.80 <= eta && eta < -0.70) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.423486) + 
    (-0.80 <= eta && eta < -0.70) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.445919) + 
    (-0.80 <= eta && eta < -0.70) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 41.30) * (0.471428) + 
    (-0.80 <= eta && eta < -0.70) * (41.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484680) + 
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000627) + 
    (-0.70 <= eta && eta < -0.60) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.019743) + 
    (-0.70 <= eta && eta < -0.60) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.041036) + 
    (-0.70 <= eta && eta < -0.60) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.063090) + 
    (-0.70 <= eta && eta < -0.60) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.084632) + 
    (-0.70 <= eta && eta < -0.60) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.107944) + 
    (-0.70 <= eta && eta < -0.60) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.131352) + 
    (-0.70 <= eta && eta < -0.60) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.153268) + 
    (-0.70 <= eta && eta < -0.60) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.174594) + 
    (-0.70 <= eta && eta < -0.60) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.197491) + 
    (-0.70 <= eta && eta < -0.60) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.221492) + 
    (-0.70 <= eta && eta < -0.60) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.243056) + 
    (-0.70 <= eta && eta < -0.60) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.264487) + 
    (-0.70 <= eta && eta < -0.60) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.286988) + 
    (-0.70 <= eta && eta < -0.60) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.309562) + 
    (-0.70 <= eta && eta < -0.60) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.331387) + 
    (-0.70 <= eta && eta < -0.60) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.352810) + 
    (-0.70 <= eta && eta < -0.60) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374492) + 
    (-0.70 <= eta && eta < -0.60) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.396058) + 
    (-0.70 <= eta && eta < -0.60) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417745) + 
    (-0.70 <= eta && eta < -0.60) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.439965) + 
    (-0.70 <= eta && eta < -0.60) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.463862) + 
    (-0.70 <= eta && eta < -0.60) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481927) + 
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000606) + 
    (-0.60 <= eta && eta < -0.50) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.019510) + 
    (-0.60 <= eta && eta < -0.50) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.040196) + 
    (-0.60 <= eta && eta < -0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.061442) + 
    (-0.60 <= eta && eta < -0.50) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.083297) + 
    (-0.60 <= eta && eta < -0.50) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.105846) + 
    (-0.60 <= eta && eta < -0.50) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.129010) + 
    (-0.60 <= eta && eta < -0.50) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.151227) + 
    (-0.60 <= eta && eta < -0.50) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.173018) + 
    (-0.60 <= eta && eta < -0.50) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.197095) + 
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.220621) + 
    (-0.60 <= eta && eta < -0.50) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.242223) + 
    (-0.60 <= eta && eta < -0.50) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.263700) + 
    (-0.60 <= eta && eta < -0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.286256) + 
    (-0.60 <= eta && eta < -0.50) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.308892) + 
    (-0.60 <= eta && eta < -0.50) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.330782) + 
    (-0.60 <= eta && eta < -0.50) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.352274) + 
    (-0.60 <= eta && eta < -0.50) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374028) + 
    (-0.60 <= eta && eta < -0.50) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.395670) + 
    (-0.60 <= eta && eta < -0.50) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417435) + 
    (-0.60 <= eta && eta < -0.50) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.439737) + 
    (-0.60 <= eta && eta < -0.50) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.463725) + 
    (-0.60 <= eta && eta < -0.50) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481858) + 
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000754) + 
    (-0.50 <= eta && eta < -0.40) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.022057) + 
    (-0.50 <= eta && eta < -0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.044626) + 
    (-0.50 <= eta && eta < -0.40) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.067546) + 
    (-0.50 <= eta && eta < -0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.089482) + 
    (-0.50 <= eta && eta < -0.40) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.112586) + 
    (-0.50 <= eta && eta < -0.40) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.135083) + 
    (-0.50 <= eta && eta < -0.40) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.158756) + 
    (-0.50 <= eta && eta < -0.40) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.180504) + 
    (-0.50 <= eta && eta < -0.40) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.204102) + 
    (-0.50 <= eta && eta < -0.40) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.227320) + 
    (-0.50 <= eta && eta < -0.40) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.250611) + 
    (-0.50 <= eta && eta < -0.40) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.273341) + 
    (-0.50 <= eta && eta < -0.40) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.295022) + 
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.316735) + 
    (-0.50 <= eta && eta < -0.40) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.338798) + 
    (-0.50 <= eta && eta < -0.40) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.361005) + 
    (-0.50 <= eta && eta < -0.40) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.382690) + 
    (-0.50 <= eta && eta < -0.40) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.404114) + 
    (-0.50 <= eta && eta < -0.40) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425866) + 
    (-0.50 <= eta && eta < -0.40) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448463) + 
    (-0.50 <= eta && eta < -0.40) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.475401) + 
    (-0.50 <= eta && eta < -0.40) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486661) + 
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000728) + 
    (-0.40 <= eta && eta < -0.30) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.021373) + 
    (-0.40 <= eta && eta < -0.30) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.042787) + 
    (-0.40 <= eta && eta < -0.30) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.065124) + 
    (-0.40 <= eta && eta < -0.30) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.087923) + 
    (-0.40 <= eta && eta < -0.30) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.111671) + 
    (-0.40 <= eta && eta < -0.30) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.134974) + 
    (-0.40 <= eta && eta < -0.30) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.157718) + 
    (-0.40 <= eta && eta < -0.30) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.179824) + 
    (-0.40 <= eta && eta < -0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.202629) + 
    (-0.40 <= eta && eta < -0.30) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.225577) + 
    (-0.40 <= eta && eta < -0.30) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.248957) + 
    (-0.40 <= eta && eta < -0.30) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.271791) + 
    (-0.40 <= eta && eta < -0.30) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.293585) + 
    (-0.40 <= eta && eta < -0.30) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.315422) + 
    (-0.40 <= eta && eta < -0.30) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.337621) + 
    (-0.40 <= eta && eta < -0.30) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.359975) + 
    (-0.40 <= eta && eta < -0.30) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.381809) + 
    (-0.40 <= eta && eta < -0.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.403386) + 
    (-0.40 <= eta && eta < -0.30) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425299) + 
    (-0.40 <= eta && eta < -0.30) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448066) + 
    (-0.40 <= eta && eta < -0.30) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.474824) + 
    (-0.40 <= eta && eta < -0.30) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486349) + 
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000713) + 
    (-0.30 <= eta && eta < -0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.021188) + 
    (-0.30 <= eta && eta < -0.20) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.043291) + 
    (-0.30 <= eta && eta < -0.20) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.065636) + 
    (-0.30 <= eta && eta < -0.20) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.087346) + 
    (-0.30 <= eta && eta < -0.20) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.111189) + 
    (-0.30 <= eta && eta < -0.20) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.135124) + 
    (-0.30 <= eta && eta < -0.20) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.157289) + 
    (-0.30 <= eta && eta < -0.20) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.178343) + 
    (-0.30 <= eta && eta < -0.20) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.201490) + 
    (-0.30 <= eta && eta < -0.20) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.225445) + 
    (-0.30 <= eta && eta < -0.20) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.248831) + 
    (-0.30 <= eta && eta < -0.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.271673) + 
    (-0.30 <= eta && eta < -0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.293475) + 
    (-0.30 <= eta && eta < -0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.315321) + 
    (-0.30 <= eta && eta < -0.20) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.337531) + 
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.359896) + 
    (-0.30 <= eta && eta < -0.20) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.381742) + 
    (-0.30 <= eta && eta < -0.20) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.403330) + 
    (-0.30 <= eta && eta < -0.20) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.425256) + 
    (-0.30 <= eta && eta < -0.20) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.448036) + 
    (-0.30 <= eta && eta < -0.20) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.474810) + 
    (-0.30 <= eta && eta < -0.20) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.486341) + 
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000670) + 
    (-0.20 <= eta && eta < -0.10) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.020397) + 
    (-0.20 <= eta && eta < -0.10) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.041717) + 
    (-0.20 <= eta && eta < -0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.064356) + 
    (-0.20 <= eta && eta < -0.10) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.085635) + 
    (-0.20 <= eta && eta < -0.10) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.109110) + 
    (-0.20 <= eta && eta < -0.10) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.131602) + 
    (-0.20 <= eta && eta < -0.10) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.153486) + 
    (-0.20 <= eta && eta < -0.10) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.175586) + 
    (-0.20 <= eta && eta < -0.10) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.198903) + 
    (-0.20 <= eta && eta < -0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.222061) + 
    (-0.20 <= eta && eta < -0.10) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.243601) + 
    (-0.20 <= eta && eta < -0.10) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.265002) + 
    (-0.20 <= eta && eta < -0.10) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.287466) + 
    (-0.20 <= eta && eta < -0.10) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.310000) + 
    (-0.20 <= eta && eta < -0.10) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.331781) + 
    (-0.20 <= eta && eta < -0.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.353161) + 
    (-0.20 <= eta && eta < -0.10) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.374795) + 
    (-0.20 <= eta && eta < -0.10) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.396312) + 
    (-0.20 <= eta && eta < -0.10) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.417947) + 
    (-0.20 <= eta && eta < -0.10) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.440351) + 
    (-0.20 <= eta && eta < -0.10) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 35.00) * (0.464506) + 
    (-0.20 <= eta && eta < -0.10) * (35.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482178) + 
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.000742) + 
    (-0.10 <= eta && eta < -0.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.021676) + 
    (-0.10 <= eta && eta < -0.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.043000) + 
    (-0.10 <= eta && eta < -0.00) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.064909) + 
    (-0.10 <= eta && eta < -0.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.086553) + 
    (-0.10 <= eta && eta < -0.00) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.108414) + 
    (-0.10 <= eta && eta < -0.00) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.131561) + 
    (-0.10 <= eta && eta < -0.00) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.155629) + 
    (-0.10 <= eta && eta < -0.00) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.177181) + 
    (-0.10 <= eta && eta < -0.00) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.200013) + 
    (-0.10 <= eta && eta < -0.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.223039) + 
    (-0.10 <= eta && eta < -0.00) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.245902) + 
    (-0.10 <= eta && eta < -0.00) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.268329) + 
    (-0.10 <= eta && eta < -0.00) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.289835) + 
    (-0.10 <= eta && eta < -0.00) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.311487) + 
    (-0.10 <= eta && eta < -0.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.333613) + 
    (-0.10 <= eta && eta < -0.00) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.356016) + 
    (-0.10 <= eta && eta < -0.00) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.378024) + 
    (-0.10 <= eta && eta < -0.00) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.399906) + 
    (-0.10 <= eta && eta < -0.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.421924) + 
    (-0.10 <= eta && eta < -0.00) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.90) * (0.444432) + 
    (-0.10 <= eta && eta < -0.00) * (26.90 <= pt * cosh(eta) && pt * cosh(eta) < 40.00) * (0.469406) + 
    (-0.10 <= eta && eta < -0.00) * (40.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.483535) + 
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000700) + 
    (-0.00 <= eta && eta < 0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.021257) + 
    (-0.00 <= eta && eta < 0.10) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.043040) + 
    (-0.00 <= eta && eta < 0.10) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.064435) + 
    (-0.00 <= eta && eta < 0.10) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.086042) + 
    (-0.00 <= eta && eta < 0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.109607) + 
    (-0.00 <= eta && eta < 0.10) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.132930) + 
    (-0.00 <= eta && eta < 0.10) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.153931) + 
    (-0.00 <= eta && eta < 0.10) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.175705) + 
    (-0.00 <= eta && eta < 0.10) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.197877) + 
    (-0.00 <= eta && eta < 0.10) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220292) + 
    (-0.00 <= eta && eta < 0.10) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.241706) + 
    (-0.00 <= eta && eta < 0.10) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263021) + 
    (-0.00 <= eta && eta < 0.10) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285436) + 
    (-0.00 <= eta && eta < 0.10) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.307965) + 
    (-0.00 <= eta && eta < 0.10) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.329785) + 
    (-0.00 <= eta && eta < 0.10) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351244) + 
    (-0.00 <= eta && eta < 0.10) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373004) + 
    (-0.00 <= eta && eta < 0.10) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394692) + 
    (-0.00 <= eta && eta < 0.10) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416547) + 
    (-0.00 <= eta && eta < 0.10) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.438994) + 
    (-0.00 <= eta && eta < 0.10) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.462884) + 
    (-0.00 <= eta && eta < 0.10) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481363) + 
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000617) + 
    (0.10 <= eta && eta < 0.20) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.019051) + 
    (0.10 <= eta && eta < 0.20) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.040330) + 
    (0.10 <= eta && eta < 0.20) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.062862) + 
    (0.10 <= eta && eta < 0.20) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085532) + 
    (0.10 <= eta && eta < 0.20) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.107916) + 
    (0.10 <= eta && eta < 0.20) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.131693) + 
    (0.10 <= eta && eta < 0.20) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.155745) + 
    (0.10 <= eta && eta < 0.20) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.177142) + 
    (0.10 <= eta && eta < 0.20) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.196263) + 
    (0.10 <= eta && eta < 0.20) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.218643) + 
    (0.10 <= eta && eta < 0.20) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.240921) + 
    (0.10 <= eta && eta < 0.20) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.262860) + 
    (0.10 <= eta && eta < 0.20) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.285861) + 
    (0.10 <= eta && eta < 0.20) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.308892) + 
    (0.10 <= eta && eta < 0.20) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.331107) + 
    (0.10 <= eta && eta < 0.20) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.352856) + 
    (0.10 <= eta && eta < 0.20) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.374803) + 
    (0.10 <= eta && eta < 0.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.396560) + 
    (0.10 <= eta && eta < 0.20) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.418357) + 
    (0.10 <= eta && eta < 0.20) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.440831) + 
    (0.10 <= eta && eta < 0.20) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.465032) + 
    (0.10 <= eta && eta < 0.20) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482582) + 
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000736) + 
    (0.20 <= eta && eta < 0.30) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.021548) + 
    (0.20 <= eta && eta < 0.30) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.043947) + 
    (0.20 <= eta && eta < 0.30) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.067742) + 
    (0.20 <= eta && eta < 0.30) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.090454) + 
    (0.20 <= eta && eta < 0.30) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.113386) + 
    (0.20 <= eta && eta < 0.30) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.136843) + 
    (0.20 <= eta && eta < 0.30) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.160311) + 
    (0.20 <= eta && eta < 0.30) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.185362) + 
    (0.20 <= eta && eta < 0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.207279) + 
    (0.20 <= eta && eta < 0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.230256) + 
    (0.20 <= eta && eta < 0.30) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.253794) + 
    (0.20 <= eta && eta < 0.30) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.276522) + 
    (0.20 <= eta && eta < 0.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.298150) + 
    (0.20 <= eta && eta < 0.30) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.319758) + 
    (0.20 <= eta && eta < 0.30) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.341663) + 
    (0.20 <= eta && eta < 0.30) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.363659) + 
    (0.20 <= eta && eta < 0.30) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.385086) + 
    (0.20 <= eta && eta < 0.30) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.406698) + 
    (0.20 <= eta && eta < 0.30) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.428629) + 
    (0.20 <= eta && eta < 0.30) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.60) * (0.451291) + 
    (0.20 <= eta && eta < 0.30) * (27.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.478658) + 
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000753) + 
    (0.30 <= eta && eta < 0.40) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.021899) + 
    (0.30 <= eta && eta < 0.40) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.045558) + 
    (0.30 <= eta && eta < 0.40) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.068436) + 
    (0.30 <= eta && eta < 0.40) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.091463) + 
    (0.30 <= eta && eta < 0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.115435) + 
    (0.30 <= eta && eta < 0.40) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.139105) + 
    (0.30 <= eta && eta < 0.40) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.162598) + 
    (0.30 <= eta && eta < 0.40) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.182795) + 
    (0.30 <= eta && eta < 0.40) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.205458) + 
    (0.30 <= eta && eta < 0.40) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.228047) + 
    (0.30 <= eta && eta < 0.40) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.249350) + 
    (0.30 <= eta && eta < 0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.270616) + 
    (0.30 <= eta && eta < 0.40) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.292865) + 
    (0.30 <= eta && eta < 0.40) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.315110) + 
    (0.30 <= eta && eta < 0.40) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.336546) + 
    (0.30 <= eta && eta < 0.40) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.357525) + 
    (0.30 <= eta && eta < 0.40) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.378695) + 
    (0.30 <= eta && eta < 0.40) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.400236) + 
    (0.30 <= eta && eta < 0.40) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.422294) + 
    (0.30 <= eta && eta < 0.40) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.444843) + 
    (0.30 <= eta && eta < 0.40) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.469965) + 
    (0.30 <= eta && eta < 0.40) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.484558) + 
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.70) * (0.000629) + 
    (0.40 <= eta && eta < 0.50) * (5.70 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.019727) + 
    (0.40 <= eta && eta < 0.50) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.041772) + 
    (0.40 <= eta && eta < 0.50) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.064756) + 
    (0.40 <= eta && eta < 0.50) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.087848) + 
    (0.40 <= eta && eta < 0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.112141) + 
    (0.40 <= eta && eta < 0.50) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.135408) + 
    (0.40 <= eta && eta < 0.50) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.157790) + 
    (0.40 <= eta && eta < 0.50) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.180619) + 
    (0.40 <= eta && eta < 0.50) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.203357) + 
    (0.40 <= eta && eta < 0.50) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.227838) + 
    (0.40 <= eta && eta < 0.50) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.249584) + 
    (0.40 <= eta && eta < 0.50) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.271029) + 
    (0.40 <= eta && eta < 0.50) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.293437) + 
    (0.40 <= eta && eta < 0.50) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.315808) + 
    (0.40 <= eta && eta < 0.50) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.337334) + 
    (0.40 <= eta && eta < 0.50) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.358365) + 
    (0.40 <= eta && eta < 0.50) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.379551) + 
    (0.40 <= eta && eta < 0.50) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.401066) + 
    (0.40 <= eta && eta < 0.50) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.423052) + 
    (0.40 <= eta && eta < 0.50) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.445684) + 
    (0.40 <= eta && eta < 0.50) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.471163) + 
    (0.40 <= eta && eta < 0.50) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.485174) + 
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000656) + 
    (0.50 <= eta && eta < 0.60) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.020312) + 
    (0.50 <= eta && eta < 0.60) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.041999) + 
    (0.50 <= eta && eta < 0.60) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.064530) + 
    (0.50 <= eta && eta < 0.60) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.087046) + 
    (0.50 <= eta && eta < 0.60) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.110549) + 
    (0.50 <= eta && eta < 0.60) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.133710) + 
    (0.50 <= eta && eta < 0.60) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.156249) + 
    (0.50 <= eta && eta < 0.60) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.178337) + 
    (0.50 <= eta && eta < 0.60) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.199630) + 
    (0.50 <= eta && eta < 0.60) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.220540) + 
    (0.50 <= eta && eta < 0.60) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.242222) + 
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.264089) + 
    (0.50 <= eta && eta < 0.60) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.287002) + 
    (0.50 <= eta && eta < 0.60) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.309935) + 
    (0.50 <= eta && eta < 0.60) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.332047) + 
    (0.50 <= eta && eta < 0.60) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.353688) + 
    (0.50 <= eta && eta < 0.60) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.375520) + 
    (0.50 <= eta && eta < 0.60) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.397160) + 
    (0.50 <= eta && eta < 0.60) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.418834) + 
    (0.50 <= eta && eta < 0.60) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.441179) + 
    (0.50 <= eta && eta < 0.60) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.465437) + 
    (0.50 <= eta && eta < 0.60) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482784) + 
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.000610) + 
    (0.60 <= eta && eta < 0.70) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.019574) + 
    (0.60 <= eta && eta < 0.70) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.040602) + 
    (0.60 <= eta && eta < 0.70) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.063109) + 
    (0.60 <= eta && eta < 0.70) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085463) + 
    (0.60 <= eta && eta < 0.70) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.108771) + 
    (0.60 <= eta && eta < 0.70) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.132650) + 
    (0.60 <= eta && eta < 0.70) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.154383) + 
    (0.60 <= eta && eta < 0.70) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.175861) + 
    (0.60 <= eta && eta < 0.70) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.198600) + 
    (0.60 <= eta && eta < 0.70) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.222323) + 
    (0.60 <= eta && eta < 0.70) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.244096) + 
    (0.60 <= eta && eta < 0.70) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.265662) + 
    (0.60 <= eta && eta < 0.70) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.288269) + 
    (0.60 <= eta && eta < 0.70) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.310912) + 
    (0.60 <= eta && eta < 0.70) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.332765) + 
    (0.60 <= eta && eta < 0.70) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.354177) + 
    (0.60 <= eta && eta < 0.70) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.375807) + 
    (0.60 <= eta && eta < 0.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.397278) + 
    (0.60 <= eta && eta < 0.70) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.418823) + 
    (0.60 <= eta && eta < 0.70) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.441082) + 
    (0.60 <= eta && eta < 0.70) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.465310) + 
    (0.60 <= eta && eta < 0.70) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.482651) + 
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.000768) + 
    (0.70 <= eta && eta < 0.80) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.022252) + 
    (0.70 <= eta && eta < 0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.042122) + 
    (0.70 <= eta && eta < 0.80) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.062694) + 
    (0.70 <= eta && eta < 0.80) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.085610) + 
    (0.70 <= eta && eta < 0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.107607) + 
    (0.70 <= eta && eta < 0.80) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.130626) + 
    (0.70 <= eta && eta < 0.80) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.154667) + 
    (0.70 <= eta && eta < 0.80) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.176433) + 
    (0.70 <= eta && eta < 0.80) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.198617) + 
    (0.70 <= eta && eta < 0.80) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.222620) + 
    (0.70 <= eta && eta < 0.80) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.244391) + 
    (0.70 <= eta && eta < 0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.265940) + 
    (0.70 <= eta && eta < 0.80) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.288527) + 
    (0.70 <= eta && eta < 0.80) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.311148) + 
    (0.70 <= eta && eta < 0.80) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.332977) + 
    (0.70 <= eta && eta < 0.80) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.354366) + 
    (0.70 <= eta && eta < 0.80) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.375969) + 
    (0.70 <= eta && eta < 0.80) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.397969) + 
    (0.70 <= eta && eta < 0.80) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.420136) + 
    (0.70 <= eta && eta < 0.80) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.442635) + 
    (0.70 <= eta && eta < 0.80) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.60) * (0.467206) + 
    (0.70 <= eta && eta < 0.80) * (36.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.483337) + 
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000728) + 
    (0.80 <= eta && eta < 0.90) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.021330) + 
    (0.80 <= eta && eta < 0.90) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.042978) + 
    (0.80 <= eta && eta < 0.90) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.065797) + 
    (0.80 <= eta && eta < 0.90) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.087589) + 
    (0.80 <= eta && eta < 0.90) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.109904) + 
    (0.80 <= eta && eta < 0.90) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.132432) + 
    (0.80 <= eta && eta < 0.90) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.155030) + 
    (0.80 <= eta && eta < 0.90) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.177015) + 
    (0.80 <= eta && eta < 0.90) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.198109) + 
    (0.80 <= eta && eta < 0.90) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220572) + 
    (0.80 <= eta && eta < 0.90) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.241975) + 
    (0.80 <= eta && eta < 0.90) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263275) + 
    (0.80 <= eta && eta < 0.90) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285672) + 
    (0.80 <= eta && eta < 0.90) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.308182) + 
    (0.80 <= eta && eta < 0.90) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.329981) + 
    (0.80 <= eta && eta < 0.90) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351418) + 
    (0.80 <= eta && eta < 0.90) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373154) + 
    (0.80 <= eta && eta < 0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394818) + 
    (0.80 <= eta && eta < 0.90) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416648) + 
    (0.80 <= eta && eta < 0.90) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.439069) + 
    (0.80 <= eta && eta < 0.90) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.463037) + 
    (0.80 <= eta && eta < 0.90) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481439) + 
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.000687) + 
    (0.90 <= eta && eta < 1.00) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.020501) + 
    (0.90 <= eta && eta < 1.00) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.041727) + 
    (0.90 <= eta && eta < 1.00) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.063800) + 
    (0.90 <= eta && eta < 1.00) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.084994) + 
    (0.90 <= eta && eta < 1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.107765) + 
    (0.90 <= eta && eta < 1.00) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.130564) + 
    (0.90 <= eta && eta < 1.00) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.152720) + 
    (0.90 <= eta && eta < 1.00) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.174241) + 
    (0.90 <= eta && eta < 1.00) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.197247) + 
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.220892) + 
    (0.90 <= eta && eta < 1.00) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.242281) + 
    (0.90 <= eta && eta < 1.00) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.263564) + 
    (0.90 <= eta && eta < 1.00) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.285942) + 
    (0.90 <= eta && eta < 1.00) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.308428) + 
    (0.90 <= eta && eta < 1.00) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.330204) + 
    (0.90 <= eta && eta < 1.00) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.351616) + 
    (0.90 <= eta && eta < 1.00) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.373325) + 
    (0.90 <= eta && eta < 1.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.394961) + 
    (0.90 <= eta && eta < 1.00) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.416763) + 
    (0.90 <= eta && eta < 1.00) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.439153) + 
    (0.90 <= eta && eta < 1.00) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.463089) + 
    (0.90 <= eta && eta < 1.00) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.481465) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 

  add EfficiencyFormula {-211} {-211} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-1.00 <= eta && eta < -0.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.70) * (0.999363) + 
    (-1.00 <= eta && eta < -0.90) * (5.70 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.977469) + 
    (-1.00 <= eta && eta < -0.90) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.953824) + 
    (-1.00 <= eta && eta < -0.90) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.930396) + 
    (-1.00 <= eta && eta < -0.90) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.908335) + 
    (-1.00 <= eta && eta < -0.90) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.885172) + 
    (-1.00 <= eta && eta < -0.90) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.861520) + 
    (-1.00 <= eta && eta < -0.90) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.838567) + 
    (-1.00 <= eta && eta < -0.90) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.817634) + 
    (-1.00 <= eta && eta < -0.90) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.793978) + 
    (-1.00 <= eta && eta < -0.90) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.769187) + 
    (-1.00 <= eta && eta < -0.90) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.745592) + 
    (-1.00 <= eta && eta < -0.90) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.721481) + 
    (-1.00 <= eta && eta < -0.90) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.698201) + 
    (-1.00 <= eta && eta < -0.90) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.676260) + 
    (-1.00 <= eta && eta < -0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.654166) + 
    (-1.00 <= eta && eta < -0.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.629819) + 
    (-1.00 <= eta && eta < -0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604030) + 
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.90 <= eta && eta < -0.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.60) * (0.999415) + 
    (-0.90 <= eta && eta < -0.80) * (5.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.978320) + 
    (-0.90 <= eta && eta < -0.80) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.954946) + 
    (-0.90 <= eta && eta < -0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.931253) + 
    (-0.90 <= eta && eta < -0.80) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.908742) + 
    (-0.90 <= eta && eta < -0.80) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.884654) + 
    (-0.90 <= eta && eta < -0.80) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.860913) + 
    (-0.90 <= eta && eta < -0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.837747) + 
    (-0.90 <= eta && eta < -0.80) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.815568) + 
    (-0.90 <= eta && eta < -0.80) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.791569) + 
    (-0.90 <= eta && eta < -0.80) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.767848) + 
    (-0.90 <= eta && eta < -0.80) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.747406) + 
    (-0.90 <= eta && eta < -0.80) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.723516) + 
    (-0.90 <= eta && eta < -0.80) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.699697) + 
    (-0.90 <= eta && eta < -0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.677096) + 
    (-0.90 <= eta && eta < -0.80) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.654458) + 
    (-0.90 <= eta && eta < -0.80) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.630047) + 
    (-0.90 <= eta && eta < -0.80) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603999) + 
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.80 <= eta && eta < -0.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999438) + 
    (-0.80 <= eta && eta < -0.70) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.978782) + 
    (-0.80 <= eta && eta < -0.70) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.953835) + 
    (-0.80 <= eta && eta < -0.70) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.929657) + 
    (-0.80 <= eta && eta < -0.70) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.907134) + 
    (-0.80 <= eta && eta < -0.70) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.882151) + 
    (-0.80 <= eta && eta < -0.70) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.858132) + 
    (-0.80 <= eta && eta < -0.70) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.835697) + 
    (-0.80 <= eta && eta < -0.70) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.813827) + 
    (-0.80 <= eta && eta < -0.70) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.790188) + 
    (-0.80 <= eta && eta < -0.70) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.765614) + 
    (-0.80 <= eta && eta < -0.70) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.742033) + 
    (-0.80 <= eta && eta < -0.70) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.720207) + 
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.698161) + 
    (-0.80 <= eta && eta < -0.70) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.675683) + 
    (-0.80 <= eta && eta < -0.70) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.653272) + 
    (-0.80 <= eta && eta < -0.70) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.628628) + 
    (-0.80 <= eta && eta < -0.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603794) + 
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.70 <= eta && eta < -0.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999386) + 
    (-0.70 <= eta && eta < -0.60) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.977686) + 
    (-0.70 <= eta && eta < -0.60) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.952568) + 
    (-0.70 <= eta && eta < -0.60) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.927697) + 
    (-0.70 <= eta && eta < -0.60) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.904394) + 
    (-0.70 <= eta && eta < -0.60) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.880318) + 
    (-0.70 <= eta && eta < -0.60) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.856544) + 
    (-0.70 <= eta && eta < -0.60) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.833257) + 
    (-0.70 <= eta && eta < -0.60) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.810933) + 
    (-0.70 <= eta && eta < -0.60) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.786986) + 
    (-0.70 <= eta && eta < -0.60) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.763520) + 
    (-0.70 <= eta && eta < -0.60) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.740261) + 
    (-0.70 <= eta && eta < -0.60) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.718074) + 
    (-0.70 <= eta && eta < -0.60) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.696228) + 
    (-0.70 <= eta && eta < -0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.674142) + 
    (-0.70 <= eta && eta < -0.60) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.651446) + 
    (-0.70 <= eta && eta < -0.60) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.625653) + 
    (-0.70 <= eta && eta < -0.60) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603457) + 
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.60 <= eta && eta < -0.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999412) + 
    (-0.60 <= eta && eta < -0.50) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.978142) + 
    (-0.60 <= eta && eta < -0.50) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.953041) + 
    (-0.60 <= eta && eta < -0.50) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.928653) + 
    (-0.60 <= eta && eta < -0.50) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.906049) + 
    (-0.60 <= eta && eta < -0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.882734) + 
    (-0.60 <= eta && eta < -0.50) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.858632) + 
    (-0.60 <= eta && eta < -0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.835440) + 
    (-0.60 <= eta && eta < -0.50) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.813533) + 
    (-0.60 <= eta && eta < -0.50) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.789624) + 
    (-0.60 <= eta && eta < -0.50) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.765618) + 
    (-0.60 <= eta && eta < -0.50) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.744122) + 
    (-0.60 <= eta && eta < -0.50) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.722552) + 
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.700157) + 
    (-0.60 <= eta && eta < -0.50) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.677085) + 
    (-0.60 <= eta && eta < -0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.654123) + 
    (-0.60 <= eta && eta < -0.50) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.629276) + 
    (-0.60 <= eta && eta < -0.50) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603840) + 
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.50 <= eta && eta < -0.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999451) + 
    (-0.50 <= eta && eta < -0.40) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979067) + 
    (-0.50 <= eta && eta < -0.40) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.958187) + 
    (-0.50 <= eta && eta < -0.40) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.937021) + 
    (-0.50 <= eta && eta < -0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.913728) + 
    (-0.50 <= eta && eta < -0.40) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.889767) + 
    (-0.50 <= eta && eta < -0.40) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.865592) + 
    (-0.50 <= eta && eta < -0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.842135) + 
    (-0.50 <= eta && eta < -0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.819614) + 
    (-0.50 <= eta && eta < -0.40) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.795484) + 
    (-0.50 <= eta && eta < -0.40) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.770150) + 
    (-0.50 <= eta && eta < -0.40) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.747397) + 
    (-0.50 <= eta && eta < -0.40) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.725362) + 
    (-0.50 <= eta && eta < -0.40) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.703515) + 
    (-0.50 <= eta && eta < -0.40) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.681967) + 
    (-0.50 <= eta && eta < -0.40) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.659657) + 
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.636152) + 
    (-0.50 <= eta && eta < -0.40) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604442) + 
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.40 <= eta && eta < -0.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999472) + 
    (-0.40 <= eta && eta < -0.30) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979617) + 
    (-0.40 <= eta && eta < -0.30) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.955630) + 
    (-0.40 <= eta && eta < -0.30) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.932152) + 
    (-0.40 <= eta && eta < -0.30) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.909173) + 
    (-0.40 <= eta && eta < -0.30) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.885363) + 
    (-0.40 <= eta && eta < -0.30) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.860582) + 
    (-0.40 <= eta && eta < -0.30) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.836574) + 
    (-0.40 <= eta && eta < -0.30) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.813641) + 
    (-0.40 <= eta && eta < -0.30) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.789231) + 
    (-0.40 <= eta && eta < -0.30) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.764564) + 
    (-0.40 <= eta && eta < -0.30) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.740178) + 
    (-0.40 <= eta && eta < -0.30) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715912) + 
    (-0.40 <= eta && eta < -0.30) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.693289) + 
    (-0.40 <= eta && eta < -0.30) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.670721) + 
    (-0.40 <= eta && eta < -0.30) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647669) + 
    (-0.40 <= eta && eta < -0.30) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605441) + 
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.30 <= eta && eta < -0.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999491) + 
    (-0.30 <= eta && eta < -0.20) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.979938) + 
    (-0.30 <= eta && eta < -0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.959867) + 
    (-0.30 <= eta && eta < -0.20) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.938816) + 
    (-0.30 <= eta && eta < -0.20) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.916434) + 
    (-0.30 <= eta && eta < -0.20) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.892350) + 
    (-0.30 <= eta && eta < -0.20) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.868247) + 
    (-0.30 <= eta && eta < -0.20) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.845108) + 
    (-0.30 <= eta && eta < -0.20) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.821653) + 
    (-0.30 <= eta && eta < -0.20) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.799669) + 
    (-0.30 <= eta && eta < -0.20) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.776530) + 
    (-0.30 <= eta && eta < -0.20) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.754283) + 
    (-0.30 <= eta && eta < -0.20) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.732042) + 
    (-0.30 <= eta && eta < -0.20) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.708971) + 
    (-0.30 <= eta && eta < -0.20) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.685968) + 
    (-0.30 <= eta && eta < -0.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.663406) + 
    (-0.30 <= eta && eta < -0.20) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.640318) + 
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604825) + 
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.20 <= eta && eta < -0.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999335) + 
    (-0.20 <= eta && eta < -0.10) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.976774) + 
    (-0.20 <= eta && eta < -0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.951314) + 
    (-0.20 <= eta && eta < -0.10) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.926371) + 
    (-0.20 <= eta && eta < -0.10) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.902907) + 
    (-0.20 <= eta && eta < -0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.878307) + 
    (-0.20 <= eta && eta < -0.10) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.855267) + 
    (-0.20 <= eta && eta < -0.10) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.831969) + 
    (-0.20 <= eta && eta < -0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.806045) + 
    (-0.20 <= eta && eta < -0.10) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.780629) + 
    (-0.20 <= eta && eta < -0.10) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.757736) + 
    (-0.20 <= eta && eta < -0.10) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.734043) + 
    (-0.20 <= eta && eta < -0.10) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.711369) + 
    (-0.20 <= eta && eta < -0.10) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.688232) + 
    (-0.20 <= eta && eta < -0.10) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.665130) + 
    (-0.20 <= eta && eta < -0.10) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.641972) + 
    (-0.20 <= eta && eta < -0.10) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605010) + 
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.10 <= eta && eta < -0.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.60) * (0.999436) + 
    (-0.10 <= eta && eta < -0.00) * (5.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.978964) + 
    (-0.10 <= eta && eta < -0.00) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.958369) + 
    (-0.10 <= eta && eta < -0.00) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.938389) + 
    (-0.10 <= eta && eta < -0.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.916194) + 
    (-0.10 <= eta && eta < -0.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.893059) + 
    (-0.10 <= eta && eta < -0.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.869718) + 
    (-0.10 <= eta && eta < -0.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.846143) + 
    (-0.10 <= eta && eta < -0.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.825154) + 
    (-0.10 <= eta && eta < -0.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.800818) + 
    (-0.10 <= eta && eta < -0.00) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.778480) + 
    (-0.10 <= eta && eta < -0.00) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.755699) + 
    (-0.10 <= eta && eta < -0.00) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.733208) + 
    (-0.10 <= eta && eta < -0.00) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.710395) + 
    (-0.10 <= eta && eta < -0.00) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.687940) + 
    (-0.10 <= eta && eta < -0.00) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.665437) + 
    (-0.10 <= eta && eta < -0.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.642619) + 
    (-0.10 <= eta && eta < -0.00) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605218) + 
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (-0.00 <= eta && eta < 0.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999487) + 
    (-0.00 <= eta && eta < 0.10) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.980086) + 
    (-0.00 <= eta && eta < 0.10) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.955878) + 
    (-0.00 <= eta && eta < 0.10) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.931602) + 
    (-0.00 <= eta && eta < 0.10) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.908797) + 
    (-0.00 <= eta && eta < 0.10) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.885676) + 
    (-0.00 <= eta && eta < 0.10) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.862280) + 
    (-0.00 <= eta && eta < 0.10) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.839159) + 
    (-0.00 <= eta && eta < 0.10) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.815846) + 
    (-0.00 <= eta && eta < 0.10) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.795129) + 
    (-0.00 <= eta && eta < 0.10) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.773860) + 
    (-0.00 <= eta && eta < 0.10) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.751980) + 
    (-0.00 <= eta && eta < 0.10) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.728013) + 
    (-0.00 <= eta && eta < 0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.705518) + 
    (-0.00 <= eta && eta < 0.10) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.683919) + 
    (-0.00 <= eta && eta < 0.10) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.661467) + 
    (-0.00 <= eta && eta < 0.10) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.638022) + 
    (-0.00 <= eta && eta < 0.10) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604731) + 
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.10 <= eta && eta < 0.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999390) + 
    (0.10 <= eta && eta < 0.20) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978436) + 
    (0.10 <= eta && eta < 0.20) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.954026) + 
    (0.10 <= eta && eta < 0.20) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.928980) + 
    (0.10 <= eta && eta < 0.20) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.905261) + 
    (0.10 <= eta && eta < 0.20) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.880827) + 
    (0.10 <= eta && eta < 0.20) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.855798) + 
    (0.10 <= eta && eta < 0.20) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.831816) + 
    (0.10 <= eta && eta < 0.20) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.811017) + 
    (0.10 <= eta && eta < 0.20) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.790324) + 
    (0.10 <= eta && eta < 0.20) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.766845) + 
    (0.10 <= eta && eta < 0.20) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.742732) + 
    (0.10 <= eta && eta < 0.20) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.721120) + 
    (0.10 <= eta && eta < 0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.698097) + 
    (0.10 <= eta && eta < 0.20) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.675216) + 
    (0.10 <= eta && eta < 0.20) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.652610) + 
    (0.10 <= eta && eta < 0.20) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.627707) + 
    (0.10 <= eta && eta < 0.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603627) + 
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.20 <= eta && eta < 0.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999486) + 
    (0.20 <= eta && eta < 0.30) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.979519) + 
    (0.20 <= eta && eta < 0.30) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.955406) + 
    (0.20 <= eta && eta < 0.30) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.930947) + 
    (0.20 <= eta && eta < 0.30) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.906986) + 
    (0.20 <= eta && eta < 0.30) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.882189) + 
    (0.20 <= eta && eta < 0.30) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.857189) + 
    (0.20 <= eta && eta < 0.30) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.833268) + 
    (0.20 <= eta && eta < 0.30) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.811348) + 
    (0.20 <= eta && eta < 0.30) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.787352) + 
    (0.20 <= eta && eta < 0.30) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.80) * (0.762475) + 
    (0.20 <= eta && eta < 0.30) * (8.80 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.737585) + 
    (0.20 <= eta && eta < 0.30) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.716144) + 
    (0.20 <= eta && eta < 0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.693021) + 
    (0.20 <= eta && eta < 0.30) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.670061) + 
    (0.20 <= eta && eta < 0.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.646891) + 
    (0.20 <= eta && eta < 0.30) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605262) + 
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.30 <= eta && eta < 0.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999459) + 
    (0.30 <= eta && eta < 0.40) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.979305) + 
    (0.30 <= eta && eta < 0.40) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.955111) + 
    (0.30 <= eta && eta < 0.40) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.929350) + 
    (0.30 <= eta && eta < 0.40) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.904229) + 
    (0.30 <= eta && eta < 0.40) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.880114) + 
    (0.30 <= eta && eta < 0.40) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.855471) + 
    (0.30 <= eta && eta < 0.40) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.831596) + 
    (0.30 <= eta && eta < 0.40) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.00) * (0.808837) + 
    (0.30 <= eta && eta < 0.40) * (8.00 <= pt * cosh(eta) && pt * cosh(eta) < 8.40) * (0.785228) + 
    (0.30 <= eta && eta < 0.40) * (8.40 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.756982) + 
    (0.30 <= eta && eta < 0.40) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.732645) + 
    (0.30 <= eta && eta < 0.40) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.709880) + 
    (0.30 <= eta && eta < 0.40) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.687128) + 
    (0.30 <= eta && eta < 0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.663925) + 
    (0.30 <= eta && eta < 0.40) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.640354) + 
    (0.30 <= eta && eta < 0.40) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.604740) + 
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.40 <= eta && eta < 0.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.30) * (0.999416) + 
    (0.40 <= eta && eta < 0.50) * (5.30 <= pt * cosh(eta) && pt * cosh(eta) < 5.80) * (0.977866) + 
    (0.40 <= eta && eta < 0.50) * (5.80 <= pt * cosh(eta) && pt * cosh(eta) < 6.10) * (0.956023) + 
    (0.40 <= eta && eta < 0.50) * (6.10 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.934804) + 
    (0.40 <= eta && eta < 0.50) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.910945) + 
    (0.40 <= eta && eta < 0.50) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.886157) + 
    (0.40 <= eta && eta < 0.50) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.861028) + 
    (0.40 <= eta && eta < 0.50) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.836587) + 
    (0.40 <= eta && eta < 0.50) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.813125) + 
    (0.40 <= eta && eta < 0.50) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.30) * (0.788762) + 
    (0.40 <= eta && eta < 0.50) * (8.30 <= pt * cosh(eta) && pt * cosh(eta) < 8.70) * (0.764916) + 
    (0.40 <= eta && eta < 0.50) * (8.70 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.739594) + 
    (0.40 <= eta && eta < 0.50) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.715349) + 
    (0.40 <= eta && eta < 0.50) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.691540) + 
    (0.40 <= eta && eta < 0.50) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.668971) + 
    (0.40 <= eta && eta < 0.50) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.645517) + 
    (0.40 <= eta && eta < 0.50) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605109) + 
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.50 <= eta && eta < 0.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999374) + 
    (0.50 <= eta && eta < 0.60) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.976949) + 
    (0.50 <= eta && eta < 0.60) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.951356) + 
    (0.50 <= eta && eta < 0.60) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.926120) + 
    (0.50 <= eta && eta < 0.60) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.902566) + 
    (0.50 <= eta && eta < 0.60) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.878135) + 
    (0.50 <= eta && eta < 0.60) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.853567) + 
    (0.50 <= eta && eta < 0.60) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.830007) + 
    (0.50 <= eta && eta < 0.60) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.807715) + 
    (0.50 <= eta && eta < 0.60) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.784361) + 
    (0.50 <= eta && eta < 0.60) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.760491) + 
    (0.50 <= eta && eta < 0.60) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.736736) + 
    (0.50 <= eta && eta < 0.60) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.712414) + 
    (0.50 <= eta && eta < 0.60) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.690695) + 
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.668691) + 
    (0.50 <= eta && eta < 0.60) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.645647) + 
    (0.50 <= eta && eta < 0.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605221) + 
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.60 <= eta && eta < 0.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999428) + 
    (0.60 <= eta && eta < 0.70) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978093) + 
    (0.60 <= eta && eta < 0.70) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.30) * (0.952896) + 
    (0.60 <= eta && eta < 0.70) * (6.30 <= pt * cosh(eta) && pt * cosh(eta) < 6.60) * (0.928743) + 
    (0.60 <= eta && eta < 0.70) * (6.60 <= pt * cosh(eta) && pt * cosh(eta) < 6.90) * (0.904835) + 
    (0.60 <= eta && eta < 0.70) * (6.90 <= pt * cosh(eta) && pt * cosh(eta) < 7.20) * (0.880239) + 
    (0.60 <= eta && eta < 0.70) * (7.20 <= pt * cosh(eta) && pt * cosh(eta) < 7.50) * (0.855800) + 
    (0.60 <= eta && eta < 0.70) * (7.50 <= pt * cosh(eta) && pt * cosh(eta) < 7.80) * (0.832129) + 
    (0.60 <= eta && eta < 0.70) * (7.80 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.809971) + 
    (0.60 <= eta && eta < 0.70) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.785533) + 
    (0.60 <= eta && eta < 0.70) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.762384) + 
    (0.60 <= eta && eta < 0.70) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.739140) + 
    (0.60 <= eta && eta < 0.70) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715270) + 
    (0.60 <= eta && eta < 0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.692367) + 
    (0.60 <= eta && eta < 0.70) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.669998) + 
    (0.60 <= eta && eta < 0.70) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647168) + 
    (0.60 <= eta && eta < 0.70) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605392) + 
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.70 <= eta && eta < 0.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.40) * (0.999431) + 
    (0.70 <= eta && eta < 0.80) * (5.40 <= pt * cosh(eta) && pt * cosh(eta) < 5.90) * (0.978854) + 
    (0.70 <= eta && eta < 0.80) * (5.90 <= pt * cosh(eta) && pt * cosh(eta) < 6.20) * (0.957614) + 
    (0.70 <= eta && eta < 0.80) * (6.20 <= pt * cosh(eta) && pt * cosh(eta) < 6.50) * (0.936746) + 
    (0.70 <= eta && eta < 0.80) * (6.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.80) * (0.914572) + 
    (0.70 <= eta && eta < 0.80) * (6.80 <= pt * cosh(eta) && pt * cosh(eta) < 7.10) * (0.889412) + 
    (0.70 <= eta && eta < 0.80) * (7.10 <= pt * cosh(eta) && pt * cosh(eta) < 7.40) * (0.863632) + 
    (0.70 <= eta && eta < 0.80) * (7.40 <= pt * cosh(eta) && pt * cosh(eta) < 7.70) * (0.839744) + 
    (0.70 <= eta && eta < 0.80) * (7.70 <= pt * cosh(eta) && pt * cosh(eta) < 8.10) * (0.814966) + 
    (0.70 <= eta && eta < 0.80) * (8.10 <= pt * cosh(eta) && pt * cosh(eta) < 8.50) * (0.787803) + 
    (0.70 <= eta && eta < 0.80) * (8.50 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.762097) + 
    (0.70 <= eta && eta < 0.80) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.738550) + 
    (0.70 <= eta && eta < 0.80) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.715206) + 
    (0.70 <= eta && eta < 0.80) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.692151) + 
    (0.70 <= eta && eta < 0.80) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.669829) + 
    (0.70 <= eta && eta < 0.80) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.647051) + 
    (0.70 <= eta && eta < 0.80) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605381) + 
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.80 <= eta && eta < 0.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999463) + 
    (0.80 <= eta && eta < 0.90) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.979480) + 
    (0.80 <= eta && eta < 0.90) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.955658) + 
    (0.80 <= eta && eta < 0.90) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.931795) + 
    (0.80 <= eta && eta < 0.90) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.909021) + 
    (0.80 <= eta && eta < 0.90) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.884511) + 
    (0.80 <= eta && eta < 0.90) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.859876) + 
    (0.80 <= eta && eta < 0.90) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.837206) + 
    (0.80 <= eta && eta < 0.90) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.815556) + 
    (0.80 <= eta && eta < 0.90) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.792098) + 
    (0.80 <= eta && eta < 0.90) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.767572) + 
    (0.80 <= eta && eta < 0.90) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.742973) + 
    (0.80 <= eta && eta < 0.90) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.719838) + 
    (0.80 <= eta && eta < 0.90) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.696957) + 
    (0.80 <= eta && eta < 0.90) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.674906) + 
    (0.80 <= eta && eta < 0.90) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.652868) + 
    (0.80 <= eta && eta < 0.90) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.628233) + 
    (0.80 <= eta && eta < 0.90) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603774) + 
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (0.90 <= eta && eta < 1.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 5.50) * (0.999490) + 
    (0.90 <= eta && eta < 1.00) * (5.50 <= pt * cosh(eta) && pt * cosh(eta) < 6.00) * (0.980505) + 
    (0.90 <= eta && eta < 1.00) * (6.00 <= pt * cosh(eta) && pt * cosh(eta) < 6.40) * (0.957300) + 
    (0.90 <= eta && eta < 1.00) * (6.40 <= pt * cosh(eta) && pt * cosh(eta) < 6.70) * (0.933891) + 
    (0.90 <= eta && eta < 1.00) * (6.70 <= pt * cosh(eta) && pt * cosh(eta) < 7.00) * (0.910870) + 
    (0.90 <= eta && eta < 1.00) * (7.00 <= pt * cosh(eta) && pt * cosh(eta) < 7.30) * (0.887137) + 
    (0.90 <= eta && eta < 1.00) * (7.30 <= pt * cosh(eta) && pt * cosh(eta) < 7.60) * (0.863393) + 
    (0.90 <= eta && eta < 1.00) * (7.60 <= pt * cosh(eta) && pt * cosh(eta) < 7.90) * (0.840634) + 
    (0.90 <= eta && eta < 1.00) * (7.90 <= pt * cosh(eta) && pt * cosh(eta) < 8.20) * (0.818273) + 
    (0.90 <= eta && eta < 1.00) * (8.20 <= pt * cosh(eta) && pt * cosh(eta) < 8.60) * (0.794422) + 
    (0.90 <= eta && eta < 1.00) * (8.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.769639) + 
    (0.90 <= eta && eta < 1.00) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.748157) + 
    (0.90 <= eta && eta < 1.00) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.724553) + 
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.699984) + 
    (0.90 <= eta && eta < 1.00) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.677138) + 
    (0.90 <= eta && eta < 1.00) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.654326) + 
    (0.90 <= eta && eta < 1.00) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.629520) + 
    (0.90 <= eta && eta < 1.00) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.603898) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 


  # --- protons ---

  add EfficiencyFormula {2212} {321} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.000705) + 
    (-1.00 <= eta && eta < -0.90) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.021025) + 
    (-1.00 <= eta && eta < -0.90) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.042481) + 
    (-1.00 <= eta && eta < -0.90) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.063884) + 
    (-1.00 <= eta && eta < -0.90) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.085969) + 
    (-1.00 <= eta && eta < -0.90) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.108965) + 
    (-1.00 <= eta && eta < -0.90) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.132030) + 
    (-1.00 <= eta && eta < -0.90) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.154564) + 
    (-1.00 <= eta && eta < -0.90) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.176174) + 
    (-1.00 <= eta && eta < -0.90) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.198024) + 
    (-1.00 <= eta && eta < -0.90) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.219729) + 
    (-1.00 <= eta && eta < -0.90) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.240906) + 
    (-1.00 <= eta && eta < -0.90) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.262363) + 
    (-1.00 <= eta && eta < -0.90) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.283549) + 
    (-1.00 <= eta && eta < -0.90) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.304817) + 
    (-1.00 <= eta && eta < -0.90) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.326209) + 
    (-1.00 <= eta && eta < -0.90) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.347442) + 
    (-1.00 <= eta && eta < -0.90) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.368549) + 
    (-1.00 <= eta && eta < -0.90) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 30.70) * (0.389801) + 
    (-1.00 <= eta && eta < -0.90) * (30.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.411336) + 
    (-1.00 <= eta && eta < -0.90) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.90) * (0.433234) + 
    (-1.00 <= eta && eta < -0.90) * (40.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.453570) + 
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.000649) + 
    (-0.90 <= eta && eta < -0.80) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.020231) + 
    (-0.90 <= eta && eta < -0.80) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.041709) + 
    (-0.90 <= eta && eta < -0.80) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.063327) + 
    (-0.90 <= eta && eta < -0.80) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.085708) + 
    (-0.90 <= eta && eta < -0.80) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.107361) + 
    (-0.90 <= eta && eta < -0.80) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.129120) + 
    (-0.90 <= eta && eta < -0.80) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.152075) + 
    (-0.90 <= eta && eta < -0.80) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.174114) + 
    (-0.90 <= eta && eta < -0.80) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.196407) + 
    (-0.90 <= eta && eta < -0.80) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.218543) + 
    (-0.90 <= eta && eta < -0.80) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.240121) + 
    (-0.90 <= eta && eta < -0.80) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.261958) + 
    (-0.90 <= eta && eta < -0.80) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.283483) + 
    (-0.90 <= eta && eta < -0.80) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.305052) + 
    (-0.90 <= eta && eta < -0.80) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.326699) + 
    (-0.90 <= eta && eta < -0.80) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.348135) + 
    (-0.90 <= eta && eta < -0.80) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.369386) + 
    (-0.90 <= eta && eta < -0.80) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.390721) + 
    (-0.90 <= eta && eta < -0.80) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.412270) + 
    (-0.90 <= eta && eta < -0.80) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.434265) + 
    (-0.90 <= eta && eta < -0.80) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.454684) + 
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000698) + 
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.020156) + 
    (-0.80 <= eta && eta < -0.70) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.040719) + 
    (-0.80 <= eta && eta < -0.70) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.062503) + 
    (-0.80 <= eta && eta < -0.70) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.083481) + 
    (-0.80 <= eta && eta < -0.70) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.103705) + 
    (-0.80 <= eta && eta < -0.70) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.124110) + 
    (-0.80 <= eta && eta < -0.70) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.145883) + 
    (-0.80 <= eta && eta < -0.70) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.168517) + 
    (-0.80 <= eta && eta < -0.70) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.189997) + 
    (-0.80 <= eta && eta < -0.70) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.211530) + 
    (-0.80 <= eta && eta < -0.70) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.234017) + 
    (-0.80 <= eta && eta < -0.70) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.255655) + 
    (-0.80 <= eta && eta < -0.70) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.277170) + 
    (-0.80 <= eta && eta < -0.70) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.298995) + 
    (-0.80 <= eta && eta < -0.70) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.320333) + 
    (-0.80 <= eta && eta < -0.70) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.341878) + 
    (-0.80 <= eta && eta < -0.70) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.50) * (0.363611) + 
    (-0.80 <= eta && eta < -0.70) * (26.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.384994) + 
    (-0.80 <= eta && eta < -0.70) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.406381) + 
    (-0.80 <= eta && eta < -0.70) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.10) * (0.428249) + 
    (-0.80 <= eta && eta < -0.70) * (38.10 <= pt * cosh(eta) && pt * cosh(eta) < 47.90) * (0.450901) + 
    (-0.80 <= eta && eta < -0.70) * (47.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462545) + 
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000644) + 
    (-0.70 <= eta && eta < -0.60) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019416) + 
    (-0.70 <= eta && eta < -0.60) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.039800) + 
    (-0.70 <= eta && eta < -0.60) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.061565) + 
    (-0.70 <= eta && eta < -0.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.082602) + 
    (-0.70 <= eta && eta < -0.60) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.102923) + 
    (-0.70 <= eta && eta < -0.60) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.123448) + 
    (-0.70 <= eta && eta < -0.60) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.145363) + 
    (-0.70 <= eta && eta < -0.60) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.168148) + 
    (-0.70 <= eta && eta < -0.60) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.189770) + 
    (-0.70 <= eta && eta < -0.60) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.211440) + 
    (-0.70 <= eta && eta < -0.60) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.232802) + 
    (-0.70 <= eta && eta < -0.60) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.253505) + 
    (-0.70 <= eta && eta < -0.60) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.274359) + 
    (-0.70 <= eta && eta < -0.60) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.295736) + 
    (-0.70 <= eta && eta < -0.60) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.317652) + 
    (-0.70 <= eta && eta < -0.60) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.339116) + 
    (-0.70 <= eta && eta < -0.60) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.360427) + 
    (-0.70 <= eta && eta < -0.60) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.381782) + 
    (-0.70 <= eta && eta < -0.60) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.403198) + 
    (-0.70 <= eta && eta < -0.60) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.424964) + 
    (-0.70 <= eta && eta < -0.60) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.447490) + 
    (-0.70 <= eta && eta < -0.60) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461199) + 
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000623) + 
    (-0.60 <= eta && eta < -0.50) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019050) + 
    (-0.60 <= eta && eta < -0.50) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.039230) + 
    (-0.60 <= eta && eta < -0.50) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.060853) + 
    (-0.60 <= eta && eta < -0.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.081798) + 
    (-0.60 <= eta && eta < -0.50) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.102059) + 
    (-0.60 <= eta && eta < -0.50) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.122547) + 
    (-0.60 <= eta && eta < -0.50) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.144443) + 
    (-0.60 <= eta && eta < -0.50) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.167226) + 
    (-0.60 <= eta && eta < -0.50) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.188860) + 
    (-0.60 <= eta && eta < -0.50) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.210554) + 
    (-0.60 <= eta && eta < -0.50) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.231950) + 
    (-0.60 <= eta && eta < -0.50) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.252693) + 
    (-0.60 <= eta && eta < -0.50) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.273595) + 
    (-0.60 <= eta && eta < -0.50) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.295028) + 
    (-0.60 <= eta && eta < -0.50) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.317006) + 
    (-0.60 <= eta && eta < -0.50) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.338536) + 
    (-0.60 <= eta && eta < -0.50) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.359915) + 
    (-0.60 <= eta && eta < -0.50) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.381344) + 
    (-0.60 <= eta && eta < -0.50) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.402835) + 
    (-0.60 <= eta && eta < -0.50) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.424681) + 
    (-0.60 <= eta && eta < -0.50) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.60) * (0.447175) + 
    (-0.60 <= eta && eta < -0.50) * (45.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460966) + 
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000665) + 
    (-0.50 <= eta && eta < -0.40) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.020045) + 
    (-0.50 <= eta && eta < -0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040989) + 
    (-0.50 <= eta && eta < -0.40) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063254) + 
    (-0.50 <= eta && eta < -0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084688) + 
    (-0.50 <= eta && eta < -0.40) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105322) + 
    (-0.50 <= eta && eta < -0.40) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.126101) + 
    (-0.50 <= eta && eta < -0.40) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146555) + 
    (-0.50 <= eta && eta < -0.40) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167958) + 
    (-0.50 <= eta && eta < -0.40) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189862) + 
    (-0.50 <= eta && eta < -0.40) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211794) + 
    (-0.50 <= eta && eta < -0.40) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233391) + 
    (-0.50 <= eta && eta < -0.40) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254294) + 
    (-0.50 <= eta && eta < -0.40) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275318) + 
    (-0.50 <= eta && eta < -0.40) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296835) + 
    (-0.50 <= eta && eta < -0.40) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318854) + 
    (-0.50 <= eta && eta < -0.40) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340377) + 
    (-0.50 <= eta && eta < -0.40) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361703) + 
    (-0.50 <= eta && eta < -0.40) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.383028) + 
    (-0.50 <= eta && eta < -0.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404363) + 
    (-0.50 <= eta && eta < -0.40) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426192) + 
    (-0.50 <= eta && eta < -0.40) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448789) + 
    (-0.50 <= eta && eta < -0.40) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462263) + 
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000645) + 
    (-0.40 <= eta && eta < -0.30) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019306) + 
    (-0.40 <= eta && eta < -0.30) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.039824) + 
    (-0.40 <= eta && eta < -0.30) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.061801) + 
    (-0.40 <= eta && eta < -0.30) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.083051) + 
    (-0.40 <= eta && eta < -0.30) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.103567) + 
    (-0.40 <= eta && eta < -0.30) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.124275) + 
    (-0.40 <= eta && eta < -0.30) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.144695) + 
    (-0.40 <= eta && eta < -0.30) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.166094) + 
    (-0.40 <= eta && eta < -0.30) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.188023) + 
    (-0.40 <= eta && eta < -0.30) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.210005) + 
    (-0.40 <= eta && eta < -0.30) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.231670) + 
    (-0.40 <= eta && eta < -0.30) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.252656) + 
    (-0.40 <= eta && eta < -0.30) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.273778) + 
    (-0.40 <= eta && eta < -0.30) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.295407) + 
    (-0.40 <= eta && eta < -0.30) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.316782) + 
    (-0.40 <= eta && eta < -0.30) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.337879) + 
    (-0.40 <= eta && eta < -0.30) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.359062) + 
    (-0.40 <= eta && eta < -0.30) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.380434) + 
    (-0.40 <= eta && eta < -0.30) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.402053) + 
    (-0.40 <= eta && eta < -0.30) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.423693) + 
    (-0.40 <= eta && eta < -0.30) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.446017) + 
    (-0.40 <= eta && eta < -0.30) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460865) + 
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000624) + 
    (-0.30 <= eta && eta < -0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019209) + 
    (-0.30 <= eta && eta < -0.20) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.039736) + 
    (-0.30 <= eta && eta < -0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.061691) + 
    (-0.30 <= eta && eta < -0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.082927) + 
    (-0.30 <= eta && eta < -0.20) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.103434) + 
    (-0.30 <= eta && eta < -0.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.124136) + 
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.144554) + 
    (-0.30 <= eta && eta < -0.20) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.165953) + 
    (-0.30 <= eta && eta < -0.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.187883) + 
    (-0.30 <= eta && eta < -0.20) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.209869) + 
    (-0.30 <= eta && eta < -0.20) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.231539) + 
    (-0.30 <= eta && eta < -0.20) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.252531) + 
    (-0.30 <= eta && eta < -0.20) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.273661) + 
    (-0.30 <= eta && eta < -0.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.295298) + 
    (-0.30 <= eta && eta < -0.20) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.316682) + 
    (-0.30 <= eta && eta < -0.20) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.337790) + 
    (-0.30 <= eta && eta < -0.20) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.358983) + 
    (-0.30 <= eta && eta < -0.20) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.380366) + 
    (-0.30 <= eta && eta < -0.20) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.401997) + 
    (-0.30 <= eta && eta < -0.20) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.423648) + 
    (-0.30 <= eta && eta < -0.20) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.445986) + 
    (-0.30 <= eta && eta < -0.20) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.460842) + 
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.000667) + 
    (-0.20 <= eta && eta < -0.10) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.019658) + 
    (-0.20 <= eta && eta < -0.10) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.040176) + 
    (-0.20 <= eta && eta < -0.10) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.062033) + 
    (-0.20 <= eta && eta < -0.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.083131) + 
    (-0.20 <= eta && eta < -0.10) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.103490) + 
    (-0.20 <= eta && eta < -0.10) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.124039) + 
    (-0.20 <= eta && eta < -0.10) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.145966) + 
    (-0.20 <= eta && eta < -0.10) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.168751) + 
    (-0.20 <= eta && eta < -0.10) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.190365) + 
    (-0.20 <= eta && eta < -0.10) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.212019) + 
    (-0.20 <= eta && eta < -0.10) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.233359) + 
    (-0.20 <= eta && eta < -0.10) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.254035) + 
    (-0.20 <= eta && eta < -0.10) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.274858) + 
    (-0.20 <= eta && eta < -0.10) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.296199) + 
    (-0.20 <= eta && eta < -0.10) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.318075) + 
    (-0.20 <= eta && eta < -0.10) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.339496) + 
    (-0.20 <= eta && eta < -0.10) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.360760) + 
    (-0.20 <= eta && eta < -0.10) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.382069) + 
    (-0.20 <= eta && eta < -0.10) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.403435) + 
    (-0.20 <= eta && eta < -0.10) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.425149) + 
    (-0.20 <= eta && eta < -0.10) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.447621) + 
    (-0.20 <= eta && eta < -0.10) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461296) + 
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.000679) + 
    (-0.10 <= eta && eta < -0.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.019921) + 
    (-0.10 <= eta && eta < -0.00) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.039923) + 
    (-0.10 <= eta && eta < -0.00) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.061117) + 
    (-0.10 <= eta && eta < -0.00) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.083207) + 
    (-0.10 <= eta && eta < -0.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.104676) + 
    (-0.10 <= eta && eta < -0.00) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.126328) + 
    (-0.10 <= eta && eta < -0.00) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.149235) + 
    (-0.10 <= eta && eta < -0.00) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.171280) + 
    (-0.10 <= eta && eta < -0.00) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.192190) + 
    (-0.10 <= eta && eta < -0.00) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.213160) + 
    (-0.10 <= eta && eta < -0.00) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.235081) + 
    (-0.10 <= eta && eta < -0.00) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.257304) + 
    (-0.10 <= eta && eta < -0.00) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.279240) + 
    (-0.10 <= eta && eta < -0.00) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.300402) + 
    (-0.10 <= eta && eta < -0.00) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.321151) + 
    (-0.10 <= eta && eta < -0.00) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.342173) + 
    (-0.10 <= eta && eta < -0.00) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.363463) + 
    (-0.10 <= eta && eta < -0.00) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.384883) + 
    (-0.10 <= eta && eta < -0.00) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.406518) + 
    (-0.10 <= eta && eta < -0.00) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 38.90) * (0.428380) + 
    (-0.10 <= eta && eta < -0.00) * (38.90 <= pt * cosh(eta) && pt * cosh(eta) < 49.00) * (0.451065) + 
    (-0.10 <= eta && eta < -0.00) * (49.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461887) + 
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000652) + 
    (-0.00 <= eta && eta < 0.10) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.019853) + 
    (-0.00 <= eta && eta < 0.10) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041350) + 
    (-0.00 <= eta && eta < 0.10) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063083) + 
    (-0.00 <= eta && eta < 0.10) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.085621) + 
    (-0.00 <= eta && eta < 0.10) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.107436) + 
    (-0.00 <= eta && eta < 0.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.127685) + 
    (-0.00 <= eta && eta < 0.10) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.149221) + 
    (-0.00 <= eta && eta < 0.10) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.171557) + 
    (-0.00 <= eta && eta < 0.10) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.192724) + 
    (-0.00 <= eta && eta < 0.10) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.213928) + 
    (-0.00 <= eta && eta < 0.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.236064) + 
    (-0.00 <= eta && eta < 0.10) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.258471) + 
    (-0.00 <= eta && eta < 0.10) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.280550) + 
    (-0.00 <= eta && eta < 0.10) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.301813) + 
    (-0.00 <= eta && eta < 0.10) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.322623) + 
    (-0.00 <= eta && eta < 0.10) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.343665) + 
    (-0.00 <= eta && eta < 0.10) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.364932) + 
    (-0.00 <= eta && eta < 0.10) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.386281) + 
    (-0.00 <= eta && eta < 0.10) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.407794) + 
    (-0.00 <= eta && eta < 0.10) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 38.80) * (0.429658) + 
    (-0.00 <= eta && eta < 0.10) * (38.80 <= pt * cosh(eta) && pt * cosh(eta) < 49.20) * (0.452424) + 
    (-0.00 <= eta && eta < 0.10) * (49.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463042) + 
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000688) + 
    (0.10 <= eta && eta < 0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.020283) + 
    (0.10 <= eta && eta < 0.20) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.041390) + 
    (0.10 <= eta && eta < 0.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063752) + 
    (0.10 <= eta && eta < 0.20) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.085249) + 
    (0.10 <= eta && eta < 0.20) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105921) + 
    (0.10 <= eta && eta < 0.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.126724) + 
    (0.10 <= eta && eta < 0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.147189) + 
    (0.10 <= eta && eta < 0.20) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.168592) + 
    (0.10 <= eta && eta < 0.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.190487) + 
    (0.10 <= eta && eta < 0.20) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.212403) + 
    (0.10 <= eta && eta < 0.20) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233976) + 
    (0.10 <= eta && eta < 0.20) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254850) + 
    (0.10 <= eta && eta < 0.20) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275841) + 
    (0.10 <= eta && eta < 0.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.297319) + 
    (0.10 <= eta && eta < 0.20) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.319295) + 
    (0.10 <= eta && eta < 0.20) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340773) + 
    (0.10 <= eta && eta < 0.20) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.362051) + 
    (0.10 <= eta && eta < 0.20) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.383326) + 
    (0.10 <= eta && eta < 0.20) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404609) + 
    (0.10 <= eta && eta < 0.20) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426384) + 
    (0.10 <= eta && eta < 0.20) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.90) * (0.449033) + 
    (0.10 <= eta && eta < 0.20) * (45.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462443) + 
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000702) + 
    (0.20 <= eta && eta < 0.30) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.021415) + 
    (0.20 <= eta && eta < 0.30) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.041984) + 
    (0.20 <= eta && eta < 0.30) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.063242) + 
    (0.20 <= eta && eta < 0.30) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.085050) + 
    (0.20 <= eta && eta < 0.30) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.106045) + 
    (0.20 <= eta && eta < 0.30) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.127172) + 
    (0.20 <= eta && eta < 0.30) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.147941) + 
    (0.20 <= eta && eta < 0.30) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.169640) + 
    (0.20 <= eta && eta < 0.30) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.191806) + 
    (0.20 <= eta && eta < 0.30) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.213956) + 
    (0.20 <= eta && eta < 0.30) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.235721) + 
    (0.20 <= eta && eta < 0.30) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.256741) + 
    (0.20 <= eta && eta < 0.30) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.277837) + 
    (0.20 <= eta && eta < 0.30) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.299379) + 
    (0.20 <= eta && eta < 0.30) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.321372) + 
    (0.20 <= eta && eta < 0.30) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.342817) + 
    (0.20 <= eta && eta < 0.30) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.364014) + 
    (0.20 <= eta && eta < 0.30) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.385555) + 
    (0.20 <= eta && eta < 0.30) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.407162) + 
    (0.20 <= eta && eta < 0.30) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.429002) + 
    (0.20 <= eta && eta < 0.30) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.451798) + 
    (0.20 <= eta && eta < 0.30) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.464291) + 
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000726) + 
    (0.30 <= eta && eta < 0.40) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.020369) + 
    (0.30 <= eta && eta < 0.40) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.041656) + 
    (0.30 <= eta && eta < 0.40) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.064296) + 
    (0.30 <= eta && eta < 0.40) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.086047) + 
    (0.30 <= eta && eta < 0.40) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.106941) + 
    (0.30 <= eta && eta < 0.40) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.127942) + 
    (0.30 <= eta && eta < 0.40) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.148575) + 
    (0.30 <= eta && eta < 0.40) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.170125) + 
    (0.30 <= eta && eta < 0.40) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.192139) + 
    (0.30 <= eta && eta < 0.40) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.214143) + 
    (0.30 <= eta && eta < 0.40) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.235774) + 
    (0.30 <= eta && eta < 0.40) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.256674) + 
    (0.30 <= eta && eta < 0.40) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.277665) + 
    (0.30 <= eta && eta < 0.40) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.299113) + 
    (0.30 <= eta && eta < 0.40) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.321029) + 
    (0.30 <= eta && eta < 0.40) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.342418) + 
    (0.30 <= eta && eta < 0.40) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.363580) + 
    (0.30 <= eta && eta < 0.40) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.385108) + 
    (0.30 <= eta && eta < 0.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.406727) + 
    (0.30 <= eta && eta < 0.40) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.428608) + 
    (0.30 <= eta && eta < 0.40) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.451376) + 
    (0.30 <= eta && eta < 0.40) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463881) + 
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.60) * (0.000671) + 
    (0.40 <= eta && eta < 0.50) * (9.60 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.020318) + 
    (0.40 <= eta && eta < 0.50) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.040535) + 
    (0.40 <= eta && eta < 0.50) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.061663) + 
    (0.40 <= eta && eta < 0.50) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.083458) + 
    (0.40 <= eta && eta < 0.50) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.104511) + 
    (0.40 <= eta && eta < 0.50) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.125739) + 
    (0.40 <= eta && eta < 0.50) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.146635) + 
    (0.40 <= eta && eta < 0.50) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.168484) + 
    (0.40 <= eta && eta < 0.50) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.190813) + 
    (0.40 <= eta && eta < 0.50) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.213130) + 
    (0.40 <= eta && eta < 0.50) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.235057) + 
    (0.40 <= eta && eta < 0.50) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.256227) + 
    (0.40 <= eta && eta < 0.50) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.277468) + 
    (0.40 <= eta && eta < 0.50) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.299144) + 
    (0.40 <= eta && eta < 0.50) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.321260) + 
    (0.40 <= eta && eta < 0.50) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.342810) + 
    (0.40 <= eta && eta < 0.50) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.364092) + 
    (0.40 <= eta && eta < 0.50) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.385699) + 
    (0.40 <= eta && eta < 0.50) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.407349) + 
    (0.40 <= eta && eta < 0.50) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.429205) + 
    (0.40 <= eta && eta < 0.50) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.70) * (0.451986) + 
    (0.40 <= eta && eta < 0.50) * (46.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.464527) + 
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.000629) + 
    (0.50 <= eta && eta < 0.60) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.019027) + 
    (0.50 <= eta && eta < 0.60) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.039568) + 
    (0.50 <= eta && eta < 0.60) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.061693) + 
    (0.50 <= eta && eta < 0.60) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.083113) + 
    (0.50 <= eta && eta < 0.60) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.103800) + 
    (0.50 <= eta && eta < 0.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.124675) + 
    (0.50 <= eta && eta < 0.60) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.145250) + 
    (0.50 <= eta && eta < 0.60) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.166798) + 
    (0.50 <= eta && eta < 0.60) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.188860) + 
    (0.50 <= eta && eta < 0.60) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.210955) + 
    (0.50 <= eta && eta < 0.60) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.232711) + 
    (0.50 <= eta && eta < 0.60) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.253762) + 
    (0.50 <= eta && eta < 0.60) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.274929) + 
    (0.50 <= eta && eta < 0.60) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.296580) + 
    (0.50 <= eta && eta < 0.60) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.317952) + 
    (0.50 <= eta && eta < 0.60) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.339024) + 
    (0.50 <= eta && eta < 0.60) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.360155) + 
    (0.50 <= eta && eta < 0.60) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.381450) + 
    (0.50 <= eta && eta < 0.60) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.402963) + 
    (0.50 <= eta && eta < 0.60) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.424674) + 
    (0.50 <= eta && eta < 0.60) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.70) * (0.447121) + 
    (0.50 <= eta && eta < 0.60) * (44.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.461667) + 
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000675) + 
    (0.60 <= eta && eta < 0.70) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019837) + 
    (0.60 <= eta && eta < 0.70) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040690) + 
    (0.60 <= eta && eta < 0.70) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.062882) + 
    (0.60 <= eta && eta < 0.70) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084270) + 
    (0.60 <= eta && eta < 0.70) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.104874) + 
    (0.60 <= eta && eta < 0.70) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.125635) + 
    (0.60 <= eta && eta < 0.70) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146081) + 
    (0.60 <= eta && eta < 0.70) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167483) + 
    (0.60 <= eta && eta < 0.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189393) + 
    (0.60 <= eta && eta < 0.70) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211339) + 
    (0.60 <= eta && eta < 0.70) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.232954) + 
    (0.60 <= eta && eta < 0.70) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.253877) + 
    (0.60 <= eta && eta < 0.70) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.274927) + 
    (0.60 <= eta && eta < 0.70) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296472) + 
    (0.60 <= eta && eta < 0.70) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318523) + 
    (0.60 <= eta && eta < 0.70) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340081) + 
    (0.60 <= eta && eta < 0.70) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361442) + 
    (0.60 <= eta && eta < 0.70) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.382805) + 
    (0.60 <= eta && eta < 0.70) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404178) + 
    (0.60 <= eta && eta < 0.70) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426049) + 
    (0.60 <= eta && eta < 0.70) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448688) + 
    (0.60 <= eta && eta < 0.70) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462189) + 
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.000676) + 
    (0.70 <= eta && eta < 0.80) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.019969) + 
    (0.70 <= eta && eta < 0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.040896) + 
    (0.70 <= eta && eta < 0.80) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.063138) + 
    (0.70 <= eta && eta < 0.80) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.084558) + 
    (0.70 <= eta && eta < 0.80) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.105182) + 
    (0.70 <= eta && eta < 0.80) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.125956) + 
    (0.70 <= eta && eta < 0.80) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.146407) + 
    (0.70 <= eta && eta < 0.80) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.167810) + 
    (0.70 <= eta && eta < 0.80) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.189716) + 
    (0.70 <= eta && eta < 0.80) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.211653) + 
    (0.70 <= eta && eta < 0.80) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.233255) + 
    (0.70 <= eta && eta < 0.80) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.254164) + 
    (0.70 <= eta && eta < 0.80) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.275196) + 
    (0.70 <= eta && eta < 0.80) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.296722) + 
    (0.70 <= eta && eta < 0.80) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.318751) + 
    (0.70 <= eta && eta < 0.80) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.340285) + 
    (0.70 <= eta && eta < 0.80) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.361622) + 
    (0.70 <= eta && eta < 0.80) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.382958) + 
    (0.70 <= eta && eta < 0.80) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.404306) + 
    (0.70 <= eta && eta < 0.80) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.426148) + 
    (0.70 <= eta && eta < 0.80) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.448757) + 
    (0.70 <= eta && eta < 0.80) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462240) + 
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000666) + 
    (0.80 <= eta && eta < 0.90) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.019973) + 
    (0.80 <= eta && eta < 0.90) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041538) + 
    (0.80 <= eta && eta < 0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063316) + 
    (0.80 <= eta && eta < 0.90) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.085883) + 
    (0.80 <= eta && eta < 0.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.107717) + 
    (0.80 <= eta && eta < 0.90) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.127976) + 
    (0.80 <= eta && eta < 0.90) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.149517) + 
    (0.80 <= eta && eta < 0.90) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.171852) + 
    (0.80 <= eta && eta < 0.90) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.193015) + 
    (0.80 <= eta && eta < 0.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.214211) + 
    (0.80 <= eta && eta < 0.90) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.236336) + 
    (0.80 <= eta && eta < 0.90) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.258728) + 
    (0.80 <= eta && eta < 0.90) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.280790) + 
    (0.80 <= eta && eta < 0.90) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.302035) + 
    (0.80 <= eta && eta < 0.90) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.322825) + 
    (0.80 <= eta && eta < 0.90) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.343847) + 
    (0.80 <= eta && eta < 0.90) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.365091) + 
    (0.80 <= eta && eta < 0.90) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.386417) + 
    (0.80 <= eta && eta < 0.90) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.407905) + 
    (0.80 <= eta && eta < 0.90) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 38.80) * (0.429743) + 
    (0.80 <= eta && eta < 0.90) * (38.80 <= pt * cosh(eta) && pt * cosh(eta) < 49.20) * (0.452482) + 
    (0.80 <= eta && eta < 0.90) * (49.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.463087) + 
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.000631) + 
    (0.90 <= eta && eta < 1.00) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.020111) + 
    (0.90 <= eta && eta < 1.00) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.041754) + 
    (0.90 <= eta && eta < 1.00) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.063581) + 
    (0.90 <= eta && eta < 1.00) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.084516) + 
    (0.90 <= eta && eta < 1.00) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.104653) + 
    (0.90 <= eta && eta < 1.00) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.124945) + 
    (0.90 <= eta && eta < 1.00) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.146583) + 
    (0.90 <= eta && eta < 1.00) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.169067) + 
    (0.90 <= eta && eta < 1.00) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.190404) + 
    (0.90 <= eta && eta < 1.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.211798) + 
    (0.90 <= eta && eta < 1.00) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.234147) + 
    (0.90 <= eta && eta < 1.00) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.255662) + 
    (0.90 <= eta && eta < 1.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.277068) + 
    (0.90 <= eta && eta < 1.00) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.298797) + 
    (0.90 <= eta && eta < 1.00) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.320057) + 
    (0.90 <= eta && eta < 1.00) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.341541) + 
    (0.90 <= eta && eta < 1.00) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.363233) + 
    (0.90 <= eta && eta < 1.00) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.384597) + 
    (0.90 <= eta && eta < 1.00) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.80) * (0.405989) + 
    (0.90 <= eta && eta < 1.00) * (32.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.427888) + 
    (0.90 <= eta && eta < 1.00) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 48.00) * (0.450605) + 
    (0.90 <= eta && eta < 1.00) * (48.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.462220) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 

  add EfficiencyFormula {2212} {2212} {
    (eta<-1.00 || eta>=1.00 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (-1.00 <= eta && eta < -0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-1.00 <= eta && eta < -0.90) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.50) * (0.999393) + 
    (-1.00 <= eta && eta < -0.90) * (9.50 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.979966) + 
    (-1.00 <= eta && eta < -0.90) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.959319) + 
    (-1.00 <= eta && eta < -0.90) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.938208) + 
    (-1.00 <= eta && eta < -0.90) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.916518) + 
    (-1.00 <= eta && eta < -0.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.893591) + 
    (-1.00 <= eta && eta < -0.90) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.870367) + 
    (-1.00 <= eta && eta < -0.90) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.847581) + 
    (-1.00 <= eta && eta < -0.90) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.825759) + 
    (-1.00 <= eta && eta < -0.90) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.803286) + 
    (-1.00 <= eta && eta < -0.90) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.780814) + 
    (-1.00 <= eta && eta < -0.90) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.759031) + 
    (-1.00 <= eta && eta < -0.90) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.737077) + 
    (-1.00 <= eta && eta < -0.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.714871) + 
    (-1.00 <= eta && eta < -0.90) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.692935) + 
    (-1.00 <= eta && eta < -0.90) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.671054) + 
    (-1.00 <= eta && eta < -0.90) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.648575) + 
    (-1.00 <= eta && eta < -0.90) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.00) * (0.621769) + 
    (-1.00 <= eta && eta < -0.90) * (32.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.605753) + 
    (-0.90 <= eta && eta < -0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-0.90 <= eta && eta < -0.80) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.999362) + 
    (-0.90 <= eta && eta < -0.80) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.978828) + 
    (-0.90 <= eta && eta < -0.80) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.957392) + 
    (-0.90 <= eta && eta < -0.80) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.935603) + 
    (-0.90 <= eta && eta < -0.80) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.913368) + 
    (-0.90 <= eta && eta < -0.80) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.890006) + 
    (-0.90 <= eta && eta < -0.80) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.866467) + 
    (-0.90 <= eta && eta < -0.80) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.843486) + 
    (-0.90 <= eta && eta < -0.80) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.821573) + 
    (-0.90 <= eta && eta < -0.80) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.799102) + 
    (-0.90 <= eta && eta < -0.80) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.776721) + 
    (-0.90 <= eta && eta < -0.80) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.755112) + 
    (-0.90 <= eta && eta < -0.80) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.733415) + 
    (-0.90 <= eta && eta < -0.80) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.711555) + 
    (-0.90 <= eta && eta < -0.80) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.689217) + 
    (-0.90 <= eta && eta < -0.80) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.666846) + 
    (-0.90 <= eta && eta < -0.80) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.643904) + 
    (-0.90 <= eta && eta < -0.80) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609976) + 
    (-0.80 <= eta && eta < -0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-0.80 <= eta && eta < -0.70) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.999417) + 
    (-0.80 <= eta && eta < -0.70) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.980438) + 
    (-0.80 <= eta && eta < -0.70) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.959304) + 
    (-0.80 <= eta && eta < -0.70) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.937529) + 
    (-0.80 <= eta && eta < -0.70) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.915124) + 
    (-0.80 <= eta && eta < -0.70) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.891475) + 
    (-0.80 <= eta && eta < -0.70) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.867589) + 
    (-0.80 <= eta && eta < -0.70) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.844246) + 
    (-0.80 <= eta && eta < -0.70) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.821987) + 
    (-0.80 <= eta && eta < -0.70) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.799175) + 
    (-0.80 <= eta && eta < -0.70) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.776484) + 
    (-0.80 <= eta && eta < -0.70) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.754614) + 
    (-0.80 <= eta && eta < -0.70) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.732703) + 
    (-0.80 <= eta && eta < -0.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.710687) + 
    (-0.80 <= eta && eta < -0.70) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.688267) + 
    (-0.80 <= eta && eta < -0.70) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.665905) + 
    (-0.80 <= eta && eta < -0.70) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.643095) + 
    (-0.80 <= eta && eta < -0.70) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609665) + 
    (-0.70 <= eta && eta < -0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-0.70 <= eta && eta < -0.60) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.999371) + 
    (-0.70 <= eta && eta < -0.60) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.979255) + 
    (-0.70 <= eta && eta < -0.60) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.957127) + 
    (-0.70 <= eta && eta < -0.60) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.934823) + 
    (-0.70 <= eta && eta < -0.60) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.912064) + 
    (-0.70 <= eta && eta < -0.60) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.888193) + 
    (-0.70 <= eta && eta < -0.60) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.864204) + 
    (-0.70 <= eta && eta < -0.60) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.840855) + 
    (-0.70 <= eta && eta < -0.60) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.818667) + 
    (-0.70 <= eta && eta < -0.60) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.795994) + 
    (-0.70 <= eta && eta < -0.60) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.773499) + 
    (-0.70 <= eta && eta < -0.60) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.751868) + 
    (-0.70 <= eta && eta < -0.60) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.730241) + 
    (-0.70 <= eta && eta < -0.60) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.708551) + 
    (-0.70 <= eta && eta < -0.60) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.686502) + 
    (-0.70 <= eta && eta < -0.60) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.664546) + 
    (-0.70 <= eta && eta < -0.60) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.641866) + 
    (-0.70 <= eta && eta < -0.60) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609387) + 
    (-0.60 <= eta && eta < -0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-0.60 <= eta && eta < -0.50) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.999416) + 
    (-0.60 <= eta && eta < -0.50) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.979498) + 
    (-0.60 <= eta && eta < -0.50) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.957852) + 
    (-0.60 <= eta && eta < -0.50) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.935720) + 
    (-0.60 <= eta && eta < -0.50) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.913076) + 
    (-0.60 <= eta && eta < -0.50) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.889276) + 
    (-0.60 <= eta && eta < -0.50) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.865320) + 
    (-0.60 <= eta && eta < -0.50) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.841972) + 
    (-0.60 <= eta && eta < -0.50) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.819759) + 
    (-0.60 <= eta && eta < -0.50) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.797040) + 
    (-0.60 <= eta && eta < -0.50) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.774479) + 
    (-0.60 <= eta && eta < -0.50) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.752769) + 
    (-0.60 <= eta && eta < -0.50) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.731048) + 
    (-0.60 <= eta && eta < -0.50) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.709251) + 
    (-0.60 <= eta && eta < -0.50) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.687080) + 
    (-0.60 <= eta && eta < -0.50) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.664991) + 
    (-0.60 <= eta && eta < -0.50) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.642162) + 
    (-0.60 <= eta && eta < -0.50) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609449) + 
    (-0.50 <= eta && eta < -0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-0.50 <= eta && eta < -0.40) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.999377) + 
    (-0.50 <= eta && eta < -0.40) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.978482) + 
    (-0.50 <= eta && eta < -0.40) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.955953) + 
    (-0.50 <= eta && eta < -0.40) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.933137) + 
    (-0.50 <= eta && eta < -0.40) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.909949) + 
    (-0.50 <= eta && eta < -0.40) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.885719) + 
    (-0.50 <= eta && eta < -0.40) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.861459) + 
    (-0.50 <= eta && eta < -0.40) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.837929) + 
    (-0.50 <= eta && eta < -0.40) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.815641) + 
    (-0.50 <= eta && eta < -0.40) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.792940) + 
    (-0.50 <= eta && eta < -0.40) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.770487) + 
    (-0.50 <= eta && eta < -0.40) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.748964) + 
    (-0.50 <= eta && eta < -0.40) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.727511) + 
    (-0.50 <= eta && eta < -0.40) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.706066) + 
    (-0.50 <= eta && eta < -0.40) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.684338) + 
    (-0.50 <= eta && eta < -0.40) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.662237) + 
    (-0.50 <= eta && eta < -0.40) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.639038) + 
    (-0.50 <= eta && eta < -0.40) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608741) + 
    (-0.40 <= eta && eta < -0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-0.40 <= eta && eta < -0.30) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.999391) + 
    (-0.40 <= eta && eta < -0.30) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.979179) + 
    (-0.40 <= eta && eta < -0.30) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.957431) + 
    (-0.40 <= eta && eta < -0.30) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.934965) + 
    (-0.40 <= eta && eta < -0.30) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.912006) + 
    (-0.40 <= eta && eta < -0.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.887916) + 
    (-0.40 <= eta && eta < -0.30) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.863717) + 
    (-0.40 <= eta && eta < -0.30) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.840182) + 
    (-0.40 <= eta && eta < -0.30) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.817839) + 
    (-0.40 <= eta && eta < -0.30) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.795040) + 
    (-0.40 <= eta && eta < -0.30) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.772450) + 
    (-0.40 <= eta && eta < -0.30) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.750764) + 
    (-0.40 <= eta && eta < -0.30) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.729121) + 
    (-0.40 <= eta && eta < -0.30) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.707458) + 
    (-0.40 <= eta && eta < -0.30) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.685485) + 
    (-0.40 <= eta && eta < -0.30) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.663109) + 
    (-0.40 <= eta && eta < -0.30) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.639895) + 
    (-0.40 <= eta && eta < -0.30) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608937) + 
    (-0.30 <= eta && eta < -0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-0.30 <= eta && eta < -0.20) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.999403) + 
    (-0.30 <= eta && eta < -0.20) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.979783) + 
    (-0.30 <= eta && eta < -0.20) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.957542) + 
    (-0.30 <= eta && eta < -0.20) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.935103) + 
    (-0.30 <= eta && eta < -0.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.912162) + 
    (-0.30 <= eta && eta < -0.20) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.888083) + 
    (-0.30 <= eta && eta < -0.20) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.863888) + 
    (-0.30 <= eta && eta < -0.20) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.840353) + 
    (-0.30 <= eta && eta < -0.20) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.818007) + 
    (-0.30 <= eta && eta < -0.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.795199) + 
    (-0.30 <= eta && eta < -0.20) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.772600) + 
    (-0.30 <= eta && eta < -0.20) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.750902) + 
    (-0.30 <= eta && eta < -0.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.729244) + 
    (-0.30 <= eta && eta < -0.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.707564) + 
    (-0.30 <= eta && eta < -0.20) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.685572) + 
    (-0.30 <= eta && eta < -0.20) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.663176) + 
    (-0.30 <= eta && eta < -0.20) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.639938) + 
    (-0.30 <= eta && eta < -0.20) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608946) + 
    (-0.20 <= eta && eta < -0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-0.20 <= eta && eta < -0.10) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.999359) + 
    (-0.20 <= eta && eta < -0.10) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.00) * (0.978520) + 
    (-0.20 <= eta && eta < -0.10) * (10.00 <= pt * cosh(eta) && pt * cosh(eta) < 10.60) * (0.956649) + 
    (-0.20 <= eta && eta < -0.10) * (10.60 <= pt * cosh(eta) && pt * cosh(eta) < 11.10) * (0.934232) + 
    (-0.20 <= eta && eta < -0.10) * (11.10 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.911399) + 
    (-0.20 <= eta && eta < -0.10) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.887482) + 
    (-0.20 <= eta && eta < -0.10) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.863473) + 
    (-0.20 <= eta && eta < -0.10) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.840125) + 
    (-0.20 <= eta && eta < -0.10) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.817953) + 
    (-0.20 <= eta && eta < -0.10) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.795311) + 
    (-0.20 <= eta && eta < -0.10) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.772859) + 
    (-0.20 <= eta && eta < -0.10) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.751280) + 
    (-0.20 <= eta && eta < -0.10) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.729714) + 
    (-0.20 <= eta && eta < -0.10) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.708095) + 
    (-0.20 <= eta && eta < -0.10) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.686125) + 
    (-0.20 <= eta && eta < -0.10) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.663701) + 
    (-0.20 <= eta && eta < -0.10) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.640366) + 
    (-0.20 <= eta && eta < -0.10) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609100) + 
    (-0.10 <= eta && eta < -0.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-0.10 <= eta && eta < -0.00) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.40) * (0.999439) + 
    (-0.10 <= eta && eta < -0.00) * (9.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.980491) + 
    (-0.10 <= eta && eta < -0.00) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.959593) + 
    (-0.10 <= eta && eta < -0.00) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.938333) + 
    (-0.10 <= eta && eta < -0.00) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.916456) + 
    (-0.10 <= eta && eta < -0.00) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.893323) + 
    (-0.10 <= eta && eta < -0.00) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.869897) + 
    (-0.10 <= eta && eta < -0.00) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.846930) + 
    (-0.10 <= eta && eta < -0.00) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.824955) + 
    (-0.10 <= eta && eta < -0.00) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.802353) + 
    (-0.10 <= eta && eta < -0.00) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.779781) + 
    (-0.10 <= eta && eta < -0.00) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.757937) + 
    (-0.10 <= eta && eta < -0.00) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.735958) + 
    (-0.10 <= eta && eta < -0.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.713770) + 
    (-0.10 <= eta && eta < -0.00) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.691898) + 
    (-0.10 <= eta && eta < -0.00) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.670134) + 
    (-0.10 <= eta && eta < -0.00) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.647473) + 
    (-0.10 <= eta && eta < -0.00) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610659) + 
    (-0.00 <= eta && eta < 0.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (-0.00 <= eta && eta < 0.10) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.999354) + 
    (-0.00 <= eta && eta < 0.10) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.979103) + 
    (-0.00 <= eta && eta < 0.10) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.958182) + 
    (-0.00 <= eta && eta < 0.10) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.936357) + 
    (-0.00 <= eta && eta < 0.10) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.914009) + 
    (-0.00 <= eta && eta < 0.10) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.890487) + 
    (-0.00 <= eta && eta < 0.10) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.866768) + 
    (-0.00 <= eta && eta < 0.10) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.843606) + 
    (-0.00 <= eta && eta < 0.10) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.821527) + 
    (-0.00 <= eta && eta < 0.10) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.798897) + 
    (-0.00 <= eta && eta < 0.10) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.776376) + 
    (-0.00 <= eta && eta < 0.10) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.754655) + 
    (-0.00 <= eta && eta < 0.10) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.732873) + 
    (-0.00 <= eta && eta < 0.10) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.710960) + 
    (-0.00 <= eta && eta < 0.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.688609) + 
    (-0.00 <= eta && eta < 0.10) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.666272) + 
    (-0.00 <= eta && eta < 0.10) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.643430) + 
    (-0.00 <= eta && eta < 0.10) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609805) + 
    (0.10 <= eta && eta < 0.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (0.10 <= eta && eta < 0.20) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.999442) + 
    (0.10 <= eta && eta < 0.20) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.980896) + 
    (0.10 <= eta && eta < 0.20) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.959254) + 
    (0.10 <= eta && eta < 0.20) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.936937) + 
    (0.10 <= eta && eta < 0.20) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.914015) + 
    (0.10 <= eta && eta < 0.20) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.889852) + 
    (0.10 <= eta && eta < 0.20) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.865505) + 
    (0.10 <= eta && eta < 0.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.841780) + 
    (0.10 <= eta && eta < 0.20) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.819231) + 
    (0.10 <= eta && eta < 0.20) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.796203) + 
    (0.10 <= eta && eta < 0.20) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.773383) + 
    (0.10 <= eta && eta < 0.20) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.751478) + 
    (0.10 <= eta && eta < 0.20) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.729625) + 
    (0.10 <= eta && eta < 0.20) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.707770) + 
    (0.10 <= eta && eta < 0.20) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.685626) + 
    (0.10 <= eta && eta < 0.20) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.663668) + 
    (0.10 <= eta && eta < 0.20) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.640799) + 
    (0.10 <= eta && eta < 0.20) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609016) + 
    (0.20 <= eta && eta < 0.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (0.20 <= eta && eta < 0.30) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.999456) + 
    (0.20 <= eta && eta < 0.30) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.980442) + 
    (0.20 <= eta && eta < 0.30) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.959342) + 
    (0.20 <= eta && eta < 0.30) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.938527) + 
    (0.20 <= eta && eta < 0.30) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.915372) + 
    (0.20 <= eta && eta < 0.30) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.890867) + 
    (0.20 <= eta && eta < 0.30) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.866131) + 
    (0.20 <= eta && eta < 0.30) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.842017) + 
    (0.20 <= eta && eta < 0.30) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.819110) + 
    (0.20 <= eta && eta < 0.30) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.797763) + 
    (0.20 <= eta && eta < 0.30) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.776311) + 
    (0.20 <= eta && eta < 0.30) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.753723) + 
    (0.20 <= eta && eta < 0.30) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.731218) + 
    (0.20 <= eta && eta < 0.30) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.708760) + 
    (0.20 <= eta && eta < 0.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.686083) + 
    (0.20 <= eta && eta < 0.30) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.663704) + 
    (0.20 <= eta && eta < 0.30) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.640869) + 
    (0.20 <= eta && eta < 0.30) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608858) + 
    (0.30 <= eta && eta < 0.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (0.30 <= eta && eta < 0.40) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.999408) + 
    (0.30 <= eta && eta < 0.40) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.980163) + 
    (0.30 <= eta && eta < 0.40) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.30) * (0.959038) + 
    (0.30 <= eta && eta < 0.40) * (10.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.80) * (0.936756) + 
    (0.30 <= eta && eta < 0.40) * (10.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.913584) + 
    (0.30 <= eta && eta < 0.40) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.889168) + 
    (0.30 <= eta && eta < 0.40) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.864590) + 
    (0.30 <= eta && eta < 0.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.840671) + 
    (0.30 <= eta && eta < 0.40) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.817972) + 
    (0.30 <= eta && eta < 0.40) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.794832) + 
    (0.30 <= eta && eta < 0.40) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.771941) + 
    (0.30 <= eta && eta < 0.40) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.750012) + 
    (0.30 <= eta && eta < 0.40) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.728180) + 
    (0.30 <= eta && eta < 0.40) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.706396) + 
    (0.30 <= eta && eta < 0.40) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.684379) + 
    (0.30 <= eta && eta < 0.40) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.662061) + 
    (0.30 <= eta && eta < 0.40) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.638751) + 
    (0.30 <= eta && eta < 0.40) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608560) + 
    (0.40 <= eta && eta < 0.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (0.40 <= eta && eta < 0.50) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 8.90) * (0.999373) + 
    (0.40 <= eta && eta < 0.50) * (8.90 <= pt * cosh(eta) && pt * cosh(eta) < 9.70) * (0.978769) + 
    (0.40 <= eta && eta < 0.50) * (9.70 <= pt * cosh(eta) && pt * cosh(eta) < 10.20) * (0.957448) + 
    (0.40 <= eta && eta < 0.50) * (10.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.936499) + 
    (0.40 <= eta && eta < 0.50) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.913062) + 
    (0.40 <= eta && eta < 0.50) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.888382) + 
    (0.40 <= eta && eta < 0.50) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.863567) + 
    (0.40 <= eta && eta < 0.50) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.839454) + 
    (0.40 <= eta && eta < 0.50) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.816607) + 
    (0.40 <= eta && eta < 0.50) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.795361) + 
    (0.40 <= eta && eta < 0.50) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.774050) + 
    (0.40 <= eta && eta < 0.50) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.751650) + 
    (0.40 <= eta && eta < 0.50) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.729367) + 
    (0.40 <= eta && eta < 0.50) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.707163) + 
    (0.40 <= eta && eta < 0.50) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.684773) + 
    (0.40 <= eta && eta < 0.50) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.662704) + 
    (0.40 <= eta && eta < 0.50) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.639904) + 
    (0.40 <= eta && eta < 0.50) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608646) + 
    (0.50 <= eta && eta < 0.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (0.50 <= eta && eta < 0.60) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.999417) + 
    (0.50 <= eta && eta < 0.60) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.979539) + 
    (0.50 <= eta && eta < 0.60) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.957986) + 
    (0.50 <= eta && eta < 0.60) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.935548) + 
    (0.50 <= eta && eta < 0.60) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.912442) + 
    (0.50 <= eta && eta < 0.60) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.888165) + 
    (0.50 <= eta && eta < 0.60) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.863766) + 
    (0.50 <= eta && eta < 0.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.840041) + 
    (0.50 <= eta && eta < 0.60) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.817531) + 
    (0.50 <= eta && eta < 0.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.794579) + 
    (0.50 <= eta && eta < 0.60) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.771863) + 
    (0.50 <= eta && eta < 0.60) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.750084) + 
    (0.50 <= eta && eta < 0.60) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.728379) + 
    (0.50 <= eta && eta < 0.60) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.706692) + 
    (0.50 <= eta && eta < 0.60) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.684739) + 
    (0.50 <= eta && eta < 0.60) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.30) * (0.662439) + 
    (0.50 <= eta && eta < 0.60) * (20.30 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.639080) + 
    (0.50 <= eta && eta < 0.60) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608689) + 
    (0.60 <= eta && eta < 0.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (0.60 <= eta && eta < 0.70) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.10) * (0.999335) + 
    (0.60 <= eta && eta < 0.70) * (9.10 <= pt * cosh(eta) && pt * cosh(eta) < 9.90) * (0.978710) + 
    (0.60 <= eta && eta < 0.70) * (9.90 <= pt * cosh(eta) && pt * cosh(eta) < 10.50) * (0.956332) + 
    (0.60 <= eta && eta < 0.70) * (10.50 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.933605) + 
    (0.60 <= eta && eta < 0.70) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.910475) + 
    (0.60 <= eta && eta < 0.70) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.886279) + 
    (0.60 <= eta && eta < 0.70) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.862034) + 
    (0.60 <= eta && eta < 0.70) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.838503) + 
    (0.60 <= eta && eta < 0.70) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.816201) + 
    (0.60 <= eta && eta < 0.70) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.793474) + 
    (0.60 <= eta && eta < 0.70) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.770986) + 
    (0.60 <= eta && eta < 0.70) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.749421) + 
    (0.60 <= eta && eta < 0.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.727920) + 
    (0.60 <= eta && eta < 0.70) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.706419) + 
    (0.60 <= eta && eta < 0.70) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.684629) + 
    (0.60 <= eta && eta < 0.70) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.662458) + 
    (0.60 <= eta && eta < 0.70) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.639180) + 
    (0.60 <= eta && eta < 0.70) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608771) + 
    (0.70 <= eta && eta < 0.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (0.70 <= eta && eta < 0.80) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.00) * (0.999463) + 
    (0.70 <= eta && eta < 0.80) * (9.00 <= pt * cosh(eta) && pt * cosh(eta) < 9.80) * (0.981050) + 
    (0.70 <= eta && eta < 0.80) * (9.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.40) * (0.959830) + 
    (0.70 <= eta && eta < 0.80) * (10.40 <= pt * cosh(eta) && pt * cosh(eta) < 10.90) * (0.937687) + 
    (0.70 <= eta && eta < 0.80) * (10.90 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.914867) + 
    (0.70 <= eta && eta < 0.80) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.890769) + 
    (0.70 <= eta && eta < 0.80) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.866452) + 
    (0.70 <= eta && eta < 0.80) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.842728) + 
    (0.70 <= eta && eta < 0.80) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.820157) + 
    (0.70 <= eta && eta < 0.80) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.797090) + 
    (0.70 <= eta && eta < 0.80) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.774214) + 
    (0.70 <= eta && eta < 0.80) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.752241) + 
    (0.70 <= eta && eta < 0.80) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.730308) + 
    (0.70 <= eta && eta < 0.80) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.708360) + 
    (0.70 <= eta && eta < 0.80) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.686112) + 
    (0.70 <= eta && eta < 0.80) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.664041) + 
    (0.70 <= eta && eta < 0.80) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.40) * (0.641356) + 
    (0.70 <= eta && eta < 0.80) * (23.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609149) + 
    (0.80 <= eta && eta < 0.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (0.80 <= eta && eta < 0.90) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.20) * (0.999454) + 
    (0.80 <= eta && eta < 0.90) * (9.20 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.979938) + 
    (0.80 <= eta && eta < 0.90) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.957951) + 
    (0.80 <= eta && eta < 0.90) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.936071) + 
    (0.80 <= eta && eta < 0.90) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.913686) + 
    (0.80 <= eta && eta < 0.90) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.890140) + 
    (0.80 <= eta && eta < 0.90) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.866410) + 
    (0.80 <= eta && eta < 0.90) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.843247) + 
    (0.80 <= eta && eta < 0.90) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.821175) + 
    (0.80 <= eta && eta < 0.90) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.798559) + 
    (0.80 <= eta && eta < 0.90) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.776059) + 
    (0.80 <= eta && eta < 0.90) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.754363) + 
    (0.80 <= eta && eta < 0.90) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.732610) + 
    (0.80 <= eta && eta < 0.90) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.710732) + 
    (0.80 <= eta && eta < 0.90) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.688420) + 
    (0.80 <= eta && eta < 0.90) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.666126) + 
    (0.80 <= eta && eta < 0.90) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.643332) + 
    (0.80 <= eta && eta < 0.90) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609784) + 
    (0.90 <= eta && eta < 1.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.60) * (0.000000) + 
    (0.90 <= eta && eta < 1.00) * (0.60 <= pt * cosh(eta) && pt * cosh(eta) < 9.30) * (0.999390) + 
    (0.90 <= eta && eta < 1.00) * (9.30 <= pt * cosh(eta) && pt * cosh(eta) < 10.10) * (0.979355) + 
    (0.90 <= eta && eta < 1.00) * (10.10 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.957686) + 
    (0.90 <= eta && eta < 1.00) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.20) * (0.935743) + 
    (0.90 <= eta && eta < 1.00) * (11.20 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.913316) + 
    (0.90 <= eta && eta < 1.00) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.889743) + 
    (0.90 <= eta && eta < 1.00) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.866001) + 
    (0.90 <= eta && eta < 1.00) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.842837) + 
    (0.90 <= eta && eta < 1.00) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.820773) + 
    (0.90 <= eta && eta < 1.00) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.798174) + 
    (0.90 <= eta && eta < 1.00) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.775697) + 
    (0.90 <= eta && eta < 1.00) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.754029) + 
    (0.90 <= eta && eta < 1.00) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.732311) + 
    (0.90 <= eta && eta < 1.00) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.710471) + 
    (0.90 <= eta && eta < 1.00) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.688204) + 
    (0.90 <= eta && eta < 1.00) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.665960) + 
    (0.90 <= eta && eta < 1.00) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.643220) + 
    (0.90 <= eta && eta < 1.00) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609760) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 

# Everything else is not ID'd at all.
  add EfficiencyFormula {0} {0} { 0.00 }

}


#My name is "dualRICH_aerogel" and I am described as follows: 
#   Detector Type =  1 [Barrel=0, Forward=1]
#    Eta coverage =  [1,4]
#           Radii =  [9.1608925814663991,212.7295320598304]
#               Z =  250

module IdentificationMap dualRICH_aerogel { 
  set InputArray TrackSmearing/tracks
  set OutputArray tracks 

  # --- kaons ---

  add EfficiencyFormula {321} {321} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.00 <= eta && eta < 1.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.999421) + 
    (1.00 <= eta && eta < 1.10) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.980545) + 
    (1.00 <= eta && eta < 1.10) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.959608) + 
    (1.00 <= eta && eta < 1.10) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.937259) + 
    (1.00 <= eta && eta < 1.10) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.915885) + 
    (1.00 <= eta && eta < 1.10) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.895401) + 
    (1.00 <= eta && eta < 1.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.874634) + 
    (1.00 <= eta && eta < 1.10) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.854132) + 
    (1.00 <= eta && eta < 1.10) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.832382) + 
    (1.00 <= eta && eta < 1.10) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.809991) + 
    (1.00 <= eta && eta < 1.10) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.787646) + 
    (1.00 <= eta && eta < 1.10) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.765900) + 
    (1.00 <= eta && eta < 1.10) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.745301) + 
    (1.00 <= eta && eta < 1.10) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.724020) + 
    (1.00 <= eta && eta < 1.10) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.702466) + 
    (1.00 <= eta && eta < 1.10) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.680941) + 
    (1.00 <= eta && eta < 1.10) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.658802) + 
    (1.00 <= eta && eta < 1.10) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.635344) + 
    (1.00 <= eta && eta < 1.10) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610241) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.10 <= eta && eta < 1.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.999432) + 
    (1.10 <= eta && eta < 1.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.981085) + 
    (1.10 <= eta && eta < 1.20) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.961106) + 
    (1.10 <= eta && eta < 1.20) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.939739) + 
    (1.10 <= eta && eta < 1.20) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.917253) + 
    (1.10 <= eta && eta < 1.20) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.895414) + 
    (1.10 <= eta && eta < 1.20) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.875248) + 
    (1.10 <= eta && eta < 1.20) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.853363) + 
    (1.10 <= eta && eta < 1.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.830383) + 
    (1.10 <= eta && eta < 1.20) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.808750) + 
    (1.10 <= eta && eta < 1.20) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.787127) + 
    (1.10 <= eta && eta < 1.20) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.764636) + 
    (1.10 <= eta && eta < 1.20) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.742284) + 
    (1.10 <= eta && eta < 1.20) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.719967) + 
    (1.10 <= eta && eta < 1.20) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.698081) + 
    (1.10 <= eta && eta < 1.20) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.676291) + 
    (1.10 <= eta && eta < 1.20) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.654085) + 
    (1.10 <= eta && eta < 1.20) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.629900) + 
    (1.10 <= eta && eta < 1.20) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609371) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.20 <= eta && eta < 1.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.999418) + 
    (1.20 <= eta && eta < 1.30) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.979968) + 
    (1.20 <= eta && eta < 1.30) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.958724) + 
    (1.20 <= eta && eta < 1.30) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.937510) + 
    (1.20 <= eta && eta < 1.30) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.915434) + 
    (1.20 <= eta && eta < 1.30) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.892143) + 
    (1.20 <= eta && eta < 1.30) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.868604) + 
    (1.20 <= eta && eta < 1.30) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.845570) + 
    (1.20 <= eta && eta < 1.30) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.823568) + 
    (1.20 <= eta && eta < 1.30) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.801299) + 
    (1.20 <= eta && eta < 1.30) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.779332) + 
    (1.20 <= eta && eta < 1.30) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.758243) + 
    (1.20 <= eta && eta < 1.30) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.70) * (0.737276) + 
    (1.20 <= eta && eta < 1.30) * (19.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.716295) + 
    (1.20 <= eta && eta < 1.30) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.694892) + 
    (1.20 <= eta && eta < 1.30) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.673234) + 
    (1.20 <= eta && eta < 1.30) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.651069) + 
    (1.20 <= eta && eta < 1.30) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.626031) + 
    (1.20 <= eta && eta < 1.30) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608708) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.30 <= eta && eta < 1.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.999375) + 
    (1.30 <= eta && eta < 1.40) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.979302) + 
    (1.30 <= eta && eta < 1.40) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.958306) + 
    (1.30 <= eta && eta < 1.40) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.937530) + 
    (1.30 <= eta && eta < 1.40) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.915968) + 
    (1.30 <= eta && eta < 1.40) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.893222) + 
    (1.30 <= eta && eta < 1.40) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.870202) + 
    (1.30 <= eta && eta < 1.40) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.847622) + 
    (1.30 <= eta && eta < 1.40) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.825993) + 
    (1.30 <= eta && eta < 1.40) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.804027) + 
    (1.30 <= eta && eta < 1.40) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.782280) + 
    (1.30 <= eta && eta < 1.40) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.761316) + 
    (1.30 <= eta && eta < 1.40) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.740385) + 
    (1.30 <= eta && eta < 1.40) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.719339) + 
    (1.30 <= eta && eta < 1.40) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.697756) + 
    (1.30 <= eta && eta < 1.40) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.675783) + 
    (1.30 <= eta && eta < 1.40) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.653479) + 
    (1.30 <= eta && eta < 1.40) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.629191) + 
    (1.30 <= eta && eta < 1.40) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609919) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.40 <= eta && eta < 1.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.999393) + 
    (1.40 <= eta && eta < 1.50) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.979903) + 
    (1.40 <= eta && eta < 1.50) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.959632) + 
    (1.40 <= eta && eta < 1.50) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.939507) + 
    (1.40 <= eta && eta < 1.50) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.918534) + 
    (1.40 <= eta && eta < 1.50) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.896310) + 
    (1.40 <= eta && eta < 1.50) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.873716) + 
    (1.40 <= eta && eta < 1.50) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.851455) + 
    (1.40 <= eta && eta < 1.50) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.830037) + 
    (1.40 <= eta && eta < 1.50) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.80) * (0.808192) + 
    (1.40 <= eta && eta < 1.50) * (17.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.786471) + 
    (1.40 <= eta && eta < 1.50) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.765443) + 
    (1.40 <= eta && eta < 1.50) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.744357) + 
    (1.40 <= eta && eta < 1.50) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.723059) + 
    (1.40 <= eta && eta < 1.50) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.701115) + 
    (1.40 <= eta && eta < 1.50) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.679219) + 
    (1.40 <= eta && eta < 1.50) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.657259) + 
    (1.40 <= eta && eta < 1.50) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.633793) + 
    (1.40 <= eta && eta < 1.50) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611473) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.50 <= eta && eta < 1.60) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.999377) + 
    (1.50 <= eta && eta < 1.60) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.979737) + 
    (1.50 <= eta && eta < 1.60) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.959755) + 
    (1.50 <= eta && eta < 1.60) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.939987) + 
    (1.50 <= eta && eta < 1.60) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.919393) + 
    (1.50 <= eta && eta < 1.60) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.897551) + 
    (1.50 <= eta && eta < 1.60) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.875307) + 
    (1.50 <= eta && eta < 1.60) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.853345) + 
    (1.50 <= eta && eta < 1.60) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.832165) + 
    (1.50 <= eta && eta < 1.60) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.810506) + 
    (1.50 <= eta && eta < 1.60) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.788911) + 
    (1.50 <= eta && eta < 1.60) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.767943) + 
    (1.50 <= eta && eta < 1.60) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.746852) + 
    (1.50 <= eta && eta < 1.60) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.725478) + 
    (1.50 <= eta && eta < 1.60) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.704144) + 
    (1.50 <= eta && eta < 1.60) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.682446) + 
    (1.50 <= eta && eta < 1.60) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.660133) + 
    (1.50 <= eta && eta < 1.60) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.636792) + 
    (1.50 <= eta && eta < 1.60) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.612517) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.60 <= eta && eta < 1.70) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.999419) + 
    (1.60 <= eta && eta < 1.70) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.980743) + 
    (1.60 <= eta && eta < 1.70) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.960088) + 
    (1.60 <= eta && eta < 1.70) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.937415) + 
    (1.60 <= eta && eta < 1.70) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.915125) + 
    (1.60 <= eta && eta < 1.70) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.893414) + 
    (1.60 <= eta && eta < 1.70) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.871449) + 
    (1.60 <= eta && eta < 1.70) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.849856) + 
    (1.60 <= eta && eta < 1.70) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.829090) + 
    (1.60 <= eta && eta < 1.70) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.807889) + 
    (1.60 <= eta && eta < 1.70) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.785352) + 
    (1.60 <= eta && eta < 1.70) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.762438) + 
    (1.60 <= eta && eta < 1.70) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.740068) + 
    (1.60 <= eta && eta < 1.70) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.718126) + 
    (1.60 <= eta && eta < 1.70) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.696206) + 
    (1.60 <= eta && eta < 1.70) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.674404) + 
    (1.60 <= eta && eta < 1.70) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 28.80) * (0.652342) + 
    (1.60 <= eta && eta < 1.70) * (28.80 <= pt * cosh(eta) && pt * cosh(eta) < 36.40) * (0.627759) + 
    (1.60 <= eta && eta < 1.70) * (36.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610244) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.70 <= eta && eta < 1.80) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.999427) + 
    (1.70 <= eta && eta < 1.80) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.981008) + 
    (1.70 <= eta && eta < 1.80) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.960700) + 
    (1.70 <= eta && eta < 1.80) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.938370) + 
    (1.70 <= eta && eta < 1.80) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.916372) + 
    (1.70 <= eta && eta < 1.80) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.894899) + 
    (1.70 <= eta && eta < 1.80) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.873129) + 
    (1.70 <= eta && eta < 1.80) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.851684) + 
    (1.70 <= eta && eta < 1.80) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.80) * (0.829357) + 
    (1.70 <= eta && eta < 1.80) * (17.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.806754) + 
    (1.70 <= eta && eta < 1.80) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.784514) + 
    (1.70 <= eta && eta < 1.80) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.761892) + 
    (1.70 <= eta && eta < 1.80) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.739786) + 
    (1.70 <= eta && eta < 1.80) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.718072) + 
    (1.70 <= eta && eta < 1.80) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.696341) + 
    (1.70 <= eta && eta < 1.80) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.674673) + 
    (1.70 <= eta && eta < 1.80) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.00) * (0.652680) + 
    (1.70 <= eta && eta < 1.80) * (29.00 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.628190) + 
    (1.70 <= eta && eta < 1.80) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610525) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.80 <= eta && eta < 1.90) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.999406) + 
    (1.80 <= eta && eta < 1.90) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.980628) + 
    (1.80 <= eta && eta < 1.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.960286) + 
    (1.80 <= eta && eta < 1.90) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.938026) + 
    (1.80 <= eta && eta < 1.90) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.916141) + 
    (1.80 <= eta && eta < 1.90) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.894797) + 
    (1.80 <= eta && eta < 1.90) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.873164) + 
    (1.80 <= eta && eta < 1.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.851850) + 
    (1.80 <= eta && eta < 1.90) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.829651) + 
    (1.80 <= eta && eta < 1.90) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.807163) + 
    (1.80 <= eta && eta < 1.90) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.785020) + 
    (1.80 <= eta && eta < 1.90) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.30) * (0.762473) + 
    (1.80 <= eta && eta < 1.90) * (20.30 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.740416) + 
    (1.80 <= eta && eta < 1.90) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.718723) + 
    (1.80 <= eta && eta < 1.90) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.696981) + 
    (1.80 <= eta && eta < 1.90) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.675267) + 
    (1.80 <= eta && eta < 1.90) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.653182) + 
    (1.80 <= eta && eta < 1.90) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 36.40) * (0.628794) + 
    (1.80 <= eta && eta < 1.90) * (36.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610819) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.90 <= eta && eta < 2.00) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.999379) + 
    (1.90 <= eta && eta < 2.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.980126) + 
    (1.90 <= eta && eta < 2.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.959687) + 
    (1.90 <= eta && eta < 2.00) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.937449) + 
    (1.90 <= eta && eta < 2.00) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.915644) + 
    (1.90 <= eta && eta < 2.00) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.894409) + 
    (1.90 <= eta && eta < 2.00) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.872900) + 
    (1.90 <= eta && eta < 2.00) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.851713) + 
    (1.90 <= eta && eta < 2.00) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.829643) + 
    (1.90 <= eta && eta < 2.00) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.807280) + 
    (1.90 <= eta && eta < 2.00) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.785247) + 
    (1.90 <= eta && eta < 2.00) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.762796) + 
    (1.90 <= eta && eta < 2.00) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.740812) + 
    (1.90 <= eta && eta < 2.00) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.719168) + 
    (1.90 <= eta && eta < 2.00) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.697447) + 
    (1.90 <= eta && eta < 2.00) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.675721) + 
    (1.90 <= eta && eta < 2.00) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.653583) + 
    (1.90 <= eta && eta < 2.00) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.40) * (0.629212) + 
    (1.90 <= eta && eta < 2.00) * (36.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611042) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.00 <= eta && eta < 2.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.999430) + 
    (2.00 <= eta && eta < 2.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.980274) + 
    (2.00 <= eta && eta < 2.10) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.958791) + 
    (2.00 <= eta && eta < 2.10) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.936501) + 
    (2.00 <= eta && eta < 2.10) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.914729) + 
    (2.00 <= eta && eta < 2.10) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.893572) + 
    (2.00 <= eta && eta < 2.10) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.872170) + 
    (2.00 <= eta && eta < 2.10) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.851104) + 
    (2.00 <= eta && eta < 2.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.829169) + 
    (2.00 <= eta && eta < 2.10) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.806944) + 
    (2.00 <= eta && eta < 2.10) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.785043) + 
    (2.00 <= eta && eta < 2.10) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.762717) + 
    (2.00 <= eta && eta < 2.10) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.740844) + 
    (2.00 <= eta && eta < 2.10) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.719290) + 
    (2.00 <= eta && eta < 2.10) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.697639) + 
    (2.00 <= eta && eta < 2.10) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.675956) + 
    (2.00 <= eta && eta < 2.10) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.653825) + 
    (2.00 <= eta && eta < 2.10) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.40) * (0.629545) + 
    (2.00 <= eta && eta < 2.10) * (36.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611235) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.10 <= eta && eta < 2.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.999378) + 
    (2.10 <= eta && eta < 2.20) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.980201) + 
    (2.10 <= eta && eta < 2.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.959989) + 
    (2.10 <= eta && eta < 2.20) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.938001) + 
    (2.10 <= eta && eta < 2.20) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.916426) + 
    (2.10 <= eta && eta < 2.20) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.895389) + 
    (2.10 <= eta && eta < 2.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.874052) + 
    (2.10 <= eta && eta < 2.20) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.853005) + 
    (2.10 <= eta && eta < 2.20) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.831047) + 
    (2.10 <= eta && eta < 2.20) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.808762) + 
    (2.10 <= eta && eta < 2.20) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.786769) + 
    (2.10 <= eta && eta < 2.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.764319) + 
    (2.10 <= eta && eta < 2.20) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.742297) + 
    (2.10 <= eta && eta < 2.20) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.720574) + 
    (2.10 <= eta && eta < 2.20) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.698730) + 
    (2.10 <= eta && eta < 2.20) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.676833) + 
    (2.10 <= eta && eta < 2.20) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.654464) + 
    (2.10 <= eta && eta < 2.20) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.630191) + 
    (2.10 <= eta && eta < 2.20) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611490) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.20 <= eta && eta < 2.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.999408) + 
    (2.20 <= eta && eta < 2.30) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.979864) + 
    (2.20 <= eta && eta < 2.30) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.958336) + 
    (2.20 <= eta && eta < 2.30) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.936107) + 
    (2.20 <= eta && eta < 2.30) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.914440) + 
    (2.20 <= eta && eta < 2.30) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.893403) + 
    (2.20 <= eta && eta < 2.30) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.872130) + 
    (2.20 <= eta && eta < 2.30) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.851191) + 
    (2.20 <= eta && eta < 2.30) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.829380) + 
    (2.20 <= eta && eta < 2.30) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.807270) + 
    (2.20 <= eta && eta < 2.30) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.70) * (0.785465) + 
    (2.20 <= eta && eta < 2.30) * (19.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.763219) + 
    (2.20 <= eta && eta < 2.30) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.741400) + 
    (2.20 <= eta && eta < 2.30) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.719876) + 
    (2.20 <= eta && eta < 2.30) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.698225) + 
    (2.20 <= eta && eta < 2.30) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.676506) + 
    (2.20 <= eta && eta < 2.30) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.40) * (0.654297) + 
    (2.20 <= eta && eta < 2.30) * (29.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.30) * (0.630155) + 
    (2.20 <= eta && eta < 2.30) * (36.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611542) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.30 <= eta && eta < 2.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.999433) + 
    (2.30 <= eta && eta < 2.40) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.980427) + 
    (2.30 <= eta && eta < 2.40) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.959228) + 
    (2.30 <= eta && eta < 2.40) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.937219) + 
    (2.30 <= eta && eta < 2.40) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.915695) + 
    (2.30 <= eta && eta < 2.40) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.894746) + 
    (2.30 <= eta && eta < 2.40) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.873520) + 
    (2.30 <= eta && eta < 2.40) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.852593) + 
    (2.30 <= eta && eta < 2.40) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.830766) + 
    (2.30 <= eta && eta < 2.40) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.808612) + 
    (2.30 <= eta && eta < 2.40) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.70) * (0.786740) + 
    (2.30 <= eta && eta < 2.40) * (19.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.764402) + 
    (2.30 <= eta && eta < 2.40) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.742475) + 
    (2.30 <= eta && eta < 2.40) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.720826) + 
    (2.30 <= eta && eta < 2.40) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.699033) + 
    (2.30 <= eta && eta < 2.40) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.677157) + 
    (2.30 <= eta && eta < 2.40) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.655106) + 
    (2.30 <= eta && eta < 2.40) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 35.90) * (0.631229) + 
    (2.30 <= eta && eta < 2.40) * (35.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611899) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.40 <= eta && eta < 2.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.999428) + 
    (2.40 <= eta && eta < 2.50) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.980322) + 
    (2.40 <= eta && eta < 2.50) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.959060) + 
    (2.40 <= eta && eta < 2.50) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.937010) + 
    (2.40 <= eta && eta < 2.50) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.915458) + 
    (2.40 <= eta && eta < 2.50) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.894492) + 
    (2.40 <= eta && eta < 2.50) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.873257) + 
    (2.40 <= eta && eta < 2.50) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.852328) + 
    (2.40 <= eta && eta < 2.50) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.830504) + 
    (2.40 <= eta && eta < 2.50) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.808357) + 
    (2.40 <= eta && eta < 2.50) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.70) * (0.786498) + 
    (2.40 <= eta && eta < 2.50) * (19.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.764178) + 
    (2.40 <= eta && eta < 2.50) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.742271) + 
    (2.40 <= eta && eta < 2.50) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.720646) + 
    (2.40 <= eta && eta < 2.50) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.698879) + 
    (2.40 <= eta && eta < 2.50) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.677034) + 
    (2.40 <= eta && eta < 2.50) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.655015) + 
    (2.40 <= eta && eta < 2.50) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 35.90) * (0.631176) + 
    (2.40 <= eta && eta < 2.50) * (35.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611879) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.50 <= eta && eta < 2.60) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.999439) + 
    (2.50 <= eta && eta < 2.60) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.980370) + 
    (2.50 <= eta && eta < 2.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.958748) + 
    (2.50 <= eta && eta < 2.60) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.936272) + 
    (2.50 <= eta && eta < 2.60) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.914314) + 
    (2.50 <= eta && eta < 2.60) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.892985) + 
    (2.50 <= eta && eta < 2.60) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.871426) + 
    (2.50 <= eta && eta < 2.60) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.850227) + 
    (2.50 <= eta && eta < 2.60) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.828176) + 
    (2.50 <= eta && eta < 2.60) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.805862) + 
    (2.50 <= eta && eta < 2.60) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.783903) + 
    (2.50 <= eta && eta < 2.60) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.761550) + 
    (2.50 <= eta && eta < 2.60) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.739683) + 
    (2.50 <= eta && eta < 2.60) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.718173) + 
    (2.50 <= eta && eta < 2.60) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.696603) + 
    (2.50 <= eta && eta < 2.60) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.675044) + 
    (2.50 <= eta && eta < 2.60) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.653090) + 
    (2.50 <= eta && eta < 2.60) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.628801) + 
    (2.50 <= eta && eta < 2.60) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610884) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.60 <= eta && eta < 2.70) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.999380) + 
    (2.60 <= eta && eta < 2.70) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.979989) + 
    (2.60 <= eta && eta < 2.70) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.959094) + 
    (2.60 <= eta && eta < 2.70) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.936345) + 
    (2.60 <= eta && eta < 2.70) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.914076) + 
    (2.60 <= eta && eta < 2.70) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.892440) + 
    (2.60 <= eta && eta < 2.70) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.870584) + 
    (2.60 <= eta && eta < 2.70) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.849117) + 
    (2.60 <= eta && eta < 2.70) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.80) * (0.826825) + 
    (2.60 <= eta && eta < 2.70) * (17.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.804309) + 
    (2.60 <= eta && eta < 2.70) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.782201) + 
    (2.60 <= eta && eta < 2.70) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.759752) + 
    (2.60 <= eta && eta < 2.70) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.737851) + 
    (2.60 <= eta && eta < 2.70) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.716370) + 
    (2.60 <= eta && eta < 2.70) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.694900) + 
    (2.60 <= eta && eta < 2.70) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.673023) + 
    (2.60 <= eta && eta < 2.70) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.650533) + 
    (2.60 <= eta && eta < 2.70) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.625118) + 
    (2.60 <= eta && eta < 2.70) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609492) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.70 <= eta && eta < 2.80) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.999427) + 
    (2.70 <= eta && eta < 2.80) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.980834) + 
    (2.70 <= eta && eta < 2.80) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.960037) + 
    (2.70 <= eta && eta < 2.80) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.938840) + 
    (2.70 <= eta && eta < 2.80) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.918248) + 
    (2.70 <= eta && eta < 2.80) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.896469) + 
    (2.70 <= eta && eta < 2.80) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.874329) + 
    (2.70 <= eta && eta < 2.80) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.852492) + 
    (2.70 <= eta && eta < 2.80) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.831447) + 
    (2.70 <= eta && eta < 2.80) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.809933) + 
    (2.70 <= eta && eta < 2.80) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.788483) + 
    (2.70 <= eta && eta < 2.80) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.767651) + 
    (2.70 <= eta && eta < 2.80) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.746688) + 
    (2.70 <= eta && eta < 2.80) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.725431) + 
    (2.70 <= eta && eta < 2.80) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.703428) + 
    (2.70 <= eta && eta < 2.80) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.80) * (0.681357) + 
    (2.70 <= eta && eta < 2.80) * (24.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.659083) + 
    (2.70 <= eta && eta < 2.80) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.635658) + 
    (2.70 <= eta && eta < 2.80) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.612337) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.80 <= eta && eta < 2.90) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.999399) + 
    (2.80 <= eta && eta < 2.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.980128) + 
    (2.80 <= eta && eta < 2.90) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.960173) + 
    (2.80 <= eta && eta < 2.90) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.940341) + 
    (2.80 <= eta && eta < 2.90) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.919639) + 
    (2.80 <= eta && eta < 2.90) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.897662) + 
    (2.80 <= eta && eta < 2.90) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.875274) + 
    (2.80 <= eta && eta < 2.90) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.853171) + 
    (2.80 <= eta && eta < 2.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.831865) + 
    (2.80 <= eta && eta < 2.90) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.810090) + 
    (2.80 <= eta && eta < 2.90) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.788396) + 
    (2.80 <= eta && eta < 2.90) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.767351) + 
    (2.80 <= eta && eta < 2.90) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.30) * (0.746206) + 
    (2.80 <= eta && eta < 2.90) * (20.30 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.724804) + 
    (2.80 <= eta && eta < 2.90) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.703471) + 
    (2.80 <= eta && eta < 2.90) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.681809) + 
    (2.80 <= eta && eta < 2.90) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.659577) + 
    (2.80 <= eta && eta < 2.90) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.636380) + 
    (2.80 <= eta && eta < 2.90) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.612308) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.90 <= eta && eta < 3.00) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.999384) + 
    (2.90 <= eta && eta < 3.00) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.979697) + 
    (2.90 <= eta && eta < 3.00) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.959313) + 
    (2.90 <= eta && eta < 3.00) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.939113) + 
    (2.90 <= eta && eta < 3.00) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.918086) + 
    (2.90 <= eta && eta < 3.00) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.895827) + 
    (2.90 <= eta && eta < 3.00) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.873214) + 
    (2.90 <= eta && eta < 3.00) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.850948) + 
    (2.90 <= eta && eta < 3.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.829536) + 
    (2.90 <= eta && eta < 3.00) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.80) * (0.807706) + 
    (2.90 <= eta && eta < 3.00) * (17.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.786010) + 
    (2.90 <= eta && eta < 3.00) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.765013) + 
    (2.90 <= eta && eta < 3.00) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.743964) + 
    (2.90 <= eta && eta < 3.00) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.722710) + 
    (2.90 <= eta && eta < 3.00) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.700817) + 
    (2.90 <= eta && eta < 3.00) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.678978) + 
    (2.90 <= eta && eta < 3.00) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.657079) + 
    (2.90 <= eta && eta < 3.00) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.80) * (0.633504) + 
    (2.90 <= eta && eta < 3.00) * (32.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611367) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (3.00 <= eta && eta < 3.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.999380) + 
    (3.00 <= eta && eta < 3.10) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.979522) + 
    (3.00 <= eta && eta < 3.10) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.958847) + 
    (3.00 <= eta && eta < 3.10) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.938370) + 
    (3.00 <= eta && eta < 3.10) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.917083) + 
    (3.00 <= eta && eta < 3.10) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.894585) + 
    (3.00 <= eta && eta < 3.10) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.871771) + 
    (3.00 <= eta && eta < 3.10) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.849348) + 
    (3.00 <= eta && eta < 3.10) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.827826) + 
    (3.00 <= eta && eta < 3.10) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.805925) + 
    (3.00 <= eta && eta < 3.10) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.784199) + 
    (3.00 <= eta && eta < 3.10) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.20) * (0.763214) + 
    (3.00 <= eta && eta < 3.10) * (19.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.742218) + 
    (3.00 <= eta && eta < 3.10) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.721061) + 
    (3.00 <= eta && eta < 3.10) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.699316) + 
    (3.00 <= eta && eta < 3.10) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.677678) + 
    (3.00 <= eta && eta < 3.10) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.655671) + 
    (3.00 <= eta && eta < 3.10) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.631723) + 
    (3.00 <= eta && eta < 3.10) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610745) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (3.10 <= eta && eta < 3.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.999388) + 
    (3.10 <= eta && eta < 3.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.979585) + 
    (3.10 <= eta && eta < 3.20) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.958745) + 
    (3.10 <= eta && eta < 3.20) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.938071) + 
    (3.10 <= eta && eta < 3.20) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.916580) + 
    (3.10 <= eta && eta < 3.20) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.893880) + 
    (3.10 <= eta && eta < 3.20) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.870884) + 
    (3.10 <= eta && eta < 3.20) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.848309) + 
    (3.10 <= eta && eta < 3.20) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.826669) + 
    (3.10 <= eta && eta < 3.20) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.804680) + 
    (3.10 <= eta && eta < 3.20) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.782898) + 
    (3.10 <= eta && eta < 3.20) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.761891) + 
    (3.10 <= eta && eta < 3.20) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.740908) + 
    (3.10 <= eta && eta < 3.20) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.719802) + 
    (3.10 <= eta && eta < 3.20) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.698150) + 
    (3.10 <= eta && eta < 3.20) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.676099) + 
    (3.10 <= eta && eta < 3.20) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.653709) + 
    (3.10 <= eta && eta < 3.20) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.629318) + 
    (3.10 <= eta && eta < 3.20) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609960) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (3.20 <= eta && eta < 3.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.999404) + 
    (3.20 <= eta && eta < 3.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.979860) + 
    (3.20 <= eta && eta < 3.30) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.958971) + 
    (3.20 <= eta && eta < 3.30) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.938175) + 
    (3.20 <= eta && eta < 3.30) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.916531) + 
    (3.20 <= eta && eta < 3.30) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.893662) + 
    (3.20 <= eta && eta < 3.30) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.870500) + 
    (3.20 <= eta && eta < 3.30) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.847776) + 
    (3.20 <= eta && eta < 3.30) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.826010) + 
    (3.20 <= eta && eta < 3.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.803914) + 
    (3.20 <= eta && eta < 3.30) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.782050) + 
    (3.20 <= eta && eta < 3.30) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.760991) + 
    (3.20 <= eta && eta < 3.30) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.739984) + 
    (3.20 <= eta && eta < 3.30) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.718885) + 
    (3.20 <= eta && eta < 3.30) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.697276) + 
    (3.20 <= eta && eta < 3.30) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.675311) + 
    (3.20 <= eta && eta < 3.30) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.653061) + 
    (3.20 <= eta && eta < 3.30) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.10) * (0.628603) + 
    (3.20 <= eta && eta < 3.30) * (34.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609653) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (3.30 <= eta && eta < 3.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.999429) + 
    (3.30 <= eta && eta < 3.40) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.980322) + 
    (3.30 <= eta && eta < 3.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.959491) + 
    (3.30 <= eta && eta < 3.40) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.938643) + 
    (3.30 <= eta && eta < 3.40) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.916891) + 
    (3.30 <= eta && eta < 3.40) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.893882) + 
    (3.30 <= eta && eta < 3.40) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.870568) + 
    (3.30 <= eta && eta < 3.40) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.847696) + 
    (3.30 <= eta && eta < 3.40) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.825795) + 
    (3.30 <= eta && eta < 3.40) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.803576) + 
    (3.30 <= eta && eta < 3.40) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.781608) + 
    (3.30 <= eta && eta < 3.40) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.760467) + 
    (3.30 <= eta && eta < 3.40) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.739401) + 
    (3.30 <= eta && eta < 3.40) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.718270) + 
    (3.30 <= eta && eta < 3.40) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.696659) + 
    (3.30 <= eta && eta < 3.40) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.674731) + 
    (3.30 <= eta && eta < 3.40) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.652220) + 
    (3.30 <= eta && eta < 3.40) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.627361) + 
    (3.30 <= eta && eta < 3.40) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609223) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (3.40 <= eta && eta < 3.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.999367) + 
    (3.40 <= eta && eta < 3.50) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.980033) + 
    (3.40 <= eta && eta < 3.50) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.960271) + 
    (3.40 <= eta && eta < 3.50) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.939434) + 
    (3.40 <= eta && eta < 3.50) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.917619) + 
    (3.40 <= eta && eta < 3.50) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.894498) + 
    (3.40 <= eta && eta < 3.50) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.871045) + 
    (3.40 <= eta && eta < 3.50) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.848025) + 
    (3.40 <= eta && eta < 3.50) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.825982) + 
    (3.40 <= eta && eta < 3.50) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.803624) + 
    (3.40 <= eta && eta < 3.50) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.781528) + 
    (3.40 <= eta && eta < 3.50) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.760279) + 
    (3.40 <= eta && eta < 3.50) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.70) * (0.739123) + 
    (3.40 <= eta && eta < 3.50) * (19.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.717923) + 
    (3.40 <= eta && eta < 3.50) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.696270) + 
    (3.40 <= eta && eta < 3.50) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.674333) + 
    (3.40 <= eta && eta < 3.50) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.651857) + 
    (3.40 <= eta && eta < 3.50) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.626836) + 
    (3.40 <= eta && eta < 3.50) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608985) } 

  add EfficiencyFormula {321} {-211} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.00 <= eta && eta < 1.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.000723) + 
    (1.00 <= eta && eta < 1.10) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.020359) + 
    (1.00 <= eta && eta < 1.10) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.041268) + 
    (1.00 <= eta && eta < 1.10) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.062771) + 
    (1.00 <= eta && eta < 1.10) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.083795) + 
    (1.00 <= eta && eta < 1.10) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.104271) + 
    (1.00 <= eta && eta < 1.10) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.124913) + 
    (1.00 <= eta && eta < 1.10) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.146676) + 
    (1.00 <= eta && eta < 1.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.169069) + 
    (1.00 <= eta && eta < 1.10) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.190320) + 
    (1.00 <= eta && eta < 1.10) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.211461) + 
    (1.00 <= eta && eta < 1.10) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.233274) + 
    (1.00 <= eta && eta < 1.10) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.255188) + 
    (1.00 <= eta && eta < 1.10) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.276693) + 
    (1.00 <= eta && eta < 1.10) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.298134) + 
    (1.00 <= eta && eta < 1.10) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.319579) + 
    (1.00 <= eta && eta < 1.10) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.340806) + 
    (1.00 <= eta && eta < 1.10) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.361894) + 
    (1.00 <= eta && eta < 1.10) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.383172) + 
    (1.00 <= eta && eta < 1.10) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.404620) + 
    (1.00 <= eta && eta < 1.10) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 43.40) * (0.426353) + 
    (1.00 <= eta && eta < 1.10) * (43.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.444798) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.10 <= eta && eta < 1.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.000688) + 
    (1.10 <= eta && eta < 1.20) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.019438) + 
    (1.10 <= eta && eta < 1.20) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.039279) + 
    (1.10 <= eta && eta < 1.20) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.059801) + 
    (1.10 <= eta && eta < 1.20) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.081395) + 
    (1.10 <= eta && eta < 1.20) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.102638) + 
    (1.10 <= eta && eta < 1.20) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.124115) + 
    (1.10 <= eta && eta < 1.20) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.146671) + 
    (1.10 <= eta && eta < 1.20) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.168432) + 
    (1.10 <= eta && eta < 1.20) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.190363) + 
    (1.10 <= eta && eta < 1.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.212088) + 
    (1.10 <= eta && eta < 1.20) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.233248) + 
    (1.10 <= eta && eta < 1.20) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.254563) + 
    (1.10 <= eta && eta < 1.20) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.275549) + 
    (1.10 <= eta && eta < 1.20) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.296550) + 
    (1.10 <= eta && eta < 1.20) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.317640) + 
    (1.10 <= eta && eta < 1.20) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.338609) + 
    (1.10 <= eta && eta < 1.20) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.359540) + 
    (1.10 <= eta && eta < 1.20) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.380769) + 
    (1.10 <= eta && eta < 1.20) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.402289) + 
    (1.10 <= eta && eta < 1.20) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 43.90) * (0.424061) + 
    (1.10 <= eta && eta < 1.20) * (43.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.442161) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.20 <= eta && eta < 1.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.000680) + 
    (1.20 <= eta && eta < 1.30) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.019858) + 
    (1.20 <= eta && eta < 1.30) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.040284) + 
    (1.20 <= eta && eta < 1.30) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.060404) + 
    (1.20 <= eta && eta < 1.30) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.081468) + 
    (1.20 <= eta && eta < 1.30) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.103556) + 
    (1.20 <= eta && eta < 1.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.125865) + 
    (1.20 <= eta && eta < 1.30) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.147802) + 
    (1.20 <= eta && eta < 1.30) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.20) * (0.168967) + 
    (1.20 <= eta && eta < 1.30) * (19.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.190319) + 
    (1.20 <= eta && eta < 1.30) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.211508) + 
    (1.20 <= eta && eta < 1.30) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.232191) + 
    (1.20 <= eta && eta < 1.30) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.253082) + 
    (1.20 <= eta && eta < 1.30) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.274542) + 
    (1.20 <= eta && eta < 1.30) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.295913) + 
    (1.20 <= eta && eta < 1.30) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.317209) + 
    (1.20 <= eta && eta < 1.30) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.338823) + 
    (1.20 <= eta && eta < 1.30) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.360107) + 
    (1.20 <= eta && eta < 1.30) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.00) * (0.381426) + 
    (1.20 <= eta && eta < 1.30) * (35.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.10) * (0.402823) + 
    (1.20 <= eta && eta < 1.30) * (39.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.20) * (0.424492) + 
    (1.20 <= eta && eta < 1.30) * (45.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.440832) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.30 <= eta && eta < 1.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.000701) + 
    (1.30 <= eta && eta < 1.40) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.019989) + 
    (1.30 <= eta && eta < 1.40) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.041052) + 
    (1.30 <= eta && eta < 1.40) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.062014) + 
    (1.30 <= eta && eta < 1.40) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.082706) + 
    (1.30 <= eta && eta < 1.40) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.104326) + 
    (1.30 <= eta && eta < 1.40) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.126133) + 
    (1.30 <= eta && eta < 1.40) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.147579) + 
    (1.30 <= eta && eta < 1.40) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.168289) + 
    (1.30 <= eta && eta < 1.40) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.189215) + 
    (1.30 <= eta && eta < 1.40) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.211120) + 
    (1.30 <= eta && eta < 1.40) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.232423) + 
    (1.30 <= eta && eta < 1.40) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.252847) + 
    (1.30 <= eta && eta < 1.40) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.273879) + 
    (1.30 <= eta && eta < 1.40) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.295594) + 
    (1.30 <= eta && eta < 1.40) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.317131) + 
    (1.30 <= eta && eta < 1.40) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.338314) + 
    (1.30 <= eta && eta < 1.40) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.60) * (0.359675) + 
    (1.30 <= eta && eta < 1.40) * (32.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.70) * (0.381014) + 
    (1.30 <= eta && eta < 1.40) * (35.70 <= pt * cosh(eta) && pt * cosh(eta) < 39.90) * (0.402342) + 
    (1.30 <= eta && eta < 1.40) * (39.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.10) * (0.424107) + 
    (1.30 <= eta && eta < 1.40) * (46.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.439317) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.40 <= eta && eta < 1.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.000671) + 
    (1.40 <= eta && eta < 1.50) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.019241) + 
    (1.40 <= eta && eta < 1.50) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.039529) + 
    (1.40 <= eta && eta < 1.50) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.061058) + 
    (1.40 <= eta && eta < 1.50) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.082565) + 
    (1.40 <= eta && eta < 1.50) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.103768) + 
    (1.40 <= eta && eta < 1.50) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.125167) + 
    (1.40 <= eta && eta < 1.50) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.146237) + 
    (1.40 <= eta && eta < 1.50) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.167853) + 
    (1.40 <= eta && eta < 1.50) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.189600) + 
    (1.40 <= eta && eta < 1.50) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.211079) + 
    (1.40 <= eta && eta < 1.50) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.232986) + 
    (1.40 <= eta && eta < 1.50) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.254800) + 
    (1.40 <= eta && eta < 1.50) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.276071) + 
    (1.40 <= eta && eta < 1.50) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.297134) + 
    (1.40 <= eta && eta < 1.50) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.318082) + 
    (1.40 <= eta && eta < 1.50) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.339249) + 
    (1.40 <= eta && eta < 1.50) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.30) * (0.360510) + 
    (1.40 <= eta && eta < 1.50) * (33.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.381635) + 
    (1.40 <= eta && eta < 1.50) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.80) * (0.402956) + 
    (1.40 <= eta && eta < 1.50) * (40.80 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.424693) + 
    (1.40 <= eta && eta < 1.50) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.438435) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.50 <= eta && eta < 1.60) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.000673) + 
    (1.50 <= eta && eta < 1.60) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.019124) + 
    (1.50 <= eta && eta < 1.60) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.039002) + 
    (1.50 <= eta && eta < 1.60) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.060084) + 
    (1.50 <= eta && eta < 1.60) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.081176) + 
    (1.50 <= eta && eta < 1.60) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.102012) + 
    (1.50 <= eta && eta < 1.60) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.123085) + 
    (1.50 <= eta && eta < 1.60) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.143882) + 
    (1.50 <= eta && eta < 1.60) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.165267) + 
    (1.50 <= eta && eta < 1.60) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.186831) + 
    (1.50 <= eta && eta < 1.60) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.208178) + 
    (1.50 <= eta && eta < 1.60) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.230000) + 
    (1.50 <= eta && eta < 1.60) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.251777) + 
    (1.50 <= eta && eta < 1.60) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.273060) + 
    (1.50 <= eta && eta < 1.60) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.294179) + 
    (1.50 <= eta && eta < 1.60) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.90) * (0.315226) + 
    (1.50 <= eta && eta < 1.60) * (28.90 <= pt * cosh(eta) && pt * cosh(eta) < 30.90) * (0.336037) + 
    (1.50 <= eta && eta < 1.60) * (30.90 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.357142) + 
    (1.50 <= eta && eta < 1.60) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.378356) + 
    (1.50 <= eta && eta < 1.60) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.399705) + 
    (1.50 <= eta && eta < 1.60) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.421490) + 
    (1.50 <= eta && eta < 1.60) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.435963) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.60 <= eta && eta < 1.70) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.000707) + 
    (1.60 <= eta && eta < 1.70) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.019593) + 
    (1.60 <= eta && eta < 1.70) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.039389) + 
    (1.60 <= eta && eta < 1.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.060261) + 
    (1.60 <= eta && eta < 1.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.081099) + 
    (1.60 <= eta && eta < 1.70) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.101674) + 
    (1.60 <= eta && eta < 1.70) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.122490) + 
    (1.60 <= eta && eta < 1.70) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.143048) + 
    (1.60 <= eta && eta < 1.70) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.164209) + 
    (1.60 <= eta && eta < 1.70) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.185573) + 
    (1.60 <= eta && eta < 1.70) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.206751) + 
    (1.60 <= eta && eta < 1.70) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 23.40) * (0.228432) + 
    (1.60 <= eta && eta < 1.70) * (23.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.250102) + 
    (1.60 <= eta && eta < 1.70) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.271314) + 
    (1.60 <= eta && eta < 1.70) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.292398) + 
    (1.60 <= eta && eta < 1.70) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.313444) + 
    (1.60 <= eta && eta < 1.70) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.334290) + 
    (1.60 <= eta && eta < 1.70) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.355467) + 
    (1.60 <= eta && eta < 1.70) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.376793) + 
    (1.60 <= eta && eta < 1.70) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.80) * (0.398050) + 
    (1.60 <= eta && eta < 1.70) * (40.80 <= pt * cosh(eta) && pt * cosh(eta) < 46.70) * (0.419582) + 
    (1.60 <= eta && eta < 1.70) * (46.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.434240) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.70 <= eta && eta < 1.80) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.000692) + 
    (1.70 <= eta && eta < 1.80) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.020003) + 
    (1.70 <= eta && eta < 1.80) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.040686) + 
    (1.70 <= eta && eta < 1.80) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.061566) + 
    (1.70 <= eta && eta < 1.80) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.082303) + 
    (1.70 <= eta && eta < 1.80) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.102719) + 
    (1.70 <= eta && eta < 1.80) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.123343) + 
    (1.70 <= eta && eta < 1.80) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.144946) + 
    (1.70 <= eta && eta < 1.80) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.167051) + 
    (1.70 <= eta && eta < 1.80) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.188072) + 
    (1.70 <= eta && eta < 1.80) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.208898) + 
    (1.70 <= eta && eta < 1.80) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.230217) + 
    (1.70 <= eta && eta < 1.80) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.251534) + 
    (1.70 <= eta && eta < 1.80) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.272414) + 
    (1.70 <= eta && eta < 1.80) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 27.70) * (0.293188) + 
    (1.70 <= eta && eta < 1.80) * (27.70 <= pt * cosh(eta) && pt * cosh(eta) < 29.40) * (0.313951) + 
    (1.70 <= eta && eta < 1.80) * (29.40 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.335044) + 
    (1.70 <= eta && eta < 1.80) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.356349) + 
    (1.70 <= eta && eta < 1.80) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 37.20) * (0.377645) + 
    (1.70 <= eta && eta < 1.80) * (37.20 <= pt * cosh(eta) && pt * cosh(eta) < 41.40) * (0.399035) + 
    (1.70 <= eta && eta < 1.80) * (41.40 <= pt * cosh(eta) && pt * cosh(eta) < 47.50) * (0.420636) + 
    (1.70 <= eta && eta < 1.80) * (47.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.434156) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.80 <= eta && eta < 1.90) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.000706) + 
    (1.80 <= eta && eta < 1.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.020190) + 
    (1.80 <= eta && eta < 1.90) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.040805) + 
    (1.80 <= eta && eta < 1.90) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.061565) + 
    (1.80 <= eta && eta < 1.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.082168) + 
    (1.80 <= eta && eta < 1.90) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.102452) + 
    (1.80 <= eta && eta < 1.90) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.122948) + 
    (1.80 <= eta && eta < 1.90) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.144427) + 
    (1.80 <= eta && eta < 1.90) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.166418) + 
    (1.80 <= eta && eta < 1.90) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.70) * (0.187346) + 
    (1.80 <= eta && eta < 1.90) * (21.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.208094) + 
    (1.80 <= eta && eta < 1.90) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.229350) + 
    (1.80 <= eta && eta < 1.90) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.250620) + 
    (1.80 <= eta && eta < 1.90) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.271471) + 
    (1.80 <= eta && eta < 1.90) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.292233) + 
    (1.80 <= eta && eta < 1.90) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.313003) + 
    (1.80 <= eta && eta < 1.90) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.334119) + 
    (1.80 <= eta && eta < 1.90) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.10) * (0.355468) + 
    (1.80 <= eta && eta < 1.90) * (34.10 <= pt * cosh(eta) && pt * cosh(eta) < 37.30) * (0.376826) + 
    (1.80 <= eta && eta < 1.90) * (37.30 <= pt * cosh(eta) && pt * cosh(eta) < 41.50) * (0.398300) + 
    (1.80 <= eta && eta < 1.90) * (41.50 <= pt * cosh(eta) && pt * cosh(eta) < 47.50) * (0.419840) + 
    (1.80 <= eta && eta < 1.90) * (47.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.433334) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.90 <= eta && eta < 2.00) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.000726) + 
    (1.90 <= eta && eta < 2.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.020479) + 
    (1.90 <= eta && eta < 2.00) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.041079) + 
    (1.90 <= eta && eta < 2.00) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.061756) + 
    (1.90 <= eta && eta < 2.00) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.082251) + 
    (1.90 <= eta && eta < 2.00) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.102420) + 
    (1.90 <= eta && eta < 2.00) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.122800) + 
    (1.90 <= eta && eta < 2.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.144161) + 
    (1.90 <= eta && eta < 2.00) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.166041) + 
    (1.90 <= eta && eta < 2.00) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.186873) + 
    (1.90 <= eta && eta < 2.00) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.207538) + 
    (1.90 <= eta && eta < 2.00) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.228721) + 
    (1.90 <= eta && eta < 2.00) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.249934) + 
    (1.90 <= eta && eta < 2.00) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.270743) + 
    (1.90 <= eta && eta < 2.00) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.291479) + 
    (1.90 <= eta && eta < 2.00) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.60) * (0.312237) + 
    (1.90 <= eta && eta < 2.00) * (29.60 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.333360) + 
    (1.90 <= eta && eta < 2.00) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.354732) + 
    (1.90 <= eta && eta < 2.00) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 37.30) * (0.375809) + 
    (1.90 <= eta && eta < 2.00) * (37.30 <= pt * cosh(eta) && pt * cosh(eta) < 41.40) * (0.396913) + 
    (1.90 <= eta && eta < 2.00) * (41.40 <= pt * cosh(eta) && pt * cosh(eta) < 47.30) * (0.418399) + 
    (1.90 <= eta && eta < 2.00) * (47.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.432293) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.00 <= eta && eta < 2.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.000675) + 
    (2.00 <= eta && eta < 2.10) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.019521) + 
    (2.00 <= eta && eta < 2.10) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.039607) + 
    (2.00 <= eta && eta < 2.10) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.059943) + 
    (2.00 <= eta && eta < 2.10) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.081447) + 
    (2.00 <= eta && eta < 2.10) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.102754) + 
    (2.00 <= eta && eta < 2.10) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.20) * (0.123034) + 
    (2.00 <= eta && eta < 2.10) * (19.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.144287) + 
    (2.00 <= eta && eta < 2.10) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.166057) + 
    (2.00 <= eta && eta < 2.10) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.187897) + 
    (2.00 <= eta && eta < 2.10) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.209448) + 
    (2.00 <= eta && eta < 2.10) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.230382) + 
    (2.00 <= eta && eta < 2.10) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.251344) + 
    (2.00 <= eta && eta < 2.10) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 26.70) * (0.272676) + 
    (2.00 <= eta && eta < 2.10) * (26.70 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.293790) + 
    (2.00 <= eta && eta < 2.10) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.314730) + 
    (2.00 <= eta && eta < 2.10) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.335893) + 
    (2.00 <= eta && eta < 2.10) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.357125) + 
    (2.00 <= eta && eta < 2.10) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.378596) + 
    (2.00 <= eta && eta < 2.10) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 42.40) * (0.400170) + 
    (2.00 <= eta && eta < 2.10) * (42.40 <= pt * cosh(eta) && pt * cosh(eta) < 48.70) * (0.421818) + 
    (2.00 <= eta && eta < 2.10) * (48.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.433583) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.10 <= eta && eta < 2.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.000720) + 
    (2.10 <= eta && eta < 2.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.020280) + 
    (2.10 <= eta && eta < 2.20) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.040610) + 
    (2.10 <= eta && eta < 2.20) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.061035) + 
    (2.10 <= eta && eta < 2.20) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.082548) + 
    (2.10 <= eta && eta < 2.20) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.103815) + 
    (2.10 <= eta && eta < 2.20) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.124026) + 
    (2.10 <= eta && eta < 2.20) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.145187) + 
    (2.10 <= eta && eta < 2.20) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.166851) + 
    (2.10 <= eta && eta < 2.20) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.188577) + 
    (2.10 <= eta && eta < 2.20) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.210014) + 
    (2.10 <= eta && eta < 2.20) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.230838) + 
    (2.10 <= eta && eta < 2.20) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.251693) + 
    (2.10 <= eta && eta < 2.20) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.272923) + 
    (2.10 <= eta && eta < 2.20) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.40) * (0.294608) + 
    (2.10 <= eta && eta < 2.20) * (28.40 <= pt * cosh(eta) && pt * cosh(eta) < 30.20) * (0.315972) + 
    (2.10 <= eta && eta < 2.20) * (30.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.336881) + 
    (2.10 <= eta && eta < 2.20) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.357875) + 
    (2.10 <= eta && eta < 2.20) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.379129) + 
    (2.10 <= eta && eta < 2.20) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 42.60) * (0.400510) + 
    (2.10 <= eta && eta < 2.20) * (42.60 <= pt * cosh(eta) && pt * cosh(eta) < 49.00) * (0.422156) + 
    (2.10 <= eta && eta < 2.20) * (49.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.433561) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.20 <= eta && eta < 2.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.000691) + 
    (2.20 <= eta && eta < 2.30) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.019733) + 
    (2.20 <= eta && eta < 2.30) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.039771) + 
    (2.20 <= eta && eta < 2.30) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.060002) + 
    (2.20 <= eta && eta < 2.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.081377) + 
    (2.20 <= eta && eta < 2.30) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.102553) + 
    (2.20 <= eta && eta < 2.30) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.122711) + 
    (2.20 <= eta && eta < 2.30) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.143846) + 
    (2.20 <= eta && eta < 2.30) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.165506) + 
    (2.20 <= eta && eta < 2.30) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.187249) + 
    (2.20 <= eta && eta < 2.30) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.208720) + 
    (2.20 <= eta && eta < 2.30) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.229590) + 
    (2.20 <= eta && eta < 2.30) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.250504) + 
    (2.20 <= eta && eta < 2.30) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.271804) + 
    (2.20 <= eta && eta < 2.30) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.292903) + 
    (2.20 <= eta && eta < 2.30) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.313845) + 
    (2.20 <= eta && eta < 2.30) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.335028) + 
    (2.20 <= eta && eta < 2.30) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.80) * (0.356297) + 
    (2.20 <= eta && eta < 2.30) * (34.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.377513) + 
    (2.20 <= eta && eta < 2.30) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.398751) + 
    (2.20 <= eta && eta < 2.30) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 48.50) * (0.420375) + 
    (2.20 <= eta && eta < 2.30) * (48.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.432512) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.30 <= eta && eta < 2.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.000666) + 
    (2.30 <= eta && eta < 2.40) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.019270) + 
    (2.30 <= eta && eta < 2.40) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.039057) + 
    (2.30 <= eta && eta < 2.40) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.059121) + 
    (2.30 <= eta && eta < 2.40) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.080375) + 
    (2.30 <= eta && eta < 2.40) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.101472) + 
    (2.30 <= eta && eta < 2.40) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.121584) + 
    (2.30 <= eta && eta < 2.40) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.142694) + 
    (2.30 <= eta && eta < 2.40) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.164350) + 
    (2.30 <= eta && eta < 2.40) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.186107) + 
    (2.30 <= eta && eta < 2.40) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.207606) + 
    (2.30 <= eta && eta < 2.40) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.228515) + 
    (2.30 <= eta && eta < 2.40) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.249479) + 
    (2.30 <= eta && eta < 2.40) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.270839) + 
    (2.30 <= eta && eta < 2.40) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.292006) + 
    (2.30 <= eta && eta < 2.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.313021) + 
    (2.30 <= eta && eta < 2.40) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.334284) + 
    (2.30 <= eta && eta < 2.40) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.355239) + 
    (2.30 <= eta && eta < 2.40) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 37.90) * (0.376286) + 
    (2.30 <= eta && eta < 2.40) * (37.90 <= pt * cosh(eta) && pt * cosh(eta) < 42.10) * (0.397538) + 
    (2.30 <= eta && eta < 2.40) * (42.10 <= pt * cosh(eta) && pt * cosh(eta) < 48.20) * (0.419128) + 
    (2.30 <= eta && eta < 2.40) * (48.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.431773) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.40 <= eta && eta < 2.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.000671) + 
    (2.40 <= eta && eta < 2.50) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.019357) + 
    (2.40 <= eta && eta < 2.50) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.039192) + 
    (2.40 <= eta && eta < 2.50) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.059287) + 
    (2.40 <= eta && eta < 2.50) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.080564) + 
    (2.40 <= eta && eta < 2.50) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.101676) + 
    (2.40 <= eta && eta < 2.50) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.121797) + 
    (2.40 <= eta && eta < 2.50) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.142912) + 
    (2.40 <= eta && eta < 2.50) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.164569) + 
    (2.40 <= eta && eta < 2.50) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.186323) + 
    (2.40 <= eta && eta < 2.50) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.207817) + 
    (2.40 <= eta && eta < 2.50) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.228719) + 
    (2.40 <= eta && eta < 2.50) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.249673) + 
    (2.40 <= eta && eta < 2.50) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.271022) + 
    (2.40 <= eta && eta < 2.50) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.292176) + 
    (2.40 <= eta && eta < 2.50) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.313177) + 
    (2.40 <= eta && eta < 2.50) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.334425) + 
    (2.40 <= eta && eta < 2.50) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.80) * (0.355764) + 
    (2.40 <= eta && eta < 2.50) * (34.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.377053) + 
    (2.40 <= eta && eta < 2.50) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.398367) + 
    (2.40 <= eta && eta < 2.50) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 48.50) * (0.420069) + 
    (2.40 <= eta && eta < 2.50) * (48.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.432252) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.50 <= eta && eta < 2.60) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.000673) + 
    (2.50 <= eta && eta < 2.60) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.019568) + 
    (2.50 <= eta && eta < 2.60) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.039845) + 
    (2.50 <= eta && eta < 2.60) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.060381) + 
    (2.50 <= eta && eta < 2.60) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.080832) + 
    (2.50 <= eta && eta < 2.60) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.101016) + 
    (2.50 <= eta && eta < 2.60) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.121450) + 
    (2.50 <= eta && eta < 2.60) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.142896) + 
    (2.50 <= eta && eta < 2.60) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.164882) + 
    (2.50 <= eta && eta < 2.60) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.70) * (0.185828) + 
    (2.50 <= eta && eta < 2.60) * (21.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.206612) + 
    (2.50 <= eta && eta < 2.60) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.227920) + 
    (2.50 <= eta && eta < 2.60) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.249258) + 
    (2.50 <= eta && eta < 2.60) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.270187) + 
    (2.50 <= eta && eta < 2.60) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.291037) + 
    (2.50 <= eta && eta < 2.60) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.311903) + 
    (2.50 <= eta && eta < 2.60) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.332624) + 
    (2.50 <= eta && eta < 2.60) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.353737) + 
    (2.50 <= eta && eta < 2.60) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.375063) + 
    (2.50 <= eta && eta < 2.60) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 41.20) * (0.396391) + 
    (2.50 <= eta && eta < 2.60) * (41.20 <= pt * cosh(eta) && pt * cosh(eta) < 47.00) * (0.417902) + 
    (2.50 <= eta && eta < 2.60) * (47.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.432206) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.60 <= eta && eta < 2.70) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.000739) + 
    (2.60 <= eta && eta < 2.70) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.020090) + 
    (2.60 <= eta && eta < 2.70) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.039984) + 
    (2.60 <= eta && eta < 2.70) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.060848) + 
    (2.60 <= eta && eta < 2.70) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.081627) + 
    (2.60 <= eta && eta < 2.70) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.102117) + 
    (2.60 <= eta && eta < 2.70) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.122833) + 
    (2.60 <= eta && eta < 2.70) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.143287) + 
    (2.60 <= eta && eta < 2.70) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.164340) + 
    (2.60 <= eta && eta < 2.70) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.185599) + 
    (2.60 <= eta && eta < 2.70) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.206678) + 
    (2.60 <= eta && eta < 2.70) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.228265) + 
    (2.60 <= eta && eta < 2.70) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.249854) + 
    (2.60 <= eta && eta < 2.70) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.270996) + 
    (2.60 <= eta && eta < 2.70) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.292023) + 
    (2.60 <= eta && eta < 2.70) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.313027) + 
    (2.60 <= eta && eta < 2.70) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 31.20) * (0.333844) + 
    (2.60 <= eta && eta < 2.70) * (31.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.355009) + 
    (2.60 <= eta && eta < 2.70) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.376338) + 
    (2.60 <= eta && eta < 2.70) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 40.90) * (0.397617) + 
    (2.60 <= eta && eta < 2.70) * (40.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.419190) + 
    (2.60 <= eta && eta < 2.70) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.433757) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.70 <= eta && eta < 2.80) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.000706) + 
    (2.70 <= eta && eta < 2.80) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.019651) + 
    (2.70 <= eta && eta < 2.80) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.039646) + 
    (2.70 <= eta && eta < 2.80) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.060732) + 
    (2.70 <= eta && eta < 2.80) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.081770) + 
    (2.70 <= eta && eta < 2.80) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 17.80) * (0.102523) + 
    (2.70 <= eta && eta < 2.80) * (17.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.123497) + 
    (2.70 <= eta && eta < 2.80) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.144188) + 
    (2.70 <= eta && eta < 2.80) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.30) * (0.165463) + 
    (2.70 <= eta && eta < 2.80) * (20.30 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.186919) + 
    (2.70 <= eta && eta < 2.80) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.208163) + 
    (2.70 <= eta && eta < 2.80) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.30) * (0.229888) + 
    (2.70 <= eta && eta < 2.80) * (23.30 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.251580) + 
    (2.70 <= eta && eta < 2.80) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.272789) + 
    (2.70 <= eta && eta < 2.80) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.293848) + 
    (2.70 <= eta && eta < 2.80) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.00) * (0.314849) + 
    (2.70 <= eta && eta < 2.80) * (29.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.335628) + 
    (2.70 <= eta && eta < 2.80) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 33.50) * (0.356716) + 
    (2.70 <= eta && eta < 2.80) * (33.50 <= pt * cosh(eta) && pt * cosh(eta) < 36.60) * (0.377930) + 
    (2.70 <= eta && eta < 2.80) * (36.60 <= pt * cosh(eta) && pt * cosh(eta) < 40.80) * (0.399297) + 
    (2.70 <= eta && eta < 2.80) * (40.80 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.420954) + 
    (2.70 <= eta && eta < 2.80) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.435367) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.80 <= eta && eta < 2.90) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.000742) + 
    (2.80 <= eta && eta < 2.90) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.020385) + 
    (2.80 <= eta && eta < 2.90) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.040943) + 
    (2.80 <= eta && eta < 2.90) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.062496) + 
    (2.80 <= eta && eta < 2.90) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.083904) + 
    (2.80 <= eta && eta < 2.90) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.104944) + 
    (2.80 <= eta && eta < 2.90) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.126142) + 
    (2.80 <= eta && eta < 2.90) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.146998) + 
    (2.80 <= eta && eta < 2.90) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.168389) + 
    (2.80 <= eta && eta < 2.90) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.189914) + 
    (2.80 <= eta && eta < 2.90) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.211182) + 
    (2.80 <= eta && eta < 2.90) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.232890) + 
    (2.80 <= eta && eta < 2.90) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.254524) + 
    (2.80 <= eta && eta < 2.90) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.275642) + 
    (2.80 <= eta && eta < 2.90) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.296577) + 
    (2.80 <= eta && eta < 2.90) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.90) * (0.317424) + 
    (2.80 <= eta && eta < 2.90) * (28.90 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.338518) + 
    (2.80 <= eta && eta < 2.90) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 33.50) * (0.359737) + 
    (2.80 <= eta && eta < 2.90) * (33.50 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.380854) + 
    (2.80 <= eta && eta < 2.90) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 41.00) * (0.402203) + 
    (2.80 <= eta && eta < 2.90) * (41.00 <= pt * cosh(eta) && pt * cosh(eta) < 47.30) * (0.423852) + 
    (2.80 <= eta && eta < 2.90) * (47.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.437441) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.90 <= eta && eta < 3.00) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.000680) + 
    (2.90 <= eta && eta < 3.00) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.019408) + 
    (2.90 <= eta && eta < 3.00) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.039788) + 
    (2.90 <= eta && eta < 3.00) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.061380) + 
    (2.90 <= eta && eta < 3.00) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.082930) + 
    (2.90 <= eta && eta < 3.00) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.104161) + 
    (2.90 <= eta && eta < 3.00) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.125576) + 
    (2.90 <= eta && eta < 3.00) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.146654) + 
    (2.90 <= eta && eta < 3.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.168271) + 
    (2.90 <= eta && eta < 3.00) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.190012) + 
    (2.90 <= eta && eta < 3.00) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.211480) + 
    (2.90 <= eta && eta < 3.00) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.233371) + 
    (2.90 <= eta && eta < 3.00) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.255165) + 
    (2.90 <= eta && eta < 3.00) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.276414) + 
    (2.90 <= eta && eta < 3.00) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.297452) + 
    (2.90 <= eta && eta < 3.00) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.318373) + 
    (2.90 <= eta && eta < 3.00) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.339511) + 
    (2.90 <= eta && eta < 3.00) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.30) * (0.360741) + 
    (2.90 <= eta && eta < 3.00) * (33.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.381833) + 
    (2.90 <= eta && eta < 3.00) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.80) * (0.403120) + 
    (2.90 <= eta && eta < 3.00) * (40.80 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.424822) + 
    (2.90 <= eta && eta < 3.00) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.438540) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (3.00 <= eta && eta < 3.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.000689) + 
    (3.00 <= eta && eta < 3.10) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.019678) + 
    (3.00 <= eta && eta < 3.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.040386) + 
    (3.00 <= eta && eta < 3.10) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.062288) + 
    (3.00 <= eta && eta < 3.10) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.084104) + 
    (3.00 <= eta && eta < 3.10) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.105555) + 
    (3.00 <= eta && eta < 3.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.127155) + 
    (3.00 <= eta && eta < 3.10) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.148379) + 
    (3.00 <= eta && eta < 3.10) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.170110) + 
    (3.00 <= eta && eta < 3.10) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.191935) + 
    (3.00 <= eta && eta < 3.10) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.213453) + 
    (3.00 <= eta && eta < 3.10) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.90) * (0.235364) + 
    (3.00 <= eta && eta < 3.10) * (22.90 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.257149) + 
    (3.00 <= eta && eta < 3.10) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.278361) + 
    (3.00 <= eta && eta < 3.10) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.90) * (0.299337) + 
    (3.00 <= eta && eta < 3.10) * (26.90 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.320759) + 
    (3.00 <= eta && eta < 3.10) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.342209) + 
    (3.00 <= eta && eta < 3.10) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.363512) + 
    (3.00 <= eta && eta < 3.10) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.384827) + 
    (3.00 <= eta && eta < 3.10) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 41.20) * (0.406224) + 
    (3.00 <= eta && eta < 3.10) * (41.20 <= pt * cosh(eta) && pt * cosh(eta) < 48.00) * (0.428061) + 
    (3.00 <= eta && eta < 3.10) * (48.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.440624) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (3.10 <= eta && eta < 3.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.000689) + 
    (3.10 <= eta && eta < 3.20) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.019759) + 
    (3.10 <= eta && eta < 3.20) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.040696) + 
    (3.10 <= eta && eta < 3.20) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.061575) + 
    (3.10 <= eta && eta < 3.20) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.082212) + 
    (3.10 <= eta && eta < 3.20) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.103793) + 
    (3.10 <= eta && eta < 3.20) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.125578) + 
    (3.10 <= eta && eta < 3.20) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.147013) + 
    (3.10 <= eta && eta < 3.20) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.167723) + 
    (3.10 <= eta && eta < 3.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.188656) + 
    (3.10 <= eta && eta < 3.20) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.210576) + 
    (3.10 <= eta && eta < 3.20) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.231900) + 
    (3.10 <= eta && eta < 3.20) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.252348) + 
    (3.10 <= eta && eta < 3.20) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.273410) + 
    (3.10 <= eta && eta < 3.20) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.295159) + 
    (3.10 <= eta && eta < 3.20) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.316734) + 
    (3.10 <= eta && eta < 3.20) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.337956) + 
    (3.10 <= eta && eta < 3.20) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.358943) + 
    (3.10 <= eta && eta < 3.20) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 35.60) * (0.380060) + 
    (3.10 <= eta && eta < 3.20) * (35.60 <= pt * cosh(eta) && pt * cosh(eta) < 39.70) * (0.401361) + 
    (3.10 <= eta && eta < 3.20) * (39.70 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.422890) + 
    (3.10 <= eta && eta < 3.20) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.438647) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (3.20 <= eta && eta < 3.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.000680) + 
    (3.20 <= eta && eta < 3.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.019670) + 
    (3.20 <= eta && eta < 3.30) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.040742) + 
    (3.20 <= eta && eta < 3.30) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.061790) + 
    (3.20 <= eta && eta < 3.30) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.082596) + 
    (3.20 <= eta && eta < 3.30) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.104346) + 
    (3.20 <= eta && eta < 3.30) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.126286) + 
    (3.20 <= eta && eta < 3.30) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.147856) + 
    (3.20 <= eta && eta < 3.30) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.168678) + 
    (3.20 <= eta && eta < 3.30) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.189706) + 
    (3.20 <= eta && eta < 3.30) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.211705) + 
    (3.20 <= eta && eta < 3.30) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.233084) + 
    (3.20 <= eta && eta < 3.30) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.253566) + 
    (3.20 <= eta && eta < 3.30) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.80) * (0.274643) + 
    (3.20 <= eta && eta < 3.30) * (24.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.296385) + 
    (3.20 <= eta && eta < 3.30) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.317933) + 
    (3.20 <= eta && eta < 3.30) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.339107) + 
    (3.20 <= eta && eta < 3.30) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.360440) + 
    (3.20 <= eta && eta < 3.30) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 35.60) * (0.381731) + 
    (3.20 <= eta && eta < 3.30) * (35.60 <= pt * cosh(eta) && pt * cosh(eta) < 39.80) * (0.402990) + 
    (3.20 <= eta && eta < 3.30) * (39.80 <= pt * cosh(eta) && pt * cosh(eta) < 46.00) * (0.424662) + 
    (3.20 <= eta && eta < 3.30) * (46.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.439913) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (3.30 <= eta && eta < 3.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.000662) + 
    (3.30 <= eta && eta < 3.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.019427) + 
    (3.30 <= eta && eta < 3.40) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.039437) + 
    (3.30 <= eta && eta < 3.40) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.059209) + 
    (3.30 <= eta && eta < 3.40) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.079967) + 
    (3.30 <= eta && eta < 3.40) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.101789) + 
    (3.30 <= eta && eta < 3.40) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.123881) + 
    (3.30 <= eta && eta < 3.40) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.145654) + 
    (3.30 <= eta && eta < 3.40) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.166700) + 
    (3.30 <= eta && eta < 3.40) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.187972) + 
    (3.30 <= eta && eta < 3.40) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.209116) + 
    (3.30 <= eta && eta < 3.40) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.229789) + 
    (3.30 <= eta && eta < 3.40) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.250700) + 
    (3.30 <= eta && eta < 3.40) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.271379) + 
    (3.30 <= eta && eta < 3.40) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.292170) + 
    (3.30 <= eta && eta < 3.40) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.313151) + 
    (3.30 <= eta && eta < 3.40) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.334116) + 
    (3.30 <= eta && eta < 3.40) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.355155) + 
    (3.30 <= eta && eta < 3.40) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.376266) + 
    (3.30 <= eta && eta < 3.40) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 38.40) * (0.397689) + 
    (3.30 <= eta && eta < 3.40) * (38.40 <= pt * cosh(eta) && pt * cosh(eta) < 44.00) * (0.419478) + 
    (3.30 <= eta && eta < 3.40) * (44.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.438009) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (3.40 <= eta && eta < 3.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.000723) + 
    (3.40 <= eta && eta < 3.50) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.019741) + 
    (3.40 <= eta && eta < 3.50) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.039030) + 
    (3.40 <= eta && eta < 3.50) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.058859) + 
    (3.40 <= eta && eta < 3.50) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.079716) + 
    (3.40 <= eta && eta < 3.50) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.101660) + 
    (3.40 <= eta && eta < 3.50) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.123883) + 
    (3.40 <= eta && eta < 3.50) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.145781) + 
    (3.40 <= eta && eta < 3.50) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.20) * (0.166943) + 
    (3.40 <= eta && eta < 3.50) * (19.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.188323) + 
    (3.40 <= eta && eta < 3.50) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.209563) + 
    (3.40 <= eta && eta < 3.50) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.230317) + 
    (3.40 <= eta && eta < 3.50) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.251296) + 
    (3.40 <= eta && eta < 3.50) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.272029) + 
    (3.40 <= eta && eta < 3.50) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.292858) + 
    (3.40 <= eta && eta < 3.50) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.313860) + 
    (3.40 <= eta && eta < 3.50) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.334830) + 
    (3.40 <= eta && eta < 3.50) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.355856) + 
    (3.40 <= eta && eta < 3.50) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.376934) + 
    (3.40 <= eta && eta < 3.50) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 38.30) * (0.398303) + 
    (3.40 <= eta && eta < 3.50) * (38.30 <= pt * cosh(eta) && pt * cosh(eta) < 43.90) * (0.420015) + 
    (3.40 <= eta && eta < 3.50) * (43.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.438585) } 

  add EfficiencyFormula {321} {2212} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.00 <= eta && eta < 1.10) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.000792) + 
    (1.00 <= eta && eta < 1.10) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.020173) + 
    (1.00 <= eta && eta < 1.10) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.040724) + 
    (1.00 <= eta && eta < 1.10) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.061613) + 
    (1.00 <= eta && eta < 1.10) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.082325) + 
    (1.00 <= eta && eta < 1.10) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 26.70) * (0.103121) + 
    (1.00 <= eta && eta < 1.10) * (26.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.124122) + 
    (1.00 <= eta && eta < 1.10) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.144827) + 
    (1.00 <= eta && eta < 1.10) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.165697) + 
    (1.00 <= eta && eta < 1.10) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.187160) + 
    (1.00 <= eta && eta < 1.10) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.30) * (0.208776) + 
    (1.00 <= eta && eta < 1.10) * (33.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.230151) + 
    (1.00 <= eta && eta < 1.10) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.60) * (0.250948) + 
    (1.00 <= eta && eta < 1.10) * (36.60 <= pt * cosh(eta) && pt * cosh(eta) < 38.50) * (0.271436) + 
    (1.00 <= eta && eta < 1.10) * (38.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.292204) + 
    (1.00 <= eta && eta < 1.10) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 43.20) * (0.313037) + 
    (1.00 <= eta && eta < 1.10) * (43.20 <= pt * cosh(eta) && pt * cosh(eta) < 46.20) * (0.333932) + 
    (1.00 <= eta && eta < 1.10) * (46.20 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.354921) + 
    (1.00 <= eta && eta < 1.10) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.365503) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.10 <= eta && eta < 1.20) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.000749) + 
    (1.10 <= eta && eta < 1.20) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.019212) + 
    (1.10 <= eta && eta < 1.20) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.038690) + 
    (1.10 <= eta && eta < 1.20) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.80) * (0.058609) + 
    (1.10 <= eta && eta < 1.20) * (24.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.079317) + 
    (1.10 <= eta && eta < 1.20) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.100267) + 
    (1.10 <= eta && eta < 1.20) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.120667) + 
    (1.10 <= eta && eta < 1.20) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.80) * (0.141681) + 
    (1.10 <= eta && eta < 1.20) * (29.80 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.162883) + 
    (1.10 <= eta && eta < 1.20) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.183892) + 
    (1.10 <= eta && eta < 1.20) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.205123) + 
    (1.10 <= eta && eta < 1.20) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 35.60) * (0.226190) + 
    (1.10 <= eta && eta < 1.20) * (35.60 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.247348) + 
    (1.10 <= eta && eta < 1.20) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.268698) + 
    (1.10 <= eta && eta < 1.20) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 41.60) * (0.289693) + 
    (1.10 <= eta && eta < 1.20) * (41.60 <= pt * cosh(eta) && pt * cosh(eta) < 44.20) * (0.310690) + 
    (1.10 <= eta && eta < 1.20) * (44.20 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.331732) + 
    (1.10 <= eta && eta < 1.20) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.350284) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.20 <= eta && eta < 1.30) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.000735) + 
    (1.20 <= eta && eta < 1.30) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.019251) + 
    (1.20 <= eta && eta < 1.30) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.039416) + 
    (1.20 <= eta && eta < 1.30) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.060388) + 
    (1.20 <= eta && eta < 1.30) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.90) * (0.081469) + 
    (1.20 <= eta && eta < 1.30) * (26.90 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.102762) + 
    (1.20 <= eta && eta < 1.30) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.124282) + 
    (1.20 <= eta && eta < 1.30) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.145489) + 
    (1.20 <= eta && eta < 1.30) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.166775) + 
    (1.20 <= eta && eta < 1.30) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.187791) + 
    (1.20 <= eta && eta < 1.30) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.208241) + 
    (1.20 <= eta && eta < 1.30) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.229170) + 
    (1.20 <= eta && eta < 1.30) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.250125) + 
    (1.20 <= eta && eta < 1.30) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 40.60) * (0.270672) + 
    (1.20 <= eta && eta < 1.30) * (40.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.90) * (0.291380) + 
    (1.20 <= eta && eta < 1.30) * (42.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.60) * (0.312442) + 
    (1.20 <= eta && eta < 1.30) * (45.60 <= pt * cosh(eta) && pt * cosh(eta) < 48.70) * (0.333395) + 
    (1.20 <= eta && eta < 1.30) * (48.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.347295) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.30 <= eta && eta < 1.40) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.000754) + 
    (1.30 <= eta && eta < 1.40) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.019346) + 
    (1.30 <= eta && eta < 1.40) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.039034) + 
    (1.30 <= eta && eta < 1.40) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.059442) + 
    (1.30 <= eta && eta < 1.40) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.079966) + 
    (1.30 <= eta && eta < 1.40) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.100735) + 
    (1.30 <= eta && eta < 1.40) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.121778) + 
    (1.30 <= eta && eta < 1.40) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.142573) + 
    (1.30 <= eta && eta < 1.40) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.163507) + 
    (1.30 <= eta && eta < 1.40) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.184950) + 
    (1.30 <= eta && eta < 1.40) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 35.80) * (0.206480) + 
    (1.30 <= eta && eta < 1.40) * (35.80 <= pt * cosh(eta) && pt * cosh(eta) < 37.50) * (0.227721) + 
    (1.30 <= eta && eta < 1.40) * (37.50 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.248910) + 
    (1.30 <= eta && eta < 1.40) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 41.50) * (0.270149) + 
    (1.30 <= eta && eta < 1.40) * (41.50 <= pt * cosh(eta) && pt * cosh(eta) < 43.90) * (0.291367) + 
    (1.30 <= eta && eta < 1.40) * (43.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.60) * (0.312385) + 
    (1.30 <= eta && eta < 1.40) * (46.60 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.333232) + 
    (1.30 <= eta && eta < 1.40) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.343952) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.40 <= eta && eta < 1.50) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.000753) + 
    (1.40 <= eta && eta < 1.50) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.019627) + 
    (1.40 <= eta && eta < 1.50) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.039570) + 
    (1.40 <= eta && eta < 1.50) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.059653) + 
    (1.40 <= eta && eta < 1.50) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.079793) + 
    (1.40 <= eta && eta < 1.50) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.100160) + 
    (1.40 <= eta && eta < 1.50) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.50) * (0.120805) + 
    (1.40 <= eta && eta < 1.50) * (30.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.141231) + 
    (1.40 <= eta && eta < 1.50) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.161827) + 
    (1.40 <= eta && eta < 1.50) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.182963) + 
    (1.40 <= eta && eta < 1.50) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.30) * (0.204230) + 
    (1.40 <= eta && eta < 1.50) * (36.30 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.225258) + 
    (1.40 <= eta && eta < 1.50) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.90) * (0.246285) + 
    (1.40 <= eta && eta < 1.50) * (39.90 <= pt * cosh(eta) && pt * cosh(eta) < 42.00) * (0.267412) + 
    (1.40 <= eta && eta < 1.50) * (42.00 <= pt * cosh(eta) && pt * cosh(eta) < 44.40) * (0.288573) + 
    (1.40 <= eta && eta < 1.50) * (44.40 <= pt * cosh(eta) && pt * cosh(eta) < 47.10) * (0.309587) + 
    (1.40 <= eta && eta < 1.50) * (47.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.329525) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.50 <= eta && eta < 1.60) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.70) * (0.000735) + 
    (1.50 <= eta && eta < 1.60) * (21.70 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.019189) + 
    (1.50 <= eta && eta < 1.60) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.039219) + 
    (1.50 <= eta && eta < 1.60) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.059595) + 
    (1.50 <= eta && eta < 1.60) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.40) * (0.080172) + 
    (1.50 <= eta && eta < 1.60) * (28.40 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.101005) + 
    (1.50 <= eta && eta < 1.60) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.121330) + 
    (1.50 <= eta && eta < 1.60) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 32.40) * (0.142193) + 
    (1.50 <= eta && eta < 1.60) * (32.40 <= pt * cosh(eta) && pt * cosh(eta) < 33.80) * (0.163178) + 
    (1.50 <= eta && eta < 1.60) * (33.80 <= pt * cosh(eta) && pt * cosh(eta) < 35.30) * (0.183926) + 
    (1.50 <= eta && eta < 1.60) * (35.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.204812) + 
    (1.50 <= eta && eta < 1.60) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.225483) + 
    (1.50 <= eta && eta < 1.60) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 40.50) * (0.246179) + 
    (1.50 <= eta && eta < 1.60) * (40.50 <= pt * cosh(eta) && pt * cosh(eta) < 42.60) * (0.267007) + 
    (1.50 <= eta && eta < 1.60) * (42.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.287906) + 
    (1.50 <= eta && eta < 1.60) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 47.70) * (0.308703) + 
    (1.50 <= eta && eta < 1.60) * (47.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.326541) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.60 <= eta && eta < 1.70) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.000751) + 
    (1.60 <= eta && eta < 1.70) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.019346) + 
    (1.60 <= eta && eta < 1.70) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.039161) + 
    (1.60 <= eta && eta < 1.70) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.059970) + 
    (1.60 <= eta && eta < 1.70) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.80) * (0.081053) + 
    (1.60 <= eta && eta < 1.70) * (28.80 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.101651) + 
    (1.60 <= eta && eta < 1.70) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.122495) + 
    (1.60 <= eta && eta < 1.70) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 32.90) * (0.143833) + 
    (1.60 <= eta && eta < 1.70) * (32.90 <= pt * cosh(eta) && pt * cosh(eta) < 34.30) * (0.164504) + 
    (1.60 <= eta && eta < 1.70) * (34.30 <= pt * cosh(eta) && pt * cosh(eta) < 35.80) * (0.184936) + 
    (1.60 <= eta && eta < 1.70) * (35.80 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.205509) + 
    (1.60 <= eta && eta < 1.70) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 39.20) * (0.226469) + 
    (1.60 <= eta && eta < 1.70) * (39.20 <= pt * cosh(eta) && pt * cosh(eta) < 41.10) * (0.247385) + 
    (1.60 <= eta && eta < 1.70) * (41.10 <= pt * cosh(eta) && pt * cosh(eta) < 43.30) * (0.268326) + 
    (1.60 <= eta && eta < 1.70) * (43.30 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.289271) + 
    (1.60 <= eta && eta < 1.70) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 48.50) * (0.310020) + 
    (1.60 <= eta && eta < 1.70) * (48.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.325235) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.70 <= eta && eta < 1.80) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.000751) + 
    (1.70 <= eta && eta < 1.80) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.019708) + 
    (1.70 <= eta && eta < 1.80) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.040014) + 
    (1.70 <= eta && eta < 1.80) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 27.70) * (0.060761) + 
    (1.70 <= eta && eta < 1.80) * (27.70 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.081705) + 
    (1.70 <= eta && eta < 1.80) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.50) * (0.102896) + 
    (1.70 <= eta && eta < 1.80) * (30.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.124309) + 
    (1.70 <= eta && eta < 1.80) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 33.30) * (0.145412) + 
    (1.70 <= eta && eta < 1.80) * (33.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.80) * (0.166545) + 
    (1.70 <= eta && eta < 1.80) * (34.80 <= pt * cosh(eta) && pt * cosh(eta) < 36.30) * (0.187371) + 
    (1.70 <= eta && eta < 1.80) * (36.30 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.208231) + 
    (1.70 <= eta && eta < 1.80) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.80) * (0.229385) + 
    (1.70 <= eta && eta < 1.80) * (39.80 <= pt * cosh(eta) && pt * cosh(eta) < 41.80) * (0.250379) + 
    (1.70 <= eta && eta < 1.80) * (41.80 <= pt * cosh(eta) && pt * cosh(eta) < 44.00) * (0.271317) + 
    (1.70 <= eta && eta < 1.80) * (44.00 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.292147) + 
    (1.70 <= eta && eta < 1.80) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 49.40) * (0.313073) + 
    (1.70 <= eta && eta < 1.80) * (49.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.325254) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.80 <= eta && eta < 1.90) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.000731) + 
    (1.80 <= eta && eta < 1.90) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.019328) + 
    (1.80 <= eta && eta < 1.90) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.039326) + 
    (1.80 <= eta && eta < 1.90) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.059822) + 
    (1.80 <= eta && eta < 1.90) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.080561) + 
    (1.80 <= eta && eta < 1.90) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.60) * (0.101587) + 
    (1.80 <= eta && eta < 1.90) * (30.60 <= pt * cosh(eta) && pt * cosh(eta) < 32.00) * (0.122870) + 
    (1.80 <= eta && eta < 1.90) * (32.00 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.143877) + 
    (1.80 <= eta && eta < 1.90) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.164942) + 
    (1.80 <= eta && eta < 1.90) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.186393) + 
    (1.80 <= eta && eta < 1.90) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.207825) + 
    (1.80 <= eta && eta < 1.90) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 40.00) * (0.228884) + 
    (1.80 <= eta && eta < 1.90) * (40.00 <= pt * cosh(eta) && pt * cosh(eta) < 42.00) * (0.249798) + 
    (1.80 <= eta && eta < 1.90) * (42.00 <= pt * cosh(eta) && pt * cosh(eta) < 44.20) * (0.270672) + 
    (1.80 <= eta && eta < 1.90) * (44.20 <= pt * cosh(eta) && pt * cosh(eta) < 46.70) * (0.291456) + 
    (1.80 <= eta && eta < 1.90) * (46.70 <= pt * cosh(eta) && pt * cosh(eta) < 49.60) * (0.312351) + 
    (1.80 <= eta && eta < 1.90) * (49.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.323866) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.90 <= eta && eta < 2.00) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.000768) + 
    (1.90 <= eta && eta < 2.00) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.80) * (0.019875) + 
    (1.90 <= eta && eta < 2.00) * (24.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.50) * (0.039975) + 
    (1.90 <= eta && eta < 2.00) * (26.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.060453) + 
    (1.90 <= eta && eta < 2.00) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.40) * (0.081116) + 
    (1.90 <= eta && eta < 2.00) * (29.40 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.102036) + 
    (1.90 <= eta && eta < 2.00) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.123197) + 
    (1.90 <= eta && eta < 2.00) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.144078) + 
    (1.90 <= eta && eta < 2.00) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.165017) + 
    (1.90 <= eta && eta < 2.00) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.186345) + 
    (1.90 <= eta && eta < 2.00) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.40) * (0.207662) + 
    (1.90 <= eta && eta < 2.00) * (38.40 <= pt * cosh(eta) && pt * cosh(eta) < 40.20) * (0.228619) + 
    (1.90 <= eta && eta < 2.00) * (40.20 <= pt * cosh(eta) && pt * cosh(eta) < 42.20) * (0.249444) + 
    (1.90 <= eta && eta < 2.00) * (42.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.40) * (0.270242) + 
    (1.90 <= eta && eta < 2.00) * (44.40 <= pt * cosh(eta) && pt * cosh(eta) < 46.90) * (0.290963) + 
    (1.90 <= eta && eta < 2.00) * (46.90 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.311813) + 
    (1.90 <= eta && eta < 2.00) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.322648) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.00 <= eta && eta < 2.10) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.000763) + 
    (2.00 <= eta && eta < 2.10) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.019751) + 
    (2.00 <= eta && eta < 2.10) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.039688) + 
    (2.00 <= eta && eta < 2.10) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.060014) + 
    (2.00 <= eta && eta < 2.10) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.080541) + 
    (2.00 <= eta && eta < 2.10) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.90) * (0.101343) + 
    (2.00 <= eta && eta < 2.10) * (30.90 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.122402) + 
    (2.00 <= eta && eta < 2.10) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.143200) + 
    (2.00 <= eta && eta < 2.10) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.164074) + 
    (2.00 <= eta && eta < 2.10) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.185352) + 
    (2.00 <= eta && eta < 2.10) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.50) * (0.206635) + 
    (2.00 <= eta && eta < 2.10) * (38.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.30) * (0.227574) + 
    (2.00 <= eta && eta < 2.10) * (40.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.248396) + 
    (2.00 <= eta && eta < 2.10) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.269205) + 
    (2.00 <= eta && eta < 2.10) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 47.00) * (0.289951) + 
    (2.00 <= eta && eta < 2.10) * (47.00 <= pt * cosh(eta) && pt * cosh(eta) < 49.90) * (0.310839) + 
    (2.00 <= eta && eta < 2.10) * (49.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.321032) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.10 <= eta && eta < 2.20) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.000727) + 
    (2.10 <= eta && eta < 2.20) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.019127) + 
    (2.10 <= eta && eta < 2.20) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.038725) + 
    (2.10 <= eta && eta < 2.20) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.058821) + 
    (2.10 <= eta && eta < 2.20) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.079190) + 
    (2.10 <= eta && eta < 2.20) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.90) * (0.099884) + 
    (2.10 <= eta && eta < 2.20) * (30.90 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.120876) + 
    (2.10 <= eta && eta < 2.20) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.141641) + 
    (2.10 <= eta && eta < 2.20) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.162508) + 
    (2.10 <= eta && eta < 2.20) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.183801) + 
    (2.10 <= eta && eta < 2.20) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.50) * (0.205121) + 
    (2.10 <= eta && eta < 2.20) * (38.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.30) * (0.226112) + 
    (2.10 <= eta && eta < 2.20) * (40.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.246999) + 
    (2.10 <= eta && eta < 2.20) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.267886) + 
    (2.10 <= eta && eta < 2.20) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 47.00) * (0.288721) + 
    (2.10 <= eta && eta < 2.20) * (47.00 <= pt * cosh(eta) && pt * cosh(eta) < 49.90) * (0.309707) + 
    (2.10 <= eta && eta < 2.20) * (49.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.319951) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.20 <= eta && eta < 2.30) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.000745) + 
    (2.20 <= eta && eta < 2.30) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.019403) + 
    (2.20 <= eta && eta < 2.30) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.70) * (0.039056) + 
    (2.20 <= eta && eta < 2.30) * (26.70 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.059148) + 
    (2.20 <= eta && eta < 2.30) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.60) * (0.079483) + 
    (2.20 <= eta && eta < 2.30) * (29.60 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.100128) + 
    (2.20 <= eta && eta < 2.30) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 32.40) * (0.121061) + 
    (2.20 <= eta && eta < 2.30) * (32.40 <= pt * cosh(eta) && pt * cosh(eta) < 33.80) * (0.141765) + 
    (2.20 <= eta && eta < 2.30) * (33.80 <= pt * cosh(eta) && pt * cosh(eta) < 35.30) * (0.162570) + 
    (2.20 <= eta && eta < 2.30) * (35.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.183803) + 
    (2.20 <= eta && eta < 2.30) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.205065) + 
    (2.20 <= eta && eta < 2.30) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 40.40) * (0.226004) + 
    (2.20 <= eta && eta < 2.30) * (40.40 <= pt * cosh(eta) && pt * cosh(eta) < 42.40) * (0.246846) + 
    (2.20 <= eta && eta < 2.30) * (42.40 <= pt * cosh(eta) && pt * cosh(eta) < 44.60) * (0.267693) + 
    (2.20 <= eta && eta < 2.30) * (44.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.10) * (0.288496) + 
    (2.20 <= eta && eta < 2.30) * (47.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.309457) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.30 <= eta && eta < 2.40) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.000768) + 
    (2.30 <= eta && eta < 2.40) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.019754) + 
    (2.30 <= eta && eta < 2.40) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.039501) + 
    (2.30 <= eta && eta < 2.40) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.059615) + 
    (2.30 <= eta && eta < 2.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.079935) + 
    (2.30 <= eta && eta < 2.40) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.100543) + 
    (2.30 <= eta && eta < 2.40) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.121426) + 
    (2.30 <= eta && eta < 2.40) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.90) * (0.142073) + 
    (2.30 <= eta && eta < 2.40) * (33.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.40) * (0.162818) + 
    (2.30 <= eta && eta < 2.40) * (35.40 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.183988) + 
    (2.30 <= eta && eta < 2.40) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 38.70) * (0.205189) + 
    (2.30 <= eta && eta < 2.40) * (38.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.50) * (0.226070) + 
    (2.30 <= eta && eta < 2.40) * (40.50 <= pt * cosh(eta) && pt * cosh(eta) < 42.50) * (0.246859) + 
    (2.30 <= eta && eta < 2.40) * (42.50 <= pt * cosh(eta) && pt * cosh(eta) < 44.70) * (0.267658) + 
    (2.30 <= eta && eta < 2.40) * (44.70 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.288418) + 
    (2.30 <= eta && eta < 2.40) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.308989) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.40 <= eta && eta < 2.50) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.000773) + 
    (2.40 <= eta && eta < 2.50) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.019842) + 
    (2.40 <= eta && eta < 2.50) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.039636) + 
    (2.40 <= eta && eta < 2.50) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.059781) + 
    (2.40 <= eta && eta < 2.50) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.080124) + 
    (2.40 <= eta && eta < 2.50) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.100746) + 
    (2.40 <= eta && eta < 2.50) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.121639) + 
    (2.40 <= eta && eta < 2.50) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.90) * (0.142291) + 
    (2.40 <= eta && eta < 2.50) * (33.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.40) * (0.163037) + 
    (2.40 <= eta && eta < 2.50) * (35.40 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.184205) + 
    (2.40 <= eta && eta < 2.50) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 38.70) * (0.205401) + 
    (2.40 <= eta && eta < 2.50) * (38.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.50) * (0.226275) + 
    (2.40 <= eta && eta < 2.50) * (40.50 <= pt * cosh(eta) && pt * cosh(eta) < 42.50) * (0.247054) + 
    (2.40 <= eta && eta < 2.50) * (42.50 <= pt * cosh(eta) && pt * cosh(eta) < 44.70) * (0.267842) + 
    (2.40 <= eta && eta < 2.50) * (44.70 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.288590) + 
    (2.40 <= eta && eta < 2.50) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.309148) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.50 <= eta && eta < 2.60) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.000746) + 
    (2.50 <= eta && eta < 2.60) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.019539) + 
    (2.50 <= eta && eta < 2.60) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.039556) + 
    (2.50 <= eta && eta < 2.60) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.060020) + 
    (2.50 <= eta && eta < 2.60) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.080706) + 
    (2.50 <= eta && eta < 2.60) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 30.70) * (0.101669) + 
    (2.50 <= eta && eta < 2.60) * (30.70 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.122884) + 
    (2.50 <= eta && eta < 2.60) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 33.50) * (0.143824) + 
    (2.50 <= eta && eta < 2.60) * (33.50 <= pt * cosh(eta) && pt * cosh(eta) < 35.00) * (0.164826) + 
    (2.50 <= eta && eta < 2.60) * (35.00 <= pt * cosh(eta) && pt * cosh(eta) < 36.60) * (0.186217) + 
    (2.50 <= eta && eta < 2.60) * (36.60 <= pt * cosh(eta) && pt * cosh(eta) < 38.30) * (0.207595) + 
    (2.50 <= eta && eta < 2.60) * (38.30 <= pt * cosh(eta) && pt * cosh(eta) < 40.10) * (0.228608) + 
    (2.50 <= eta && eta < 2.60) * (40.10 <= pt * cosh(eta) && pt * cosh(eta) < 42.10) * (0.249484) + 
    (2.50 <= eta && eta < 2.60) * (42.10 <= pt * cosh(eta) && pt * cosh(eta) < 44.30) * (0.270328) + 
    (2.50 <= eta && eta < 2.60) * (44.30 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.291089) + 
    (2.50 <= eta && eta < 2.60) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 49.70) * (0.311971) + 
    (2.50 <= eta && eta < 2.60) * (49.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.323152) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.60 <= eta && eta < 2.70) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.000749) + 
    (2.60 <= eta && eta < 2.70) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.019720) + 
    (2.60 <= eta && eta < 2.70) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.040133) + 
    (2.60 <= eta && eta < 2.70) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.60) * (0.060997) + 
    (2.60 <= eta && eta < 2.70) * (27.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.00) * (0.082052) + 
    (2.60 <= eta && eta < 2.70) * (29.00 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.103347) + 
    (2.60 <= eta && eta < 2.70) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.124852) + 
    (2.60 <= eta && eta < 2.70) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.146033) + 
    (2.60 <= eta && eta < 2.70) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.167231) + 
    (2.60 <= eta && eta < 2.70) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.188109) + 
    (2.60 <= eta && eta < 2.70) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 37.90) * (0.209008) + 
    (2.60 <= eta && eta < 2.70) * (37.90 <= pt * cosh(eta) && pt * cosh(eta) < 39.70) * (0.230189) + 
    (2.60 <= eta && eta < 2.70) * (39.70 <= pt * cosh(eta) && pt * cosh(eta) < 41.70) * (0.251197) + 
    (2.60 <= eta && eta < 2.70) * (41.70 <= pt * cosh(eta) && pt * cosh(eta) < 43.90) * (0.272136) + 
    (2.60 <= eta && eta < 2.70) * (43.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.40) * (0.292956) + 
    (2.60 <= eta && eta < 2.70) * (46.40 <= pt * cosh(eta) && pt * cosh(eta) < 49.30) * (0.313858) + 
    (2.60 <= eta && eta < 2.70) * (49.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.326348) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.70 <= eta && eta < 2.80) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.000734) + 
    (2.70 <= eta && eta < 2.80) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.019136) + 
    (2.70 <= eta && eta < 2.80) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.039036) + 
    (2.70 <= eta && eta < 2.80) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.059994) + 
    (2.70 <= eta && eta < 2.80) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.60) * (0.081246) + 
    (2.70 <= eta && eta < 2.80) * (28.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.90) * (0.102010) + 
    (2.70 <= eta && eta < 2.80) * (29.90 <= pt * cosh(eta) && pt * cosh(eta) < 31.20) * (0.122240) + 
    (2.70 <= eta && eta < 2.80) * (31.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.60) * (0.142986) + 
    (2.70 <= eta && eta < 2.80) * (32.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.163845) + 
    (2.70 <= eta && eta < 2.80) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 35.50) * (0.184464) + 
    (2.70 <= eta && eta < 2.80) * (35.50 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.205221) + 
    (2.70 <= eta && eta < 2.80) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 38.80) * (0.225768) + 
    (2.70 <= eta && eta < 2.80) * (38.80 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.246347) + 
    (2.70 <= eta && eta < 2.80) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 42.80) * (0.267066) + 
    (2.70 <= eta && eta < 2.80) * (42.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.20) * (0.287866) + 
    (2.70 <= eta && eta < 2.80) * (45.20 <= pt * cosh(eta) && pt * cosh(eta) < 48.00) * (0.308948) + 
    (2.70 <= eta && eta < 2.80) * (48.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.326040) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.80 <= eta && eta < 2.90) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.000755) + 
    (2.80 <= eta && eta < 2.90) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.019582) + 
    (2.80 <= eta && eta < 2.90) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.039299) + 
    (2.80 <= eta && eta < 2.90) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.059144) + 
    (2.80 <= eta && eta < 2.90) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.079825) + 
    (2.80 <= eta && eta < 2.90) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.100784) + 
    (2.80 <= eta && eta < 2.90) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.121242) + 
    (2.80 <= eta && eta < 2.90) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.142240) + 
    (2.80 <= eta && eta < 2.90) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.163357) + 
    (2.80 <= eta && eta < 2.90) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.184227) + 
    (2.80 <= eta && eta < 2.90) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.205225) + 
    (2.80 <= eta && eta < 2.90) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.40) * (0.225993) + 
    (2.80 <= eta && eta < 2.90) * (38.40 <= pt * cosh(eta) && pt * cosh(eta) < 40.30) * (0.246773) + 
    (2.80 <= eta && eta < 2.90) * (40.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.40) * (0.267668) + 
    (2.80 <= eta && eta < 2.90) * (42.40 <= pt * cosh(eta) && pt * cosh(eta) < 44.80) * (0.288618) + 
    (2.80 <= eta && eta < 2.90) * (44.80 <= pt * cosh(eta) && pt * cosh(eta) < 47.50) * (0.309448) + 
    (2.80 <= eta && eta < 2.90) * (47.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.327951) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.90 <= eta && eta < 3.00) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.000763) + 
    (2.90 <= eta && eta < 3.00) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.019795) + 
    (2.90 <= eta && eta < 3.00) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.039828) + 
    (2.90 <= eta && eta < 3.00) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.059973) + 
    (2.90 <= eta && eta < 3.00) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.080154) + 
    (2.90 <= eta && eta < 3.00) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.100548) + 
    (2.90 <= eta && eta < 3.00) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.50) * (0.121212) + 
    (2.90 <= eta && eta < 3.00) * (30.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.141647) + 
    (2.90 <= eta && eta < 3.00) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.162245) + 
    (2.90 <= eta && eta < 3.00) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.183377) + 
    (2.90 <= eta && eta < 3.00) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.30) * (0.204635) + 
    (2.90 <= eta && eta < 3.00) * (36.30 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.225649) + 
    (2.90 <= eta && eta < 3.00) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.90) * (0.246659) + 
    (2.90 <= eta && eta < 3.00) * (39.90 <= pt * cosh(eta) && pt * cosh(eta) < 42.00) * (0.267765) + 
    (2.90 <= eta && eta < 3.00) * (42.00 <= pt * cosh(eta) && pt * cosh(eta) < 44.40) * (0.288902) + 
    (2.90 <= eta && eta < 3.00) * (44.40 <= pt * cosh(eta) && pt * cosh(eta) < 47.10) * (0.309890) + 
    (2.90 <= eta && eta < 3.00) * (47.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.329801) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (3.00 <= eta && eta < 3.10) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.000757) + 
    (3.00 <= eta && eta < 3.10) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.40) * (0.019793) + 
    (3.00 <= eta && eta < 3.10) * (23.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.040032) + 
    (3.00 <= eta && eta < 3.10) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.060403) + 
    (3.00 <= eta && eta < 3.10) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.70) * (0.080804) + 
    (3.00 <= eta && eta < 3.10) * (27.70 <= pt * cosh(eta) && pt * cosh(eta) < 29.00) * (0.101403) + 
    (3.00 <= eta && eta < 3.10) * (29.00 <= pt * cosh(eta) && pt * cosh(eta) < 30.30) * (0.122252) + 
    (3.00 <= eta && eta < 3.10) * (30.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.142847) + 
    (3.00 <= eta && eta < 3.10) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 33.00) * (0.163582) + 
    (3.00 <= eta && eta < 3.10) * (33.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.184829) + 
    (3.00 <= eta && eta < 3.10) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 36.10) * (0.206176) + 
    (3.00 <= eta && eta < 3.10) * (36.10 <= pt * cosh(eta) && pt * cosh(eta) < 37.80) * (0.227253) + 
    (3.00 <= eta && eta < 3.10) * (37.80 <= pt * cosh(eta) && pt * cosh(eta) < 39.70) * (0.248298) + 
    (3.00 <= eta && eta < 3.10) * (39.70 <= pt * cosh(eta) && pt * cosh(eta) < 41.80) * (0.269416) + 
    (3.00 <= eta && eta < 3.10) * (41.80 <= pt * cosh(eta) && pt * cosh(eta) < 44.20) * (0.290537) + 
    (3.00 <= eta && eta < 3.10) * (44.20 <= pt * cosh(eta) && pt * cosh(eta) < 46.90) * (0.311485) + 
    (3.00 <= eta && eta < 3.10) * (46.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.331973) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (3.10 <= eta && eta < 3.20) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.000741) + 
    (3.10 <= eta && eta < 3.20) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.019121) + 
    (3.10 <= eta && eta < 3.20) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.038687) + 
    (3.10 <= eta && eta < 3.20) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.059011) + 
    (3.10 <= eta && eta < 3.20) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.079478) + 
    (3.10 <= eta && eta < 3.20) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.100208) + 
    (3.10 <= eta && eta < 3.20) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.121227) + 
    (3.10 <= eta && eta < 3.20) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.142009) + 
    (3.10 <= eta && eta < 3.20) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.162941) + 
    (3.10 <= eta && eta < 3.20) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.184390) + 
    (3.10 <= eta && eta < 3.20) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 35.80) * (0.205933) + 
    (3.10 <= eta && eta < 3.20) * (35.80 <= pt * cosh(eta) && pt * cosh(eta) < 37.50) * (0.227192) + 
    (3.10 <= eta && eta < 3.20) * (37.50 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.248406) + 
    (3.10 <= eta && eta < 3.20) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 41.50) * (0.269674) + 
    (3.10 <= eta && eta < 3.20) * (41.50 <= pt * cosh(eta) && pt * cosh(eta) < 43.80) * (0.290490) + 
    (3.10 <= eta && eta < 3.20) * (43.80 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.311207) + 
    (3.10 <= eta && eta < 3.20) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 49.70) * (0.332209) + 
    (3.10 <= eta && eta < 3.20) * (49.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.343307) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (3.20 <= eta && eta < 3.30) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.000768) + 
    (3.20 <= eta && eta < 3.30) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.019633) + 
    (3.20 <= eta && eta < 3.30) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.039580) + 
    (3.20 <= eta && eta < 3.30) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.060212) + 
    (3.20 <= eta && eta < 3.30) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.080924) + 
    (3.20 <= eta && eta < 3.30) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 28.60) * (0.101848) + 
    (3.20 <= eta && eta < 3.30) * (28.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.90) * (0.123017) + 
    (3.20 <= eta && eta < 3.30) * (29.90 <= pt * cosh(eta) && pt * cosh(eta) < 31.20) * (0.143910) + 
    (3.20 <= eta && eta < 3.30) * (31.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.60) * (0.164917) + 
    (3.20 <= eta && eta < 3.30) * (32.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.10) * (0.186410) + 
    (3.20 <= eta && eta < 3.30) * (34.10 <= pt * cosh(eta) && pt * cosh(eta) < 35.70) * (0.207967) + 
    (3.20 <= eta && eta < 3.30) * (35.70 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.229213) + 
    (3.20 <= eta && eta < 3.30) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 39.30) * (0.250388) + 
    (3.20 <= eta && eta < 3.30) * (39.30 <= pt * cosh(eta) && pt * cosh(eta) < 41.40) * (0.271594) + 
    (3.20 <= eta && eta < 3.30) * (41.40 <= pt * cosh(eta) && pt * cosh(eta) < 43.80) * (0.292761) + 
    (3.20 <= eta && eta < 3.30) * (43.80 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.313711) + 
    (3.20 <= eta && eta < 3.30) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 49.70) * (0.334474) + 
    (3.20 <= eta && eta < 3.30) * (49.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.345440) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (3.30 <= eta && eta < 3.40) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.000788) + 
    (3.30 <= eta && eta < 3.40) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.90) * (0.020005) + 
    (3.30 <= eta && eta < 3.40) * (22.90 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.040258) + 
    (3.30 <= eta && eta < 3.40) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.061145) + 
    (3.30 <= eta && eta < 3.40) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.082064) + 
    (3.30 <= eta && eta < 3.40) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.103157) + 
    (3.30 <= eta && eta < 3.40) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.80) * (0.124460) + 
    (3.30 <= eta && eta < 3.40) * (29.80 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.145453) + 
    (3.30 <= eta && eta < 3.40) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.166532) + 
    (3.30 <= eta && eta < 3.40) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.90) * (0.187357) + 
    (3.30 <= eta && eta < 3.40) * (33.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.40) * (0.207638) + 
    (3.30 <= eta && eta < 3.40) * (35.40 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.228413) + 
    (3.30 <= eta && eta < 3.40) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 38.90) * (0.249237) + 
    (3.30 <= eta && eta < 3.40) * (38.90 <= pt * cosh(eta) && pt * cosh(eta) < 40.90) * (0.269678) + 
    (3.30 <= eta && eta < 3.40) * (40.90 <= pt * cosh(eta) && pt * cosh(eta) < 43.20) * (0.290306) + 
    (3.30 <= eta && eta < 3.40) * (43.20 <= pt * cosh(eta) && pt * cosh(eta) < 45.90) * (0.311314) + 
    (3.30 <= eta && eta < 3.40) * (45.90 <= pt * cosh(eta) && pt * cosh(eta) < 49.00) * (0.332244) + 
    (3.30 <= eta && eta < 3.40) * (49.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.345265) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (3.40 <= eta && eta < 3.50) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.000743) + 
    (3.40 <= eta && eta < 3.50) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.019342) + 
    (3.40 <= eta && eta < 3.50) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.039448) + 
    (3.40 <= eta && eta < 3.50) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.060332) + 
    (3.40 <= eta && eta < 3.50) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.081319) + 
    (3.40 <= eta && eta < 3.50) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.102517) + 
    (3.40 <= eta && eta < 3.50) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.60) * (0.123948) + 
    (3.40 <= eta && eta < 3.50) * (29.60 <= pt * cosh(eta) && pt * cosh(eta) < 30.90) * (0.145075) + 
    (3.40 <= eta && eta < 3.50) * (30.90 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.166290) + 
    (3.40 <= eta && eta < 3.50) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.187247) + 
    (3.40 <= eta && eta < 3.50) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.207650) + 
    (3.40 <= eta && eta < 3.50) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.228542) + 
    (3.40 <= eta && eta < 3.50) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 38.70) * (0.249470) + 
    (3.40 <= eta && eta < 3.50) * (38.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.270002) + 
    (3.40 <= eta && eta < 3.50) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 43.00) * (0.290705) + 
    (3.40 <= eta && eta < 3.50) * (43.00 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.311775) + 
    (3.40 <= eta && eta < 3.50) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 48.80) * (0.332748) + 
    (3.40 <= eta && eta < 3.50) * (48.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.346376) } 

  # --- pions ---

  add EfficiencyFormula {-211} {321} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.00 <= eta && eta < 1.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.000723) + 
    (1.00 <= eta && eta < 1.10) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.020359) + 
    (1.00 <= eta && eta < 1.10) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.041268) + 
    (1.00 <= eta && eta < 1.10) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.062771) + 
    (1.00 <= eta && eta < 1.10) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.083795) + 
    (1.00 <= eta && eta < 1.10) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.104271) + 
    (1.00 <= eta && eta < 1.10) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.124913) + 
    (1.00 <= eta && eta < 1.10) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.146676) + 
    (1.00 <= eta && eta < 1.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.169069) + 
    (1.00 <= eta && eta < 1.10) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.190320) + 
    (1.00 <= eta && eta < 1.10) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.211461) + 
    (1.00 <= eta && eta < 1.10) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.233274) + 
    (1.00 <= eta && eta < 1.10) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.255188) + 
    (1.00 <= eta && eta < 1.10) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.276693) + 
    (1.00 <= eta && eta < 1.10) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.298134) + 
    (1.00 <= eta && eta < 1.10) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.319579) + 
    (1.00 <= eta && eta < 1.10) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.340806) + 
    (1.00 <= eta && eta < 1.10) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.361894) + 
    (1.00 <= eta && eta < 1.10) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.383172) + 
    (1.00 <= eta && eta < 1.10) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.404620) + 
    (1.00 <= eta && eta < 1.10) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 43.40) * (0.426353) + 
    (1.00 <= eta && eta < 1.10) * (43.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.444798) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.10 <= eta && eta < 1.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.000688) + 
    (1.10 <= eta && eta < 1.20) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.019438) + 
    (1.10 <= eta && eta < 1.20) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.039279) + 
    (1.10 <= eta && eta < 1.20) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.059801) + 
    (1.10 <= eta && eta < 1.20) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.081395) + 
    (1.10 <= eta && eta < 1.20) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.102638) + 
    (1.10 <= eta && eta < 1.20) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.124115) + 
    (1.10 <= eta && eta < 1.20) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.146671) + 
    (1.10 <= eta && eta < 1.20) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.168432) + 
    (1.10 <= eta && eta < 1.20) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.190363) + 
    (1.10 <= eta && eta < 1.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.212088) + 
    (1.10 <= eta && eta < 1.20) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.233248) + 
    (1.10 <= eta && eta < 1.20) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.254563) + 
    (1.10 <= eta && eta < 1.20) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.275549) + 
    (1.10 <= eta && eta < 1.20) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.296550) + 
    (1.10 <= eta && eta < 1.20) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.317640) + 
    (1.10 <= eta && eta < 1.20) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.338609) + 
    (1.10 <= eta && eta < 1.20) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.359540) + 
    (1.10 <= eta && eta < 1.20) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.380769) + 
    (1.10 <= eta && eta < 1.20) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.402289) + 
    (1.10 <= eta && eta < 1.20) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 43.90) * (0.424061) + 
    (1.10 <= eta && eta < 1.20) * (43.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.442161) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.20 <= eta && eta < 1.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.000680) + 
    (1.20 <= eta && eta < 1.30) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.019858) + 
    (1.20 <= eta && eta < 1.30) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.040284) + 
    (1.20 <= eta && eta < 1.30) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.060404) + 
    (1.20 <= eta && eta < 1.30) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.081468) + 
    (1.20 <= eta && eta < 1.30) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.103556) + 
    (1.20 <= eta && eta < 1.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.125865) + 
    (1.20 <= eta && eta < 1.30) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.147802) + 
    (1.20 <= eta && eta < 1.30) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.20) * (0.168967) + 
    (1.20 <= eta && eta < 1.30) * (19.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.190319) + 
    (1.20 <= eta && eta < 1.30) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.211508) + 
    (1.20 <= eta && eta < 1.30) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.232191) + 
    (1.20 <= eta && eta < 1.30) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.253082) + 
    (1.20 <= eta && eta < 1.30) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.274542) + 
    (1.20 <= eta && eta < 1.30) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.295913) + 
    (1.20 <= eta && eta < 1.30) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.317209) + 
    (1.20 <= eta && eta < 1.30) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.338823) + 
    (1.20 <= eta && eta < 1.30) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.360107) + 
    (1.20 <= eta && eta < 1.30) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.00) * (0.381426) + 
    (1.20 <= eta && eta < 1.30) * (35.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.10) * (0.402823) + 
    (1.20 <= eta && eta < 1.30) * (39.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.20) * (0.424492) + 
    (1.20 <= eta && eta < 1.30) * (45.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.440832) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.30 <= eta && eta < 1.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.000701) + 
    (1.30 <= eta && eta < 1.40) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.019989) + 
    (1.30 <= eta && eta < 1.40) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.041052) + 
    (1.30 <= eta && eta < 1.40) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.062014) + 
    (1.30 <= eta && eta < 1.40) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.082706) + 
    (1.30 <= eta && eta < 1.40) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.104326) + 
    (1.30 <= eta && eta < 1.40) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.126133) + 
    (1.30 <= eta && eta < 1.40) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.147579) + 
    (1.30 <= eta && eta < 1.40) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.168289) + 
    (1.30 <= eta && eta < 1.40) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.189215) + 
    (1.30 <= eta && eta < 1.40) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.211120) + 
    (1.30 <= eta && eta < 1.40) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.232423) + 
    (1.30 <= eta && eta < 1.40) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.252847) + 
    (1.30 <= eta && eta < 1.40) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.273879) + 
    (1.30 <= eta && eta < 1.40) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.295594) + 
    (1.30 <= eta && eta < 1.40) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.317131) + 
    (1.30 <= eta && eta < 1.40) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.338314) + 
    (1.30 <= eta && eta < 1.40) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.60) * (0.359675) + 
    (1.30 <= eta && eta < 1.40) * (32.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.70) * (0.381014) + 
    (1.30 <= eta && eta < 1.40) * (35.70 <= pt * cosh(eta) && pt * cosh(eta) < 39.90) * (0.402342) + 
    (1.30 <= eta && eta < 1.40) * (39.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.10) * (0.424107) + 
    (1.30 <= eta && eta < 1.40) * (46.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.439317) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.40 <= eta && eta < 1.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.000671) + 
    (1.40 <= eta && eta < 1.50) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.019241) + 
    (1.40 <= eta && eta < 1.50) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.039529) + 
    (1.40 <= eta && eta < 1.50) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.061058) + 
    (1.40 <= eta && eta < 1.50) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.082565) + 
    (1.40 <= eta && eta < 1.50) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.103768) + 
    (1.40 <= eta && eta < 1.50) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.125167) + 
    (1.40 <= eta && eta < 1.50) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.146237) + 
    (1.40 <= eta && eta < 1.50) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.167853) + 
    (1.40 <= eta && eta < 1.50) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.189600) + 
    (1.40 <= eta && eta < 1.50) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.211079) + 
    (1.40 <= eta && eta < 1.50) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.232986) + 
    (1.40 <= eta && eta < 1.50) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.254800) + 
    (1.40 <= eta && eta < 1.50) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.276071) + 
    (1.40 <= eta && eta < 1.50) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.297134) + 
    (1.40 <= eta && eta < 1.50) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.318082) + 
    (1.40 <= eta && eta < 1.50) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.339249) + 
    (1.40 <= eta && eta < 1.50) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.30) * (0.360510) + 
    (1.40 <= eta && eta < 1.50) * (33.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.381635) + 
    (1.40 <= eta && eta < 1.50) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.80) * (0.402956) + 
    (1.40 <= eta && eta < 1.50) * (40.80 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.424693) + 
    (1.40 <= eta && eta < 1.50) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.438435) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.50 <= eta && eta < 1.60) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.000673) + 
    (1.50 <= eta && eta < 1.60) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.019124) + 
    (1.50 <= eta && eta < 1.60) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.039002) + 
    (1.50 <= eta && eta < 1.60) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.060084) + 
    (1.50 <= eta && eta < 1.60) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.081176) + 
    (1.50 <= eta && eta < 1.60) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.102012) + 
    (1.50 <= eta && eta < 1.60) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.123085) + 
    (1.50 <= eta && eta < 1.60) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.143882) + 
    (1.50 <= eta && eta < 1.60) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.165267) + 
    (1.50 <= eta && eta < 1.60) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.186831) + 
    (1.50 <= eta && eta < 1.60) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.208178) + 
    (1.50 <= eta && eta < 1.60) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.230000) + 
    (1.50 <= eta && eta < 1.60) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.251777) + 
    (1.50 <= eta && eta < 1.60) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.273060) + 
    (1.50 <= eta && eta < 1.60) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.294179) + 
    (1.50 <= eta && eta < 1.60) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.90) * (0.315226) + 
    (1.50 <= eta && eta < 1.60) * (28.90 <= pt * cosh(eta) && pt * cosh(eta) < 30.90) * (0.336037) + 
    (1.50 <= eta && eta < 1.60) * (30.90 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.357142) + 
    (1.50 <= eta && eta < 1.60) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.378356) + 
    (1.50 <= eta && eta < 1.60) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.399705) + 
    (1.50 <= eta && eta < 1.60) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.421490) + 
    (1.50 <= eta && eta < 1.60) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.435963) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.60 <= eta && eta < 1.70) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.000707) + 
    (1.60 <= eta && eta < 1.70) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.019593) + 
    (1.60 <= eta && eta < 1.70) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.039389) + 
    (1.60 <= eta && eta < 1.70) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.060261) + 
    (1.60 <= eta && eta < 1.70) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.081099) + 
    (1.60 <= eta && eta < 1.70) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.101674) + 
    (1.60 <= eta && eta < 1.70) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.122490) + 
    (1.60 <= eta && eta < 1.70) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.143048) + 
    (1.60 <= eta && eta < 1.70) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.164209) + 
    (1.60 <= eta && eta < 1.70) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.185573) + 
    (1.60 <= eta && eta < 1.70) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.206751) + 
    (1.60 <= eta && eta < 1.70) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 23.40) * (0.228432) + 
    (1.60 <= eta && eta < 1.70) * (23.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.250102) + 
    (1.60 <= eta && eta < 1.70) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.271314) + 
    (1.60 <= eta && eta < 1.70) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.292398) + 
    (1.60 <= eta && eta < 1.70) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.313444) + 
    (1.60 <= eta && eta < 1.70) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.334290) + 
    (1.60 <= eta && eta < 1.70) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.355467) + 
    (1.60 <= eta && eta < 1.70) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.376793) + 
    (1.60 <= eta && eta < 1.70) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.80) * (0.398050) + 
    (1.60 <= eta && eta < 1.70) * (40.80 <= pt * cosh(eta) && pt * cosh(eta) < 46.70) * (0.419582) + 
    (1.60 <= eta && eta < 1.70) * (46.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.434240) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.70 <= eta && eta < 1.80) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.000692) + 
    (1.70 <= eta && eta < 1.80) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.020003) + 
    (1.70 <= eta && eta < 1.80) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.040686) + 
    (1.70 <= eta && eta < 1.80) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.061566) + 
    (1.70 <= eta && eta < 1.80) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.082303) + 
    (1.70 <= eta && eta < 1.80) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.102719) + 
    (1.70 <= eta && eta < 1.80) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.123343) + 
    (1.70 <= eta && eta < 1.80) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.144946) + 
    (1.70 <= eta && eta < 1.80) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.167051) + 
    (1.70 <= eta && eta < 1.80) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.188072) + 
    (1.70 <= eta && eta < 1.80) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.208898) + 
    (1.70 <= eta && eta < 1.80) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.230217) + 
    (1.70 <= eta && eta < 1.80) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.251534) + 
    (1.70 <= eta && eta < 1.80) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.272414) + 
    (1.70 <= eta && eta < 1.80) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 27.70) * (0.293188) + 
    (1.70 <= eta && eta < 1.80) * (27.70 <= pt * cosh(eta) && pt * cosh(eta) < 29.40) * (0.313951) + 
    (1.70 <= eta && eta < 1.80) * (29.40 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.335044) + 
    (1.70 <= eta && eta < 1.80) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.356349) + 
    (1.70 <= eta && eta < 1.80) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 37.20) * (0.377645) + 
    (1.70 <= eta && eta < 1.80) * (37.20 <= pt * cosh(eta) && pt * cosh(eta) < 41.40) * (0.399035) + 
    (1.70 <= eta && eta < 1.80) * (41.40 <= pt * cosh(eta) && pt * cosh(eta) < 47.50) * (0.420636) + 
    (1.70 <= eta && eta < 1.80) * (47.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.434156) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.80 <= eta && eta < 1.90) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.000706) + 
    (1.80 <= eta && eta < 1.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.020190) + 
    (1.80 <= eta && eta < 1.90) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.040805) + 
    (1.80 <= eta && eta < 1.90) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.061565) + 
    (1.80 <= eta && eta < 1.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.082168) + 
    (1.80 <= eta && eta < 1.90) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.102452) + 
    (1.80 <= eta && eta < 1.90) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.122948) + 
    (1.80 <= eta && eta < 1.90) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.144427) + 
    (1.80 <= eta && eta < 1.90) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.166418) + 
    (1.80 <= eta && eta < 1.90) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.70) * (0.187346) + 
    (1.80 <= eta && eta < 1.90) * (21.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.208094) + 
    (1.80 <= eta && eta < 1.90) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.229350) + 
    (1.80 <= eta && eta < 1.90) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.250620) + 
    (1.80 <= eta && eta < 1.90) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.271471) + 
    (1.80 <= eta && eta < 1.90) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.292233) + 
    (1.80 <= eta && eta < 1.90) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.313003) + 
    (1.80 <= eta && eta < 1.90) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.334119) + 
    (1.80 <= eta && eta < 1.90) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.10) * (0.355468) + 
    (1.80 <= eta && eta < 1.90) * (34.10 <= pt * cosh(eta) && pt * cosh(eta) < 37.30) * (0.376826) + 
    (1.80 <= eta && eta < 1.90) * (37.30 <= pt * cosh(eta) && pt * cosh(eta) < 41.50) * (0.398300) + 
    (1.80 <= eta && eta < 1.90) * (41.50 <= pt * cosh(eta) && pt * cosh(eta) < 47.50) * (0.419840) + 
    (1.80 <= eta && eta < 1.90) * (47.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.433334) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (1.90 <= eta && eta < 2.00) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.000726) + 
    (1.90 <= eta && eta < 2.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.020479) + 
    (1.90 <= eta && eta < 2.00) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.041079) + 
    (1.90 <= eta && eta < 2.00) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.061756) + 
    (1.90 <= eta && eta < 2.00) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.082251) + 
    (1.90 <= eta && eta < 2.00) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.102420) + 
    (1.90 <= eta && eta < 2.00) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.122800) + 
    (1.90 <= eta && eta < 2.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.144161) + 
    (1.90 <= eta && eta < 2.00) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.166041) + 
    (1.90 <= eta && eta < 2.00) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.186873) + 
    (1.90 <= eta && eta < 2.00) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.207538) + 
    (1.90 <= eta && eta < 2.00) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.228721) + 
    (1.90 <= eta && eta < 2.00) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.249934) + 
    (1.90 <= eta && eta < 2.00) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.270743) + 
    (1.90 <= eta && eta < 2.00) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.291479) + 
    (1.90 <= eta && eta < 2.00) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.60) * (0.312237) + 
    (1.90 <= eta && eta < 2.00) * (29.60 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.333360) + 
    (1.90 <= eta && eta < 2.00) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.354732) + 
    (1.90 <= eta && eta < 2.00) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 37.30) * (0.375809) + 
    (1.90 <= eta && eta < 2.00) * (37.30 <= pt * cosh(eta) && pt * cosh(eta) < 41.40) * (0.396913) + 
    (1.90 <= eta && eta < 2.00) * (41.40 <= pt * cosh(eta) && pt * cosh(eta) < 47.30) * (0.418399) + 
    (1.90 <= eta && eta < 2.00) * (47.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.432293) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.00 <= eta && eta < 2.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.000675) + 
    (2.00 <= eta && eta < 2.10) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.019521) + 
    (2.00 <= eta && eta < 2.10) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.039607) + 
    (2.00 <= eta && eta < 2.10) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.059943) + 
    (2.00 <= eta && eta < 2.10) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.081447) + 
    (2.00 <= eta && eta < 2.10) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.102754) + 
    (2.00 <= eta && eta < 2.10) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.20) * (0.123034) + 
    (2.00 <= eta && eta < 2.10) * (19.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.144287) + 
    (2.00 <= eta && eta < 2.10) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.166057) + 
    (2.00 <= eta && eta < 2.10) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.187897) + 
    (2.00 <= eta && eta < 2.10) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.209448) + 
    (2.00 <= eta && eta < 2.10) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.230382) + 
    (2.00 <= eta && eta < 2.10) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.251344) + 
    (2.00 <= eta && eta < 2.10) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 26.70) * (0.272676) + 
    (2.00 <= eta && eta < 2.10) * (26.70 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.293790) + 
    (2.00 <= eta && eta < 2.10) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.314730) + 
    (2.00 <= eta && eta < 2.10) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.335893) + 
    (2.00 <= eta && eta < 2.10) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.357125) + 
    (2.00 <= eta && eta < 2.10) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.378596) + 
    (2.00 <= eta && eta < 2.10) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 42.40) * (0.400170) + 
    (2.00 <= eta && eta < 2.10) * (42.40 <= pt * cosh(eta) && pt * cosh(eta) < 48.70) * (0.421818) + 
    (2.00 <= eta && eta < 2.10) * (48.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.433583) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.10 <= eta && eta < 2.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.000720) + 
    (2.10 <= eta && eta < 2.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.020280) + 
    (2.10 <= eta && eta < 2.20) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.040610) + 
    (2.10 <= eta && eta < 2.20) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.061035) + 
    (2.10 <= eta && eta < 2.20) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.082548) + 
    (2.10 <= eta && eta < 2.20) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.103815) + 
    (2.10 <= eta && eta < 2.20) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.124026) + 
    (2.10 <= eta && eta < 2.20) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.145187) + 
    (2.10 <= eta && eta < 2.20) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.166851) + 
    (2.10 <= eta && eta < 2.20) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.188577) + 
    (2.10 <= eta && eta < 2.20) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.210014) + 
    (2.10 <= eta && eta < 2.20) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.230838) + 
    (2.10 <= eta && eta < 2.20) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.251693) + 
    (2.10 <= eta && eta < 2.20) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.272923) + 
    (2.10 <= eta && eta < 2.20) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.40) * (0.294608) + 
    (2.10 <= eta && eta < 2.20) * (28.40 <= pt * cosh(eta) && pt * cosh(eta) < 30.20) * (0.315972) + 
    (2.10 <= eta && eta < 2.20) * (30.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.336881) + 
    (2.10 <= eta && eta < 2.20) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.357875) + 
    (2.10 <= eta && eta < 2.20) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.379129) + 
    (2.10 <= eta && eta < 2.20) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 42.60) * (0.400510) + 
    (2.10 <= eta && eta < 2.20) * (42.60 <= pt * cosh(eta) && pt * cosh(eta) < 49.00) * (0.422156) + 
    (2.10 <= eta && eta < 2.20) * (49.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.433561) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.20 <= eta && eta < 2.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.000691) + 
    (2.20 <= eta && eta < 2.30) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.019733) + 
    (2.20 <= eta && eta < 2.30) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.039771) + 
    (2.20 <= eta && eta < 2.30) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.060002) + 
    (2.20 <= eta && eta < 2.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.081377) + 
    (2.20 <= eta && eta < 2.30) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.102553) + 
    (2.20 <= eta && eta < 2.30) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.122711) + 
    (2.20 <= eta && eta < 2.30) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.143846) + 
    (2.20 <= eta && eta < 2.30) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.165506) + 
    (2.20 <= eta && eta < 2.30) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.187249) + 
    (2.20 <= eta && eta < 2.30) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.208720) + 
    (2.20 <= eta && eta < 2.30) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.229590) + 
    (2.20 <= eta && eta < 2.30) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.250504) + 
    (2.20 <= eta && eta < 2.30) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.271804) + 
    (2.20 <= eta && eta < 2.30) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.292903) + 
    (2.20 <= eta && eta < 2.30) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.313845) + 
    (2.20 <= eta && eta < 2.30) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.335028) + 
    (2.20 <= eta && eta < 2.30) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.80) * (0.356297) + 
    (2.20 <= eta && eta < 2.30) * (34.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.377513) + 
    (2.20 <= eta && eta < 2.30) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.398751) + 
    (2.20 <= eta && eta < 2.30) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 48.50) * (0.420375) + 
    (2.20 <= eta && eta < 2.30) * (48.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.432512) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.30 <= eta && eta < 2.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.000666) + 
    (2.30 <= eta && eta < 2.40) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.019270) + 
    (2.30 <= eta && eta < 2.40) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.039057) + 
    (2.30 <= eta && eta < 2.40) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.059121) + 
    (2.30 <= eta && eta < 2.40) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.080375) + 
    (2.30 <= eta && eta < 2.40) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.101472) + 
    (2.30 <= eta && eta < 2.40) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.121584) + 
    (2.30 <= eta && eta < 2.40) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.142694) + 
    (2.30 <= eta && eta < 2.40) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.164350) + 
    (2.30 <= eta && eta < 2.40) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.186107) + 
    (2.30 <= eta && eta < 2.40) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.207606) + 
    (2.30 <= eta && eta < 2.40) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.228515) + 
    (2.30 <= eta && eta < 2.40) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.249479) + 
    (2.30 <= eta && eta < 2.40) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.270839) + 
    (2.30 <= eta && eta < 2.40) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.292006) + 
    (2.30 <= eta && eta < 2.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.313021) + 
    (2.30 <= eta && eta < 2.40) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.334284) + 
    (2.30 <= eta && eta < 2.40) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.355239) + 
    (2.30 <= eta && eta < 2.40) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 37.90) * (0.376286) + 
    (2.30 <= eta && eta < 2.40) * (37.90 <= pt * cosh(eta) && pt * cosh(eta) < 42.10) * (0.397538) + 
    (2.30 <= eta && eta < 2.40) * (42.10 <= pt * cosh(eta) && pt * cosh(eta) < 48.20) * (0.419128) + 
    (2.30 <= eta && eta < 2.40) * (48.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.431773) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.40 <= eta && eta < 2.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.000671) + 
    (2.40 <= eta && eta < 2.50) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.019357) + 
    (2.40 <= eta && eta < 2.50) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.039192) + 
    (2.40 <= eta && eta < 2.50) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.059287) + 
    (2.40 <= eta && eta < 2.50) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.080564) + 
    (2.40 <= eta && eta < 2.50) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.101676) + 
    (2.40 <= eta && eta < 2.50) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.121797) + 
    (2.40 <= eta && eta < 2.50) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.142912) + 
    (2.40 <= eta && eta < 2.50) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.164569) + 
    (2.40 <= eta && eta < 2.50) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.186323) + 
    (2.40 <= eta && eta < 2.50) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.207817) + 
    (2.40 <= eta && eta < 2.50) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.228719) + 
    (2.40 <= eta && eta < 2.50) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.249673) + 
    (2.40 <= eta && eta < 2.50) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.271022) + 
    (2.40 <= eta && eta < 2.50) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.292176) + 
    (2.40 <= eta && eta < 2.50) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.313177) + 
    (2.40 <= eta && eta < 2.50) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.334425) + 
    (2.40 <= eta && eta < 2.50) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.80) * (0.355764) + 
    (2.40 <= eta && eta < 2.50) * (34.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.377053) + 
    (2.40 <= eta && eta < 2.50) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.398367) + 
    (2.40 <= eta && eta < 2.50) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 48.50) * (0.420069) + 
    (2.40 <= eta && eta < 2.50) * (48.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.432252) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.50 <= eta && eta < 2.60) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.000673) + 
    (2.50 <= eta && eta < 2.60) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.019568) + 
    (2.50 <= eta && eta < 2.60) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.039845) + 
    (2.50 <= eta && eta < 2.60) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.060381) + 
    (2.50 <= eta && eta < 2.60) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.080832) + 
    (2.50 <= eta && eta < 2.60) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.101016) + 
    (2.50 <= eta && eta < 2.60) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.121450) + 
    (2.50 <= eta && eta < 2.60) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.142896) + 
    (2.50 <= eta && eta < 2.60) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.164882) + 
    (2.50 <= eta && eta < 2.60) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.70) * (0.185828) + 
    (2.50 <= eta && eta < 2.60) * (21.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.206612) + 
    (2.50 <= eta && eta < 2.60) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.227920) + 
    (2.50 <= eta && eta < 2.60) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.249258) + 
    (2.50 <= eta && eta < 2.60) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.270187) + 
    (2.50 <= eta && eta < 2.60) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.291037) + 
    (2.50 <= eta && eta < 2.60) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.311903) + 
    (2.50 <= eta && eta < 2.60) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.332624) + 
    (2.50 <= eta && eta < 2.60) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.353737) + 
    (2.50 <= eta && eta < 2.60) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.375063) + 
    (2.50 <= eta && eta < 2.60) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 41.20) * (0.396391) + 
    (2.50 <= eta && eta < 2.60) * (41.20 <= pt * cosh(eta) && pt * cosh(eta) < 47.00) * (0.417902) + 
    (2.50 <= eta && eta < 2.60) * (47.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.432206) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.60 <= eta && eta < 2.70) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.000739) + 
    (2.60 <= eta && eta < 2.70) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.020090) + 
    (2.60 <= eta && eta < 2.70) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.039984) + 
    (2.60 <= eta && eta < 2.70) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.060848) + 
    (2.60 <= eta && eta < 2.70) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.081627) + 
    (2.60 <= eta && eta < 2.70) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.102117) + 
    (2.60 <= eta && eta < 2.70) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.122833) + 
    (2.60 <= eta && eta < 2.70) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.143287) + 
    (2.60 <= eta && eta < 2.70) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.164340) + 
    (2.60 <= eta && eta < 2.70) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.185599) + 
    (2.60 <= eta && eta < 2.70) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.206678) + 
    (2.60 <= eta && eta < 2.70) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.228265) + 
    (2.60 <= eta && eta < 2.70) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.249854) + 
    (2.60 <= eta && eta < 2.70) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.270996) + 
    (2.60 <= eta && eta < 2.70) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.292023) + 
    (2.60 <= eta && eta < 2.70) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.313027) + 
    (2.60 <= eta && eta < 2.70) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 31.20) * (0.333844) + 
    (2.60 <= eta && eta < 2.70) * (31.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.355009) + 
    (2.60 <= eta && eta < 2.70) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.376338) + 
    (2.60 <= eta && eta < 2.70) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 40.90) * (0.397617) + 
    (2.60 <= eta && eta < 2.70) * (40.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.419190) + 
    (2.60 <= eta && eta < 2.70) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.433757) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.70 <= eta && eta < 2.80) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.000706) + 
    (2.70 <= eta && eta < 2.80) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.019651) + 
    (2.70 <= eta && eta < 2.80) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.039646) + 
    (2.70 <= eta && eta < 2.80) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.060732) + 
    (2.70 <= eta && eta < 2.80) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.081770) + 
    (2.70 <= eta && eta < 2.80) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 17.80) * (0.102523) + 
    (2.70 <= eta && eta < 2.80) * (17.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.123497) + 
    (2.70 <= eta && eta < 2.80) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.144188) + 
    (2.70 <= eta && eta < 2.80) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.30) * (0.165463) + 
    (2.70 <= eta && eta < 2.80) * (20.30 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.186919) + 
    (2.70 <= eta && eta < 2.80) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.208163) + 
    (2.70 <= eta && eta < 2.80) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.30) * (0.229888) + 
    (2.70 <= eta && eta < 2.80) * (23.30 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.251580) + 
    (2.70 <= eta && eta < 2.80) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.272789) + 
    (2.70 <= eta && eta < 2.80) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.293848) + 
    (2.70 <= eta && eta < 2.80) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.00) * (0.314849) + 
    (2.70 <= eta && eta < 2.80) * (29.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.335628) + 
    (2.70 <= eta && eta < 2.80) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 33.50) * (0.356716) + 
    (2.70 <= eta && eta < 2.80) * (33.50 <= pt * cosh(eta) && pt * cosh(eta) < 36.60) * (0.377930) + 
    (2.70 <= eta && eta < 2.80) * (36.60 <= pt * cosh(eta) && pt * cosh(eta) < 40.80) * (0.399297) + 
    (2.70 <= eta && eta < 2.80) * (40.80 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.420954) + 
    (2.70 <= eta && eta < 2.80) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.435367) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.80 <= eta && eta < 2.90) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.000742) + 
    (2.80 <= eta && eta < 2.90) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.020385) + 
    (2.80 <= eta && eta < 2.90) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.040943) + 
    (2.80 <= eta && eta < 2.90) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.062496) + 
    (2.80 <= eta && eta < 2.90) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.083904) + 
    (2.80 <= eta && eta < 2.90) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.104944) + 
    (2.80 <= eta && eta < 2.90) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.126142) + 
    (2.80 <= eta && eta < 2.90) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.146998) + 
    (2.80 <= eta && eta < 2.90) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.168389) + 
    (2.80 <= eta && eta < 2.90) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.189914) + 
    (2.80 <= eta && eta < 2.90) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.211182) + 
    (2.80 <= eta && eta < 2.90) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.232890) + 
    (2.80 <= eta && eta < 2.90) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.254524) + 
    (2.80 <= eta && eta < 2.90) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.275642) + 
    (2.80 <= eta && eta < 2.90) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.296577) + 
    (2.80 <= eta && eta < 2.90) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.90) * (0.317424) + 
    (2.80 <= eta && eta < 2.90) * (28.90 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.338518) + 
    (2.80 <= eta && eta < 2.90) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 33.50) * (0.359737) + 
    (2.80 <= eta && eta < 2.90) * (33.50 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.380854) + 
    (2.80 <= eta && eta < 2.90) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 41.00) * (0.402203) + 
    (2.80 <= eta && eta < 2.90) * (41.00 <= pt * cosh(eta) && pt * cosh(eta) < 47.30) * (0.423852) + 
    (2.80 <= eta && eta < 2.90) * (47.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.437441) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (2.90 <= eta && eta < 3.00) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.000680) + 
    (2.90 <= eta && eta < 3.00) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.019408) + 
    (2.90 <= eta && eta < 3.00) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.039788) + 
    (2.90 <= eta && eta < 3.00) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.061380) + 
    (2.90 <= eta && eta < 3.00) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.082930) + 
    (2.90 <= eta && eta < 3.00) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.104161) + 
    (2.90 <= eta && eta < 3.00) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.125576) + 
    (2.90 <= eta && eta < 3.00) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.146654) + 
    (2.90 <= eta && eta < 3.00) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.168271) + 
    (2.90 <= eta && eta < 3.00) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.190012) + 
    (2.90 <= eta && eta < 3.00) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.211480) + 
    (2.90 <= eta && eta < 3.00) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.233371) + 
    (2.90 <= eta && eta < 3.00) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.255165) + 
    (2.90 <= eta && eta < 3.00) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.276414) + 
    (2.90 <= eta && eta < 3.00) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.297452) + 
    (2.90 <= eta && eta < 3.00) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.318373) + 
    (2.90 <= eta && eta < 3.00) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.339511) + 
    (2.90 <= eta && eta < 3.00) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.30) * (0.360741) + 
    (2.90 <= eta && eta < 3.00) * (33.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.381833) + 
    (2.90 <= eta && eta < 3.00) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.80) * (0.403120) + 
    (2.90 <= eta && eta < 3.00) * (40.80 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.424822) + 
    (2.90 <= eta && eta < 3.00) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.438540) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (3.00 <= eta && eta < 3.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.000689) + 
    (3.00 <= eta && eta < 3.10) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.019678) + 
    (3.00 <= eta && eta < 3.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.040386) + 
    (3.00 <= eta && eta < 3.10) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.062288) + 
    (3.00 <= eta && eta < 3.10) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.084104) + 
    (3.00 <= eta && eta < 3.10) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.105555) + 
    (3.00 <= eta && eta < 3.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.127155) + 
    (3.00 <= eta && eta < 3.10) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.148379) + 
    (3.00 <= eta && eta < 3.10) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.170110) + 
    (3.00 <= eta && eta < 3.10) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.191935) + 
    (3.00 <= eta && eta < 3.10) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.213453) + 
    (3.00 <= eta && eta < 3.10) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.90) * (0.235364) + 
    (3.00 <= eta && eta < 3.10) * (22.90 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.257149) + 
    (3.00 <= eta && eta < 3.10) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.278361) + 
    (3.00 <= eta && eta < 3.10) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.90) * (0.299337) + 
    (3.00 <= eta && eta < 3.10) * (26.90 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.320759) + 
    (3.00 <= eta && eta < 3.10) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.342209) + 
    (3.00 <= eta && eta < 3.10) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.363512) + 
    (3.00 <= eta && eta < 3.10) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.384827) + 
    (3.00 <= eta && eta < 3.10) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 41.20) * (0.406224) + 
    (3.00 <= eta && eta < 3.10) * (41.20 <= pt * cosh(eta) && pt * cosh(eta) < 48.00) * (0.428061) + 
    (3.00 <= eta && eta < 3.10) * (48.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.440624) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (3.10 <= eta && eta < 3.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.000689) + 
    (3.10 <= eta && eta < 3.20) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.019759) + 
    (3.10 <= eta && eta < 3.20) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.040696) + 
    (3.10 <= eta && eta < 3.20) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.061575) + 
    (3.10 <= eta && eta < 3.20) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.082212) + 
    (3.10 <= eta && eta < 3.20) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.103793) + 
    (3.10 <= eta && eta < 3.20) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.125578) + 
    (3.10 <= eta && eta < 3.20) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.147013) + 
    (3.10 <= eta && eta < 3.20) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.167723) + 
    (3.10 <= eta && eta < 3.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.188656) + 
    (3.10 <= eta && eta < 3.20) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.210576) + 
    (3.10 <= eta && eta < 3.20) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.231900) + 
    (3.10 <= eta && eta < 3.20) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.252348) + 
    (3.10 <= eta && eta < 3.20) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.273410) + 
    (3.10 <= eta && eta < 3.20) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.295159) + 
    (3.10 <= eta && eta < 3.20) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.316734) + 
    (3.10 <= eta && eta < 3.20) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.337956) + 
    (3.10 <= eta && eta < 3.20) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.358943) + 
    (3.10 <= eta && eta < 3.20) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 35.60) * (0.380060) + 
    (3.10 <= eta && eta < 3.20) * (35.60 <= pt * cosh(eta) && pt * cosh(eta) < 39.70) * (0.401361) + 
    (3.10 <= eta && eta < 3.20) * (39.70 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.422890) + 
    (3.10 <= eta && eta < 3.20) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.438647) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (3.20 <= eta && eta < 3.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.000680) + 
    (3.20 <= eta && eta < 3.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.019670) + 
    (3.20 <= eta && eta < 3.30) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.040742) + 
    (3.20 <= eta && eta < 3.30) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.061790) + 
    (3.20 <= eta && eta < 3.30) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.082596) + 
    (3.20 <= eta && eta < 3.30) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.104346) + 
    (3.20 <= eta && eta < 3.30) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.126286) + 
    (3.20 <= eta && eta < 3.30) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.147856) + 
    (3.20 <= eta && eta < 3.30) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.168678) + 
    (3.20 <= eta && eta < 3.30) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.189706) + 
    (3.20 <= eta && eta < 3.30) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.211705) + 
    (3.20 <= eta && eta < 3.30) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.233084) + 
    (3.20 <= eta && eta < 3.30) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.253566) + 
    (3.20 <= eta && eta < 3.30) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.80) * (0.274643) + 
    (3.20 <= eta && eta < 3.30) * (24.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.296385) + 
    (3.20 <= eta && eta < 3.30) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.317933) + 
    (3.20 <= eta && eta < 3.30) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.339107) + 
    (3.20 <= eta && eta < 3.30) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.360440) + 
    (3.20 <= eta && eta < 3.30) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 35.60) * (0.381731) + 
    (3.20 <= eta && eta < 3.30) * (35.60 <= pt * cosh(eta) && pt * cosh(eta) < 39.80) * (0.402990) + 
    (3.20 <= eta && eta < 3.30) * (39.80 <= pt * cosh(eta) && pt * cosh(eta) < 46.00) * (0.424662) + 
    (3.20 <= eta && eta < 3.30) * (46.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.439913) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (3.30 <= eta && eta < 3.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.000662) + 
    (3.30 <= eta && eta < 3.40) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.019427) + 
    (3.30 <= eta && eta < 3.40) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.039437) + 
    (3.30 <= eta && eta < 3.40) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.059209) + 
    (3.30 <= eta && eta < 3.40) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.079967) + 
    (3.30 <= eta && eta < 3.40) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.101789) + 
    (3.30 <= eta && eta < 3.40) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.123881) + 
    (3.30 <= eta && eta < 3.40) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.145654) + 
    (3.30 <= eta && eta < 3.40) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.166700) + 
    (3.30 <= eta && eta < 3.40) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.187972) + 
    (3.30 <= eta && eta < 3.40) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.209116) + 
    (3.30 <= eta && eta < 3.40) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.229789) + 
    (3.30 <= eta && eta < 3.40) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.250700) + 
    (3.30 <= eta && eta < 3.40) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.271379) + 
    (3.30 <= eta && eta < 3.40) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.292170) + 
    (3.30 <= eta && eta < 3.40) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.313151) + 
    (3.30 <= eta && eta < 3.40) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.334116) + 
    (3.30 <= eta && eta < 3.40) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.355155) + 
    (3.30 <= eta && eta < 3.40) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.376266) + 
    (3.30 <= eta && eta < 3.40) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 38.40) * (0.397689) + 
    (3.30 <= eta && eta < 3.40) * (38.40 <= pt * cosh(eta) && pt * cosh(eta) < 44.00) * (0.419478) + 
    (3.30 <= eta && eta < 3.40) * (44.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.438009) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.500000) + 
    (3.40 <= eta && eta < 3.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.000723) + 
    (3.40 <= eta && eta < 3.50) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.019741) + 
    (3.40 <= eta && eta < 3.50) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.039030) + 
    (3.40 <= eta && eta < 3.50) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.058859) + 
    (3.40 <= eta && eta < 3.50) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.079716) + 
    (3.40 <= eta && eta < 3.50) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.101660) + 
    (3.40 <= eta && eta < 3.50) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.123883) + 
    (3.40 <= eta && eta < 3.50) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.145781) + 
    (3.40 <= eta && eta < 3.50) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.20) * (0.166943) + 
    (3.40 <= eta && eta < 3.50) * (19.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.188323) + 
    (3.40 <= eta && eta < 3.50) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.209563) + 
    (3.40 <= eta && eta < 3.50) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.230317) + 
    (3.40 <= eta && eta < 3.50) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.251296) + 
    (3.40 <= eta && eta < 3.50) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.272029) + 
    (3.40 <= eta && eta < 3.50) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.292858) + 
    (3.40 <= eta && eta < 3.50) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.313860) + 
    (3.40 <= eta && eta < 3.50) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.334830) + 
    (3.40 <= eta && eta < 3.50) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.355856) + 
    (3.40 <= eta && eta < 3.50) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.376934) + 
    (3.40 <= eta && eta < 3.50) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 38.30) * (0.398303) + 
    (3.40 <= eta && eta < 3.50) * (38.30 <= pt * cosh(eta) && pt * cosh(eta) < 43.90) * (0.420015) + 
    (3.40 <= eta && eta < 3.50) * (43.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.438585) } 

  add EfficiencyFormula {-211} {-211} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.00 <= eta && eta < 1.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 10.70) * (0.999421) + 
    (1.00 <= eta && eta < 1.10) * (10.70 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.980545) + 
    (1.00 <= eta && eta < 1.10) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.959608) + 
    (1.00 <= eta && eta < 1.10) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.937259) + 
    (1.00 <= eta && eta < 1.10) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.915885) + 
    (1.00 <= eta && eta < 1.10) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.895401) + 
    (1.00 <= eta && eta < 1.10) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.874634) + 
    (1.00 <= eta && eta < 1.10) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.854132) + 
    (1.00 <= eta && eta < 1.10) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.832382) + 
    (1.00 <= eta && eta < 1.10) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.809991) + 
    (1.00 <= eta && eta < 1.10) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.787646) + 
    (1.00 <= eta && eta < 1.10) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.765900) + 
    (1.00 <= eta && eta < 1.10) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.745301) + 
    (1.00 <= eta && eta < 1.10) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.724020) + 
    (1.00 <= eta && eta < 1.10) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.702466) + 
    (1.00 <= eta && eta < 1.10) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.680941) + 
    (1.00 <= eta && eta < 1.10) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.658802) + 
    (1.00 <= eta && eta < 1.10) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.635344) + 
    (1.00 <= eta && eta < 1.10) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610241) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.10 <= eta && eta < 1.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.00) * (0.999432) + 
    (1.10 <= eta && eta < 1.20) * (11.00 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.981085) + 
    (1.10 <= eta && eta < 1.20) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.961106) + 
    (1.10 <= eta && eta < 1.20) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.939739) + 
    (1.10 <= eta && eta < 1.20) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.917253) + 
    (1.10 <= eta && eta < 1.20) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.895414) + 
    (1.10 <= eta && eta < 1.20) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.875248) + 
    (1.10 <= eta && eta < 1.20) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.853363) + 
    (1.10 <= eta && eta < 1.20) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.830383) + 
    (1.10 <= eta && eta < 1.20) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.808750) + 
    (1.10 <= eta && eta < 1.20) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.787127) + 
    (1.10 <= eta && eta < 1.20) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.764636) + 
    (1.10 <= eta && eta < 1.20) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.742284) + 
    (1.10 <= eta && eta < 1.20) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.719967) + 
    (1.10 <= eta && eta < 1.20) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.698081) + 
    (1.10 <= eta && eta < 1.20) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.676291) + 
    (1.10 <= eta && eta < 1.20) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.654085) + 
    (1.10 <= eta && eta < 1.20) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.629900) + 
    (1.10 <= eta && eta < 1.20) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609371) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.20 <= eta && eta < 1.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.30) * (0.999418) + 
    (1.20 <= eta && eta < 1.30) * (11.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.979968) + 
    (1.20 <= eta && eta < 1.30) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.958724) + 
    (1.20 <= eta && eta < 1.30) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.937510) + 
    (1.20 <= eta && eta < 1.30) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.915434) + 
    (1.20 <= eta && eta < 1.30) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.892143) + 
    (1.20 <= eta && eta < 1.30) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.868604) + 
    (1.20 <= eta && eta < 1.30) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.845570) + 
    (1.20 <= eta && eta < 1.30) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.823568) + 
    (1.20 <= eta && eta < 1.30) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.801299) + 
    (1.20 <= eta && eta < 1.30) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.779332) + 
    (1.20 <= eta && eta < 1.30) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.758243) + 
    (1.20 <= eta && eta < 1.30) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.70) * (0.737276) + 
    (1.20 <= eta && eta < 1.30) * (19.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.716295) + 
    (1.20 <= eta && eta < 1.30) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.694892) + 
    (1.20 <= eta && eta < 1.30) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.673234) + 
    (1.20 <= eta && eta < 1.30) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.651069) + 
    (1.20 <= eta && eta < 1.30) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.626031) + 
    (1.20 <= eta && eta < 1.30) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608708) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.30 <= eta && eta < 1.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.999375) + 
    (1.30 <= eta && eta < 1.40) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.979302) + 
    (1.30 <= eta && eta < 1.40) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.958306) + 
    (1.30 <= eta && eta < 1.40) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.937530) + 
    (1.30 <= eta && eta < 1.40) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.915968) + 
    (1.30 <= eta && eta < 1.40) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.893222) + 
    (1.30 <= eta && eta < 1.40) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.870202) + 
    (1.30 <= eta && eta < 1.40) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.847622) + 
    (1.30 <= eta && eta < 1.40) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.825993) + 
    (1.30 <= eta && eta < 1.40) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.804027) + 
    (1.30 <= eta && eta < 1.40) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.782280) + 
    (1.30 <= eta && eta < 1.40) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.761316) + 
    (1.30 <= eta && eta < 1.40) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.740385) + 
    (1.30 <= eta && eta < 1.40) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.719339) + 
    (1.30 <= eta && eta < 1.40) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.697756) + 
    (1.30 <= eta && eta < 1.40) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.675783) + 
    (1.30 <= eta && eta < 1.40) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.653479) + 
    (1.30 <= eta && eta < 1.40) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.629191) + 
    (1.30 <= eta && eta < 1.40) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609919) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.40 <= eta && eta < 1.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.999393) + 
    (1.40 <= eta && eta < 1.50) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.979903) + 
    (1.40 <= eta && eta < 1.50) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.959632) + 
    (1.40 <= eta && eta < 1.50) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.939507) + 
    (1.40 <= eta && eta < 1.50) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.918534) + 
    (1.40 <= eta && eta < 1.50) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.896310) + 
    (1.40 <= eta && eta < 1.50) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.873716) + 
    (1.40 <= eta && eta < 1.50) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.851455) + 
    (1.40 <= eta && eta < 1.50) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.830037) + 
    (1.40 <= eta && eta < 1.50) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.80) * (0.808192) + 
    (1.40 <= eta && eta < 1.50) * (17.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.786471) + 
    (1.40 <= eta && eta < 1.50) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.765443) + 
    (1.40 <= eta && eta < 1.50) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.744357) + 
    (1.40 <= eta && eta < 1.50) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.723059) + 
    (1.40 <= eta && eta < 1.50) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.701115) + 
    (1.40 <= eta && eta < 1.50) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.679219) + 
    (1.40 <= eta && eta < 1.50) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.657259) + 
    (1.40 <= eta && eta < 1.50) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.633793) + 
    (1.40 <= eta && eta < 1.50) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611473) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.50 <= eta && eta < 1.60) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.999377) + 
    (1.50 <= eta && eta < 1.60) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.979737) + 
    (1.50 <= eta && eta < 1.60) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.959755) + 
    (1.50 <= eta && eta < 1.60) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.939987) + 
    (1.50 <= eta && eta < 1.60) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.919393) + 
    (1.50 <= eta && eta < 1.60) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.897551) + 
    (1.50 <= eta && eta < 1.60) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.875307) + 
    (1.50 <= eta && eta < 1.60) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.853345) + 
    (1.50 <= eta && eta < 1.60) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.832165) + 
    (1.50 <= eta && eta < 1.60) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.810506) + 
    (1.50 <= eta && eta < 1.60) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.788911) + 
    (1.50 <= eta && eta < 1.60) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.767943) + 
    (1.50 <= eta && eta < 1.60) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.746852) + 
    (1.50 <= eta && eta < 1.60) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.725478) + 
    (1.50 <= eta && eta < 1.60) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.704144) + 
    (1.50 <= eta && eta < 1.60) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.682446) + 
    (1.50 <= eta && eta < 1.60) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.660133) + 
    (1.50 <= eta && eta < 1.60) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.636792) + 
    (1.50 <= eta && eta < 1.60) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.612517) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.60 <= eta && eta < 1.70) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.10) * (0.999419) + 
    (1.60 <= eta && eta < 1.70) * (12.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.980743) + 
    (1.60 <= eta && eta < 1.70) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.960088) + 
    (1.60 <= eta && eta < 1.70) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.937415) + 
    (1.60 <= eta && eta < 1.70) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.915125) + 
    (1.60 <= eta && eta < 1.70) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.893414) + 
    (1.60 <= eta && eta < 1.70) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.871449) + 
    (1.60 <= eta && eta < 1.70) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.849856) + 
    (1.60 <= eta && eta < 1.70) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.829090) + 
    (1.60 <= eta && eta < 1.70) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.807889) + 
    (1.60 <= eta && eta < 1.70) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.785352) + 
    (1.60 <= eta && eta < 1.70) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.762438) + 
    (1.60 <= eta && eta < 1.70) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.740068) + 
    (1.60 <= eta && eta < 1.70) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.718126) + 
    (1.60 <= eta && eta < 1.70) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.696206) + 
    (1.60 <= eta && eta < 1.70) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.674404) + 
    (1.60 <= eta && eta < 1.70) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 28.80) * (0.652342) + 
    (1.60 <= eta && eta < 1.70) * (28.80 <= pt * cosh(eta) && pt * cosh(eta) < 36.40) * (0.627759) + 
    (1.60 <= eta && eta < 1.70) * (36.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610244) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.70 <= eta && eta < 1.80) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.999427) + 
    (1.70 <= eta && eta < 1.80) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.981008) + 
    (1.70 <= eta && eta < 1.80) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.960700) + 
    (1.70 <= eta && eta < 1.80) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.938370) + 
    (1.70 <= eta && eta < 1.80) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.916372) + 
    (1.70 <= eta && eta < 1.80) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.894899) + 
    (1.70 <= eta && eta < 1.80) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.873129) + 
    (1.70 <= eta && eta < 1.80) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.851684) + 
    (1.70 <= eta && eta < 1.80) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.80) * (0.829357) + 
    (1.70 <= eta && eta < 1.80) * (17.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.806754) + 
    (1.70 <= eta && eta < 1.80) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.784514) + 
    (1.70 <= eta && eta < 1.80) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.761892) + 
    (1.70 <= eta && eta < 1.80) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.739786) + 
    (1.70 <= eta && eta < 1.80) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.718072) + 
    (1.70 <= eta && eta < 1.80) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.696341) + 
    (1.70 <= eta && eta < 1.80) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.674673) + 
    (1.70 <= eta && eta < 1.80) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.00) * (0.652680) + 
    (1.70 <= eta && eta < 1.80) * (29.00 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.628190) + 
    (1.70 <= eta && eta < 1.80) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610525) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.80 <= eta && eta < 1.90) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.999406) + 
    (1.80 <= eta && eta < 1.90) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.980628) + 
    (1.80 <= eta && eta < 1.90) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.960286) + 
    (1.80 <= eta && eta < 1.90) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.938026) + 
    (1.80 <= eta && eta < 1.90) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.916141) + 
    (1.80 <= eta && eta < 1.90) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.894797) + 
    (1.80 <= eta && eta < 1.90) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.873164) + 
    (1.80 <= eta && eta < 1.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.851850) + 
    (1.80 <= eta && eta < 1.90) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.829651) + 
    (1.80 <= eta && eta < 1.90) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.807163) + 
    (1.80 <= eta && eta < 1.90) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.785020) + 
    (1.80 <= eta && eta < 1.90) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.30) * (0.762473) + 
    (1.80 <= eta && eta < 1.90) * (20.30 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.740416) + 
    (1.80 <= eta && eta < 1.90) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.718723) + 
    (1.80 <= eta && eta < 1.90) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.696981) + 
    (1.80 <= eta && eta < 1.90) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.675267) + 
    (1.80 <= eta && eta < 1.90) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.653182) + 
    (1.80 <= eta && eta < 1.90) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 36.40) * (0.628794) + 
    (1.80 <= eta && eta < 1.90) * (36.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610819) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (1.90 <= eta && eta < 2.00) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.999379) + 
    (1.90 <= eta && eta < 2.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.980126) + 
    (1.90 <= eta && eta < 2.00) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.959687) + 
    (1.90 <= eta && eta < 2.00) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.937449) + 
    (1.90 <= eta && eta < 2.00) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.915644) + 
    (1.90 <= eta && eta < 2.00) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.894409) + 
    (1.90 <= eta && eta < 2.00) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.872900) + 
    (1.90 <= eta && eta < 2.00) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.851713) + 
    (1.90 <= eta && eta < 2.00) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.829643) + 
    (1.90 <= eta && eta < 2.00) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.807280) + 
    (1.90 <= eta && eta < 2.00) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.785247) + 
    (1.90 <= eta && eta < 2.00) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.762796) + 
    (1.90 <= eta && eta < 2.00) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.740812) + 
    (1.90 <= eta && eta < 2.00) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.719168) + 
    (1.90 <= eta && eta < 2.00) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.697447) + 
    (1.90 <= eta && eta < 2.00) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.675721) + 
    (1.90 <= eta && eta < 2.00) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.653583) + 
    (1.90 <= eta && eta < 2.00) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.40) * (0.629212) + 
    (1.90 <= eta && eta < 2.00) * (36.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611042) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.00 <= eta && eta < 2.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.999430) + 
    (2.00 <= eta && eta < 2.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.980274) + 
    (2.00 <= eta && eta < 2.10) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.958791) + 
    (2.00 <= eta && eta < 2.10) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.936501) + 
    (2.00 <= eta && eta < 2.10) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.914729) + 
    (2.00 <= eta && eta < 2.10) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.893572) + 
    (2.00 <= eta && eta < 2.10) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.872170) + 
    (2.00 <= eta && eta < 2.10) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.851104) + 
    (2.00 <= eta && eta < 2.10) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.829169) + 
    (2.00 <= eta && eta < 2.10) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.806944) + 
    (2.00 <= eta && eta < 2.10) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.785043) + 
    (2.00 <= eta && eta < 2.10) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.762717) + 
    (2.00 <= eta && eta < 2.10) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.740844) + 
    (2.00 <= eta && eta < 2.10) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.719290) + 
    (2.00 <= eta && eta < 2.10) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.697639) + 
    (2.00 <= eta && eta < 2.10) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.675956) + 
    (2.00 <= eta && eta < 2.10) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.653825) + 
    (2.00 <= eta && eta < 2.10) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.40) * (0.629545) + 
    (2.00 <= eta && eta < 2.10) * (36.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611235) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.10 <= eta && eta < 2.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.999378) + 
    (2.10 <= eta && eta < 2.20) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.980201) + 
    (2.10 <= eta && eta < 2.20) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.959989) + 
    (2.10 <= eta && eta < 2.20) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.938001) + 
    (2.10 <= eta && eta < 2.20) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.916426) + 
    (2.10 <= eta && eta < 2.20) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.895389) + 
    (2.10 <= eta && eta < 2.20) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.874052) + 
    (2.10 <= eta && eta < 2.20) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.853005) + 
    (2.10 <= eta && eta < 2.20) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.831047) + 
    (2.10 <= eta && eta < 2.20) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.808762) + 
    (2.10 <= eta && eta < 2.20) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.786769) + 
    (2.10 <= eta && eta < 2.20) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.764319) + 
    (2.10 <= eta && eta < 2.20) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.742297) + 
    (2.10 <= eta && eta < 2.20) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.720574) + 
    (2.10 <= eta && eta < 2.20) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.698730) + 
    (2.10 <= eta && eta < 2.20) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.676833) + 
    (2.10 <= eta && eta < 2.20) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.654464) + 
    (2.10 <= eta && eta < 2.20) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.630191) + 
    (2.10 <= eta && eta < 2.20) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611490) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.20 <= eta && eta < 2.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.999408) + 
    (2.20 <= eta && eta < 2.30) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.979864) + 
    (2.20 <= eta && eta < 2.30) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.958336) + 
    (2.20 <= eta && eta < 2.30) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.936107) + 
    (2.20 <= eta && eta < 2.30) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.914440) + 
    (2.20 <= eta && eta < 2.30) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.893403) + 
    (2.20 <= eta && eta < 2.30) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.872130) + 
    (2.20 <= eta && eta < 2.30) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.851191) + 
    (2.20 <= eta && eta < 2.30) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.829380) + 
    (2.20 <= eta && eta < 2.30) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.807270) + 
    (2.20 <= eta && eta < 2.30) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.70) * (0.785465) + 
    (2.20 <= eta && eta < 2.30) * (19.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.763219) + 
    (2.20 <= eta && eta < 2.30) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.741400) + 
    (2.20 <= eta && eta < 2.30) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.719876) + 
    (2.20 <= eta && eta < 2.30) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.698225) + 
    (2.20 <= eta && eta < 2.30) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.676506) + 
    (2.20 <= eta && eta < 2.30) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.40) * (0.654297) + 
    (2.20 <= eta && eta < 2.30) * (29.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.30) * (0.630155) + 
    (2.20 <= eta && eta < 2.30) * (36.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611542) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.30 <= eta && eta < 2.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.999433) + 
    (2.30 <= eta && eta < 2.40) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.980427) + 
    (2.30 <= eta && eta < 2.40) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.959228) + 
    (2.30 <= eta && eta < 2.40) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.937219) + 
    (2.30 <= eta && eta < 2.40) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.915695) + 
    (2.30 <= eta && eta < 2.40) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.894746) + 
    (2.30 <= eta && eta < 2.40) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.873520) + 
    (2.30 <= eta && eta < 2.40) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.852593) + 
    (2.30 <= eta && eta < 2.40) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.830766) + 
    (2.30 <= eta && eta < 2.40) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.808612) + 
    (2.30 <= eta && eta < 2.40) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.70) * (0.786740) + 
    (2.30 <= eta && eta < 2.40) * (19.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.764402) + 
    (2.30 <= eta && eta < 2.40) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.742475) + 
    (2.30 <= eta && eta < 2.40) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.720826) + 
    (2.30 <= eta && eta < 2.40) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.699033) + 
    (2.30 <= eta && eta < 2.40) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.677157) + 
    (2.30 <= eta && eta < 2.40) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.655106) + 
    (2.30 <= eta && eta < 2.40) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 35.90) * (0.631229) + 
    (2.30 <= eta && eta < 2.40) * (35.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611899) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.40 <= eta && eta < 2.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.999428) + 
    (2.40 <= eta && eta < 2.50) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.980322) + 
    (2.40 <= eta && eta < 2.50) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.959060) + 
    (2.40 <= eta && eta < 2.50) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.937010) + 
    (2.40 <= eta && eta < 2.50) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.915458) + 
    (2.40 <= eta && eta < 2.50) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.894492) + 
    (2.40 <= eta && eta < 2.50) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.873257) + 
    (2.40 <= eta && eta < 2.50) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.852328) + 
    (2.40 <= eta && eta < 2.50) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.830504) + 
    (2.40 <= eta && eta < 2.50) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.808357) + 
    (2.40 <= eta && eta < 2.50) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.70) * (0.786498) + 
    (2.40 <= eta && eta < 2.50) * (19.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.764178) + 
    (2.40 <= eta && eta < 2.50) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.742271) + 
    (2.40 <= eta && eta < 2.50) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.720646) + 
    (2.40 <= eta && eta < 2.50) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.698879) + 
    (2.40 <= eta && eta < 2.50) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.677034) + 
    (2.40 <= eta && eta < 2.50) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.655015) + 
    (2.40 <= eta && eta < 2.50) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 35.90) * (0.631176) + 
    (2.40 <= eta && eta < 2.50) * (35.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611879) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.50 <= eta && eta < 2.60) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.999439) + 
    (2.50 <= eta && eta < 2.60) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.980370) + 
    (2.50 <= eta && eta < 2.60) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.958748) + 
    (2.50 <= eta && eta < 2.60) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.936272) + 
    (2.50 <= eta && eta < 2.60) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.914314) + 
    (2.50 <= eta && eta < 2.60) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.892985) + 
    (2.50 <= eta && eta < 2.60) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.871426) + 
    (2.50 <= eta && eta < 2.60) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.850227) + 
    (2.50 <= eta && eta < 2.60) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.828176) + 
    (2.50 <= eta && eta < 2.60) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.70) * (0.805862) + 
    (2.50 <= eta && eta < 2.60) * (18.70 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.783903) + 
    (2.50 <= eta && eta < 2.60) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.761550) + 
    (2.50 <= eta && eta < 2.60) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.739683) + 
    (2.50 <= eta && eta < 2.60) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.718173) + 
    (2.50 <= eta && eta < 2.60) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.696603) + 
    (2.50 <= eta && eta < 2.60) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.675044) + 
    (2.50 <= eta && eta < 2.60) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.653090) + 
    (2.50 <= eta && eta < 2.60) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.628801) + 
    (2.50 <= eta && eta < 2.60) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610884) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.60 <= eta && eta < 2.70) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.20) * (0.999380) + 
    (2.60 <= eta && eta < 2.70) * (12.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.979989) + 
    (2.60 <= eta && eta < 2.70) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.959094) + 
    (2.60 <= eta && eta < 2.70) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.936345) + 
    (2.60 <= eta && eta < 2.70) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.914076) + 
    (2.60 <= eta && eta < 2.70) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.892440) + 
    (2.60 <= eta && eta < 2.70) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.870584) + 
    (2.60 <= eta && eta < 2.70) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.849117) + 
    (2.60 <= eta && eta < 2.70) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.80) * (0.826825) + 
    (2.60 <= eta && eta < 2.70) * (17.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.804309) + 
    (2.60 <= eta && eta < 2.70) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.782201) + 
    (2.60 <= eta && eta < 2.70) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.759752) + 
    (2.60 <= eta && eta < 2.70) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.737851) + 
    (2.60 <= eta && eta < 2.70) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.716370) + 
    (2.60 <= eta && eta < 2.70) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.694900) + 
    (2.60 <= eta && eta < 2.70) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.673023) + 
    (2.60 <= eta && eta < 2.70) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.650533) + 
    (2.60 <= eta && eta < 2.70) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.625118) + 
    (2.60 <= eta && eta < 2.70) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609492) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.70 <= eta && eta < 2.80) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.00) * (0.999427) + 
    (2.70 <= eta && eta < 2.80) * (12.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.980834) + 
    (2.70 <= eta && eta < 2.80) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.960037) + 
    (2.70 <= eta && eta < 2.80) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.938840) + 
    (2.70 <= eta && eta < 2.80) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.918248) + 
    (2.70 <= eta && eta < 2.80) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.896469) + 
    (2.70 <= eta && eta < 2.80) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.874329) + 
    (2.70 <= eta && eta < 2.80) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.852492) + 
    (2.70 <= eta && eta < 2.80) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.831447) + 
    (2.70 <= eta && eta < 2.80) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.809933) + 
    (2.70 <= eta && eta < 2.80) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.788483) + 
    (2.70 <= eta && eta < 2.80) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.767651) + 
    (2.70 <= eta && eta < 2.80) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.746688) + 
    (2.70 <= eta && eta < 2.80) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.725431) + 
    (2.70 <= eta && eta < 2.80) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.703428) + 
    (2.70 <= eta && eta < 2.80) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.80) * (0.681357) + 
    (2.70 <= eta && eta < 2.80) * (24.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.659083) + 
    (2.70 <= eta && eta < 2.80) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.635658) + 
    (2.70 <= eta && eta < 2.80) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.612337) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.80 <= eta && eta < 2.90) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.90) * (0.999399) + 
    (2.80 <= eta && eta < 2.90) * (11.90 <= pt * cosh(eta) && pt * cosh(eta) < 12.90) * (0.980128) + 
    (2.80 <= eta && eta < 2.90) * (12.90 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.960173) + 
    (2.80 <= eta && eta < 2.90) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.940341) + 
    (2.80 <= eta && eta < 2.90) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.919639) + 
    (2.80 <= eta && eta < 2.90) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.897662) + 
    (2.80 <= eta && eta < 2.90) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.875274) + 
    (2.80 <= eta && eta < 2.90) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.853171) + 
    (2.80 <= eta && eta < 2.90) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.20) * (0.831865) + 
    (2.80 <= eta && eta < 2.90) * (17.20 <= pt * cosh(eta) && pt * cosh(eta) < 17.90) * (0.810090) + 
    (2.80 <= eta && eta < 2.90) * (17.90 <= pt * cosh(eta) && pt * cosh(eta) < 18.60) * (0.788396) + 
    (2.80 <= eta && eta < 2.90) * (18.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.767351) + 
    (2.80 <= eta && eta < 2.90) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 20.30) * (0.746206) + 
    (2.80 <= eta && eta < 2.90) * (20.30 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.724804) + 
    (2.80 <= eta && eta < 2.90) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.703471) + 
    (2.80 <= eta && eta < 2.90) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.681809) + 
    (2.80 <= eta && eta < 2.90) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.659577) + 
    (2.80 <= eta && eta < 2.90) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.636380) + 
    (2.80 <= eta && eta < 2.90) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.612308) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (2.90 <= eta && eta < 3.00) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.80) * (0.999384) + 
    (2.90 <= eta && eta < 3.00) * (11.80 <= pt * cosh(eta) && pt * cosh(eta) < 12.80) * (0.979697) + 
    (2.90 <= eta && eta < 3.00) * (12.80 <= pt * cosh(eta) && pt * cosh(eta) < 13.50) * (0.959313) + 
    (2.90 <= eta && eta < 3.00) * (13.50 <= pt * cosh(eta) && pt * cosh(eta) < 14.10) * (0.939113) + 
    (2.90 <= eta && eta < 3.00) * (14.10 <= pt * cosh(eta) && pt * cosh(eta) < 14.70) * (0.918086) + 
    (2.90 <= eta && eta < 3.00) * (14.70 <= pt * cosh(eta) && pt * cosh(eta) < 15.30) * (0.895827) + 
    (2.90 <= eta && eta < 3.00) * (15.30 <= pt * cosh(eta) && pt * cosh(eta) < 15.90) * (0.873214) + 
    (2.90 <= eta && eta < 3.00) * (15.90 <= pt * cosh(eta) && pt * cosh(eta) < 16.50) * (0.850948) + 
    (2.90 <= eta && eta < 3.00) * (16.50 <= pt * cosh(eta) && pt * cosh(eta) < 17.10) * (0.829536) + 
    (2.90 <= eta && eta < 3.00) * (17.10 <= pt * cosh(eta) && pt * cosh(eta) < 17.80) * (0.807706) + 
    (2.90 <= eta && eta < 3.00) * (17.80 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.786010) + 
    (2.90 <= eta && eta < 3.00) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.765013) + 
    (2.90 <= eta && eta < 3.00) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.743964) + 
    (2.90 <= eta && eta < 3.00) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.722710) + 
    (2.90 <= eta && eta < 3.00) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.700817) + 
    (2.90 <= eta && eta < 3.00) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.678978) + 
    (2.90 <= eta && eta < 3.00) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.657079) + 
    (2.90 <= eta && eta < 3.00) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.80) * (0.633504) + 
    (2.90 <= eta && eta < 3.00) * (32.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.611367) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (3.00 <= eta && eta < 3.10) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.70) * (0.999380) + 
    (3.00 <= eta && eta < 3.10) * (11.70 <= pt * cosh(eta) && pt * cosh(eta) < 12.70) * (0.979522) + 
    (3.00 <= eta && eta < 3.10) * (12.70 <= pt * cosh(eta) && pt * cosh(eta) < 13.40) * (0.958847) + 
    (3.00 <= eta && eta < 3.10) * (13.40 <= pt * cosh(eta) && pt * cosh(eta) < 14.00) * (0.938370) + 
    (3.00 <= eta && eta < 3.10) * (14.00 <= pt * cosh(eta) && pt * cosh(eta) < 14.60) * (0.917083) + 
    (3.00 <= eta && eta < 3.10) * (14.60 <= pt * cosh(eta) && pt * cosh(eta) < 15.20) * (0.894585) + 
    (3.00 <= eta && eta < 3.10) * (15.20 <= pt * cosh(eta) && pt * cosh(eta) < 15.80) * (0.871771) + 
    (3.00 <= eta && eta < 3.10) * (15.80 <= pt * cosh(eta) && pt * cosh(eta) < 16.40) * (0.849348) + 
    (3.00 <= eta && eta < 3.10) * (16.40 <= pt * cosh(eta) && pt * cosh(eta) < 17.00) * (0.827826) + 
    (3.00 <= eta && eta < 3.10) * (17.00 <= pt * cosh(eta) && pt * cosh(eta) < 17.70) * (0.805925) + 
    (3.00 <= eta && eta < 3.10) * (17.70 <= pt * cosh(eta) && pt * cosh(eta) < 18.40) * (0.784199) + 
    (3.00 <= eta && eta < 3.10) * (18.40 <= pt * cosh(eta) && pt * cosh(eta) < 19.20) * (0.763214) + 
    (3.00 <= eta && eta < 3.10) * (19.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.742218) + 
    (3.00 <= eta && eta < 3.10) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.721061) + 
    (3.00 <= eta && eta < 3.10) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.699316) + 
    (3.00 <= eta && eta < 3.10) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.677678) + 
    (3.00 <= eta && eta < 3.10) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.655671) + 
    (3.00 <= eta && eta < 3.10) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.631723) + 
    (3.00 <= eta && eta < 3.10) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.610745) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (3.10 <= eta && eta < 3.20) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.60) * (0.999388) + 
    (3.10 <= eta && eta < 3.20) * (11.60 <= pt * cosh(eta) && pt * cosh(eta) < 12.60) * (0.979585) + 
    (3.10 <= eta && eta < 3.20) * (12.60 <= pt * cosh(eta) && pt * cosh(eta) < 13.30) * (0.958745) + 
    (3.10 <= eta && eta < 3.20) * (13.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.90) * (0.938071) + 
    (3.10 <= eta && eta < 3.20) * (13.90 <= pt * cosh(eta) && pt * cosh(eta) < 14.50) * (0.916580) + 
    (3.10 <= eta && eta < 3.20) * (14.50 <= pt * cosh(eta) && pt * cosh(eta) < 15.10) * (0.893880) + 
    (3.10 <= eta && eta < 3.20) * (15.10 <= pt * cosh(eta) && pt * cosh(eta) < 15.70) * (0.870884) + 
    (3.10 <= eta && eta < 3.20) * (15.70 <= pt * cosh(eta) && pt * cosh(eta) < 16.30) * (0.848309) + 
    (3.10 <= eta && eta < 3.20) * (16.30 <= pt * cosh(eta) && pt * cosh(eta) < 16.90) * (0.826669) + 
    (3.10 <= eta && eta < 3.20) * (16.90 <= pt * cosh(eta) && pt * cosh(eta) < 17.60) * (0.804680) + 
    (3.10 <= eta && eta < 3.20) * (17.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.30) * (0.782898) + 
    (3.10 <= eta && eta < 3.20) * (18.30 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.761891) + 
    (3.10 <= eta && eta < 3.20) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.740908) + 
    (3.10 <= eta && eta < 3.20) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.719802) + 
    (3.10 <= eta && eta < 3.20) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.698150) + 
    (3.10 <= eta && eta < 3.20) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.676099) + 
    (3.10 <= eta && eta < 3.20) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.653709) + 
    (3.10 <= eta && eta < 3.20) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.629318) + 
    (3.10 <= eta && eta < 3.20) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609960) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (3.20 <= eta && eta < 3.30) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.50) * (0.999404) + 
    (3.20 <= eta && eta < 3.30) * (11.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.50) * (0.979860) + 
    (3.20 <= eta && eta < 3.30) * (12.50 <= pt * cosh(eta) && pt * cosh(eta) < 13.20) * (0.958971) + 
    (3.20 <= eta && eta < 3.30) * (13.20 <= pt * cosh(eta) && pt * cosh(eta) < 13.80) * (0.938175) + 
    (3.20 <= eta && eta < 3.30) * (13.80 <= pt * cosh(eta) && pt * cosh(eta) < 14.40) * (0.916531) + 
    (3.20 <= eta && eta < 3.30) * (14.40 <= pt * cosh(eta) && pt * cosh(eta) < 15.00) * (0.893662) + 
    (3.20 <= eta && eta < 3.30) * (15.00 <= pt * cosh(eta) && pt * cosh(eta) < 15.60) * (0.870500) + 
    (3.20 <= eta && eta < 3.30) * (15.60 <= pt * cosh(eta) && pt * cosh(eta) < 16.20) * (0.847776) + 
    (3.20 <= eta && eta < 3.30) * (16.20 <= pt * cosh(eta) && pt * cosh(eta) < 16.80) * (0.826010) + 
    (3.20 <= eta && eta < 3.30) * (16.80 <= pt * cosh(eta) && pt * cosh(eta) < 17.50) * (0.803914) + 
    (3.20 <= eta && eta < 3.30) * (17.50 <= pt * cosh(eta) && pt * cosh(eta) < 18.20) * (0.782050) + 
    (3.20 <= eta && eta < 3.30) * (18.20 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.760991) + 
    (3.20 <= eta && eta < 3.30) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.90) * (0.739984) + 
    (3.20 <= eta && eta < 3.30) * (19.90 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.718885) + 
    (3.20 <= eta && eta < 3.30) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.697276) + 
    (3.20 <= eta && eta < 3.30) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.675311) + 
    (3.20 <= eta && eta < 3.30) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.653061) + 
    (3.20 <= eta && eta < 3.30) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.10) * (0.628603) + 
    (3.20 <= eta && eta < 3.30) * (34.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609653) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (3.30 <= eta && eta < 3.40) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.999429) + 
    (3.30 <= eta && eta < 3.40) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.980322) + 
    (3.30 <= eta && eta < 3.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 13.10) * (0.959491) + 
    (3.30 <= eta && eta < 3.40) * (13.10 <= pt * cosh(eta) && pt * cosh(eta) < 13.70) * (0.938643) + 
    (3.30 <= eta && eta < 3.40) * (13.70 <= pt * cosh(eta) && pt * cosh(eta) < 14.30) * (0.916891) + 
    (3.30 <= eta && eta < 3.40) * (14.30 <= pt * cosh(eta) && pt * cosh(eta) < 14.90) * (0.893882) + 
    (3.30 <= eta && eta < 3.40) * (14.90 <= pt * cosh(eta) && pt * cosh(eta) < 15.50) * (0.870568) + 
    (3.30 <= eta && eta < 3.40) * (15.50 <= pt * cosh(eta) && pt * cosh(eta) < 16.10) * (0.847696) + 
    (3.30 <= eta && eta < 3.40) * (16.10 <= pt * cosh(eta) && pt * cosh(eta) < 16.70) * (0.825795) + 
    (3.30 <= eta && eta < 3.40) * (16.70 <= pt * cosh(eta) && pt * cosh(eta) < 17.40) * (0.803576) + 
    (3.30 <= eta && eta < 3.40) * (17.40 <= pt * cosh(eta) && pt * cosh(eta) < 18.10) * (0.781608) + 
    (3.30 <= eta && eta < 3.40) * (18.10 <= pt * cosh(eta) && pt * cosh(eta) < 18.90) * (0.760467) + 
    (3.30 <= eta && eta < 3.40) * (18.90 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.739401) + 
    (3.30 <= eta && eta < 3.40) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.718270) + 
    (3.30 <= eta && eta < 3.40) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.696659) + 
    (3.30 <= eta && eta < 3.40) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.674731) + 
    (3.30 <= eta && eta < 3.40) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.652220) + 
    (3.30 <= eta && eta < 3.40) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.627361) + 
    (3.30 <= eta && eta < 3.40) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.609223) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 0.80) * (0.601058) + 
    (3.40 <= eta && eta < 3.50) * (0.80 <= pt * cosh(eta) && pt * cosh(eta) < 11.40) * (0.999367) + 
    (3.40 <= eta && eta < 3.50) * (11.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.30) * (0.980033) + 
    (3.40 <= eta && eta < 3.50) * (12.30 <= pt * cosh(eta) && pt * cosh(eta) < 13.00) * (0.960271) + 
    (3.40 <= eta && eta < 3.50) * (13.00 <= pt * cosh(eta) && pt * cosh(eta) < 13.60) * (0.939434) + 
    (3.40 <= eta && eta < 3.50) * (13.60 <= pt * cosh(eta) && pt * cosh(eta) < 14.20) * (0.917619) + 
    (3.40 <= eta && eta < 3.50) * (14.20 <= pt * cosh(eta) && pt * cosh(eta) < 14.80) * (0.894498) + 
    (3.40 <= eta && eta < 3.50) * (14.80 <= pt * cosh(eta) && pt * cosh(eta) < 15.40) * (0.871045) + 
    (3.40 <= eta && eta < 3.50) * (15.40 <= pt * cosh(eta) && pt * cosh(eta) < 16.00) * (0.848025) + 
    (3.40 <= eta && eta < 3.50) * (16.00 <= pt * cosh(eta) && pt * cosh(eta) < 16.60) * (0.825982) + 
    (3.40 <= eta && eta < 3.50) * (16.60 <= pt * cosh(eta) && pt * cosh(eta) < 17.30) * (0.803624) + 
    (3.40 <= eta && eta < 3.50) * (17.30 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.781528) + 
    (3.40 <= eta && eta < 3.50) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 18.80) * (0.760279) + 
    (3.40 <= eta && eta < 3.50) * (18.80 <= pt * cosh(eta) && pt * cosh(eta) < 19.70) * (0.739123) + 
    (3.40 <= eta && eta < 3.50) * (19.70 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.717923) + 
    (3.40 <= eta && eta < 3.50) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.696270) + 
    (3.40 <= eta && eta < 3.50) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.674333) + 
    (3.40 <= eta && eta < 3.50) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.651857) + 
    (3.40 <= eta && eta < 3.50) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.626836) + 
    (3.40 <= eta && eta < 3.50) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.608985) } 


  # --- protons ---

  add EfficiencyFormula {2212} {321} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.00 <= eta && eta < 1.10) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.000792) + 
    (1.00 <= eta && eta < 1.10) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.020173) + 
    (1.00 <= eta && eta < 1.10) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.040724) + 
    (1.00 <= eta && eta < 1.10) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.061613) + 
    (1.00 <= eta && eta < 1.10) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.082325) + 
    (1.00 <= eta && eta < 1.10) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 26.70) * (0.103121) + 
    (1.00 <= eta && eta < 1.10) * (26.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.124122) + 
    (1.00 <= eta && eta < 1.10) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.144827) + 
    (1.00 <= eta && eta < 1.10) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.165697) + 
    (1.00 <= eta && eta < 1.10) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.187160) + 
    (1.00 <= eta && eta < 1.10) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.30) * (0.208776) + 
    (1.00 <= eta && eta < 1.10) * (33.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.230151) + 
    (1.00 <= eta && eta < 1.10) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.60) * (0.250948) + 
    (1.00 <= eta && eta < 1.10) * (36.60 <= pt * cosh(eta) && pt * cosh(eta) < 38.50) * (0.271436) + 
    (1.00 <= eta && eta < 1.10) * (38.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.292204) + 
    (1.00 <= eta && eta < 1.10) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 43.20) * (0.313037) + 
    (1.00 <= eta && eta < 1.10) * (43.20 <= pt * cosh(eta) && pt * cosh(eta) < 46.20) * (0.333932) + 
    (1.00 <= eta && eta < 1.10) * (46.20 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.354921) + 
    (1.00 <= eta && eta < 1.10) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.365503) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.10 <= eta && eta < 1.20) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.000749) + 
    (1.10 <= eta && eta < 1.20) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.019212) + 
    (1.10 <= eta && eta < 1.20) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.038690) + 
    (1.10 <= eta && eta < 1.20) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.80) * (0.058609) + 
    (1.10 <= eta && eta < 1.20) * (24.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.079317) + 
    (1.10 <= eta && eta < 1.20) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.100267) + 
    (1.10 <= eta && eta < 1.20) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.120667) + 
    (1.10 <= eta && eta < 1.20) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.80) * (0.141681) + 
    (1.10 <= eta && eta < 1.20) * (29.80 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.162883) + 
    (1.10 <= eta && eta < 1.20) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.183892) + 
    (1.10 <= eta && eta < 1.20) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.205123) + 
    (1.10 <= eta && eta < 1.20) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 35.60) * (0.226190) + 
    (1.10 <= eta && eta < 1.20) * (35.60 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.247348) + 
    (1.10 <= eta && eta < 1.20) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.268698) + 
    (1.10 <= eta && eta < 1.20) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 41.60) * (0.289693) + 
    (1.10 <= eta && eta < 1.20) * (41.60 <= pt * cosh(eta) && pt * cosh(eta) < 44.20) * (0.310690) + 
    (1.10 <= eta && eta < 1.20) * (44.20 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.331732) + 
    (1.10 <= eta && eta < 1.20) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.350284) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.20 <= eta && eta < 1.30) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.000735) + 
    (1.20 <= eta && eta < 1.30) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.019251) + 
    (1.20 <= eta && eta < 1.30) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.039416) + 
    (1.20 <= eta && eta < 1.30) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.060388) + 
    (1.20 <= eta && eta < 1.30) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.90) * (0.081469) + 
    (1.20 <= eta && eta < 1.30) * (26.90 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.102762) + 
    (1.20 <= eta && eta < 1.30) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.124282) + 
    (1.20 <= eta && eta < 1.30) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.145489) + 
    (1.20 <= eta && eta < 1.30) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.166775) + 
    (1.20 <= eta && eta < 1.30) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.187791) + 
    (1.20 <= eta && eta < 1.30) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.208241) + 
    (1.20 <= eta && eta < 1.30) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.229170) + 
    (1.20 <= eta && eta < 1.30) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.250125) + 
    (1.20 <= eta && eta < 1.30) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 40.60) * (0.270672) + 
    (1.20 <= eta && eta < 1.30) * (40.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.90) * (0.291380) + 
    (1.20 <= eta && eta < 1.30) * (42.90 <= pt * cosh(eta) && pt * cosh(eta) < 45.60) * (0.312442) + 
    (1.20 <= eta && eta < 1.30) * (45.60 <= pt * cosh(eta) && pt * cosh(eta) < 48.70) * (0.333395) + 
    (1.20 <= eta && eta < 1.30) * (48.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.347295) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.30 <= eta && eta < 1.40) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.000754) + 
    (1.30 <= eta && eta < 1.40) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.019346) + 
    (1.30 <= eta && eta < 1.40) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.039034) + 
    (1.30 <= eta && eta < 1.40) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.059442) + 
    (1.30 <= eta && eta < 1.40) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.079966) + 
    (1.30 <= eta && eta < 1.40) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.100735) + 
    (1.30 <= eta && eta < 1.40) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.121778) + 
    (1.30 <= eta && eta < 1.40) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.142573) + 
    (1.30 <= eta && eta < 1.40) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.163507) + 
    (1.30 <= eta && eta < 1.40) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.184950) + 
    (1.30 <= eta && eta < 1.40) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 35.80) * (0.206480) + 
    (1.30 <= eta && eta < 1.40) * (35.80 <= pt * cosh(eta) && pt * cosh(eta) < 37.50) * (0.227721) + 
    (1.30 <= eta && eta < 1.40) * (37.50 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.248910) + 
    (1.30 <= eta && eta < 1.40) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 41.50) * (0.270149) + 
    (1.30 <= eta && eta < 1.40) * (41.50 <= pt * cosh(eta) && pt * cosh(eta) < 43.90) * (0.291367) + 
    (1.30 <= eta && eta < 1.40) * (43.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.60) * (0.312385) + 
    (1.30 <= eta && eta < 1.40) * (46.60 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.333232) + 
    (1.30 <= eta && eta < 1.40) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.343952) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.40 <= eta && eta < 1.50) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.000753) + 
    (1.40 <= eta && eta < 1.50) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.019627) + 
    (1.40 <= eta && eta < 1.50) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.039570) + 
    (1.40 <= eta && eta < 1.50) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.059653) + 
    (1.40 <= eta && eta < 1.50) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.079793) + 
    (1.40 <= eta && eta < 1.50) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.100160) + 
    (1.40 <= eta && eta < 1.50) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.50) * (0.120805) + 
    (1.40 <= eta && eta < 1.50) * (30.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.141231) + 
    (1.40 <= eta && eta < 1.50) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.161827) + 
    (1.40 <= eta && eta < 1.50) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.182963) + 
    (1.40 <= eta && eta < 1.50) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.30) * (0.204230) + 
    (1.40 <= eta && eta < 1.50) * (36.30 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.225258) + 
    (1.40 <= eta && eta < 1.50) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.90) * (0.246285) + 
    (1.40 <= eta && eta < 1.50) * (39.90 <= pt * cosh(eta) && pt * cosh(eta) < 42.00) * (0.267412) + 
    (1.40 <= eta && eta < 1.50) * (42.00 <= pt * cosh(eta) && pt * cosh(eta) < 44.40) * (0.288573) + 
    (1.40 <= eta && eta < 1.50) * (44.40 <= pt * cosh(eta) && pt * cosh(eta) < 47.10) * (0.309587) + 
    (1.40 <= eta && eta < 1.50) * (47.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.329525) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.50 <= eta && eta < 1.60) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.70) * (0.000735) + 
    (1.50 <= eta && eta < 1.60) * (21.70 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.019189) + 
    (1.50 <= eta && eta < 1.60) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.039219) + 
    (1.50 <= eta && eta < 1.60) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.059595) + 
    (1.50 <= eta && eta < 1.60) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.40) * (0.080172) + 
    (1.50 <= eta && eta < 1.60) * (28.40 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.101005) + 
    (1.50 <= eta && eta < 1.60) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.121330) + 
    (1.50 <= eta && eta < 1.60) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 32.40) * (0.142193) + 
    (1.50 <= eta && eta < 1.60) * (32.40 <= pt * cosh(eta) && pt * cosh(eta) < 33.80) * (0.163178) + 
    (1.50 <= eta && eta < 1.60) * (33.80 <= pt * cosh(eta) && pt * cosh(eta) < 35.30) * (0.183926) + 
    (1.50 <= eta && eta < 1.60) * (35.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.204812) + 
    (1.50 <= eta && eta < 1.60) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.225483) + 
    (1.50 <= eta && eta < 1.60) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 40.50) * (0.246179) + 
    (1.50 <= eta && eta < 1.60) * (40.50 <= pt * cosh(eta) && pt * cosh(eta) < 42.60) * (0.267007) + 
    (1.50 <= eta && eta < 1.60) * (42.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.287906) + 
    (1.50 <= eta && eta < 1.60) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 47.70) * (0.308703) + 
    (1.50 <= eta && eta < 1.60) * (47.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.326541) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.60 <= eta && eta < 1.70) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.000751) + 
    (1.60 <= eta && eta < 1.70) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.019346) + 
    (1.60 <= eta && eta < 1.70) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.039161) + 
    (1.60 <= eta && eta < 1.70) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.059970) + 
    (1.60 <= eta && eta < 1.70) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.80) * (0.081053) + 
    (1.60 <= eta && eta < 1.70) * (28.80 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.101651) + 
    (1.60 <= eta && eta < 1.70) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.122495) + 
    (1.60 <= eta && eta < 1.70) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 32.90) * (0.143833) + 
    (1.60 <= eta && eta < 1.70) * (32.90 <= pt * cosh(eta) && pt * cosh(eta) < 34.30) * (0.164504) + 
    (1.60 <= eta && eta < 1.70) * (34.30 <= pt * cosh(eta) && pt * cosh(eta) < 35.80) * (0.184936) + 
    (1.60 <= eta && eta < 1.70) * (35.80 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.205509) + 
    (1.60 <= eta && eta < 1.70) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 39.20) * (0.226469) + 
    (1.60 <= eta && eta < 1.70) * (39.20 <= pt * cosh(eta) && pt * cosh(eta) < 41.10) * (0.247385) + 
    (1.60 <= eta && eta < 1.70) * (41.10 <= pt * cosh(eta) && pt * cosh(eta) < 43.30) * (0.268326) + 
    (1.60 <= eta && eta < 1.70) * (43.30 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.289271) + 
    (1.60 <= eta && eta < 1.70) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 48.50) * (0.310020) + 
    (1.60 <= eta && eta < 1.70) * (48.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.325235) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.70 <= eta && eta < 1.80) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.000751) + 
    (1.70 <= eta && eta < 1.80) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.019708) + 
    (1.70 <= eta && eta < 1.80) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.040014) + 
    (1.70 <= eta && eta < 1.80) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 27.70) * (0.060761) + 
    (1.70 <= eta && eta < 1.80) * (27.70 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.081705) + 
    (1.70 <= eta && eta < 1.80) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.50) * (0.102896) + 
    (1.70 <= eta && eta < 1.80) * (30.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.124309) + 
    (1.70 <= eta && eta < 1.80) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 33.30) * (0.145412) + 
    (1.70 <= eta && eta < 1.80) * (33.30 <= pt * cosh(eta) && pt * cosh(eta) < 34.80) * (0.166545) + 
    (1.70 <= eta && eta < 1.80) * (34.80 <= pt * cosh(eta) && pt * cosh(eta) < 36.30) * (0.187371) + 
    (1.70 <= eta && eta < 1.80) * (36.30 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.208231) + 
    (1.70 <= eta && eta < 1.80) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.80) * (0.229385) + 
    (1.70 <= eta && eta < 1.80) * (39.80 <= pt * cosh(eta) && pt * cosh(eta) < 41.80) * (0.250379) + 
    (1.70 <= eta && eta < 1.80) * (41.80 <= pt * cosh(eta) && pt * cosh(eta) < 44.00) * (0.271317) + 
    (1.70 <= eta && eta < 1.80) * (44.00 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.292147) + 
    (1.70 <= eta && eta < 1.80) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 49.40) * (0.313073) + 
    (1.70 <= eta && eta < 1.80) * (49.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.325254) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.80 <= eta && eta < 1.90) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.000731) + 
    (1.80 <= eta && eta < 1.90) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.019328) + 
    (1.80 <= eta && eta < 1.90) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.039326) + 
    (1.80 <= eta && eta < 1.90) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.059822) + 
    (1.80 <= eta && eta < 1.90) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.080561) + 
    (1.80 <= eta && eta < 1.90) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.60) * (0.101587) + 
    (1.80 <= eta && eta < 1.90) * (30.60 <= pt * cosh(eta) && pt * cosh(eta) < 32.00) * (0.122870) + 
    (1.80 <= eta && eta < 1.90) * (32.00 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.143877) + 
    (1.80 <= eta && eta < 1.90) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.164942) + 
    (1.80 <= eta && eta < 1.90) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.186393) + 
    (1.80 <= eta && eta < 1.90) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.207825) + 
    (1.80 <= eta && eta < 1.90) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 40.00) * (0.228884) + 
    (1.80 <= eta && eta < 1.90) * (40.00 <= pt * cosh(eta) && pt * cosh(eta) < 42.00) * (0.249798) + 
    (1.80 <= eta && eta < 1.90) * (42.00 <= pt * cosh(eta) && pt * cosh(eta) < 44.20) * (0.270672) + 
    (1.80 <= eta && eta < 1.90) * (44.20 <= pt * cosh(eta) && pt * cosh(eta) < 46.70) * (0.291456) + 
    (1.80 <= eta && eta < 1.90) * (46.70 <= pt * cosh(eta) && pt * cosh(eta) < 49.60) * (0.312351) + 
    (1.80 <= eta && eta < 1.90) * (49.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.323866) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (1.90 <= eta && eta < 2.00) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.000768) + 
    (1.90 <= eta && eta < 2.00) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.80) * (0.019875) + 
    (1.90 <= eta && eta < 2.00) * (24.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.50) * (0.039975) + 
    (1.90 <= eta && eta < 2.00) * (26.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.060453) + 
    (1.90 <= eta && eta < 2.00) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.40) * (0.081116) + 
    (1.90 <= eta && eta < 2.00) * (29.40 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.102036) + 
    (1.90 <= eta && eta < 2.00) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.123197) + 
    (1.90 <= eta && eta < 2.00) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.144078) + 
    (1.90 <= eta && eta < 2.00) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.165017) + 
    (1.90 <= eta && eta < 2.00) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.186345) + 
    (1.90 <= eta && eta < 2.00) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.40) * (0.207662) + 
    (1.90 <= eta && eta < 2.00) * (38.40 <= pt * cosh(eta) && pt * cosh(eta) < 40.20) * (0.228619) + 
    (1.90 <= eta && eta < 2.00) * (40.20 <= pt * cosh(eta) && pt * cosh(eta) < 42.20) * (0.249444) + 
    (1.90 <= eta && eta < 2.00) * (42.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.40) * (0.270242) + 
    (1.90 <= eta && eta < 2.00) * (44.40 <= pt * cosh(eta) && pt * cosh(eta) < 46.90) * (0.290963) + 
    (1.90 <= eta && eta < 2.00) * (46.90 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.311813) + 
    (1.90 <= eta && eta < 2.00) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.322648) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.00 <= eta && eta < 2.10) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.000763) + 
    (2.00 <= eta && eta < 2.10) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.019751) + 
    (2.00 <= eta && eta < 2.10) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.039688) + 
    (2.00 <= eta && eta < 2.10) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.060014) + 
    (2.00 <= eta && eta < 2.10) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.080541) + 
    (2.00 <= eta && eta < 2.10) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.90) * (0.101343) + 
    (2.00 <= eta && eta < 2.10) * (30.90 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.122402) + 
    (2.00 <= eta && eta < 2.10) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.143200) + 
    (2.00 <= eta && eta < 2.10) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.164074) + 
    (2.00 <= eta && eta < 2.10) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.185352) + 
    (2.00 <= eta && eta < 2.10) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.50) * (0.206635) + 
    (2.00 <= eta && eta < 2.10) * (38.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.30) * (0.227574) + 
    (2.00 <= eta && eta < 2.10) * (40.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.248396) + 
    (2.00 <= eta && eta < 2.10) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.269205) + 
    (2.00 <= eta && eta < 2.10) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 47.00) * (0.289951) + 
    (2.00 <= eta && eta < 2.10) * (47.00 <= pt * cosh(eta) && pt * cosh(eta) < 49.90) * (0.310839) + 
    (2.00 <= eta && eta < 2.10) * (49.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.321032) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.10 <= eta && eta < 2.20) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.000727) + 
    (2.10 <= eta && eta < 2.20) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.019127) + 
    (2.10 <= eta && eta < 2.20) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.038725) + 
    (2.10 <= eta && eta < 2.20) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.058821) + 
    (2.10 <= eta && eta < 2.20) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.079190) + 
    (2.10 <= eta && eta < 2.20) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.90) * (0.099884) + 
    (2.10 <= eta && eta < 2.20) * (30.90 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.120876) + 
    (2.10 <= eta && eta < 2.20) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.141641) + 
    (2.10 <= eta && eta < 2.20) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.162508) + 
    (2.10 <= eta && eta < 2.20) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.80) * (0.183801) + 
    (2.10 <= eta && eta < 2.20) * (36.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.50) * (0.205121) + 
    (2.10 <= eta && eta < 2.20) * (38.50 <= pt * cosh(eta) && pt * cosh(eta) < 40.30) * (0.226112) + 
    (2.10 <= eta && eta < 2.20) * (40.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.246999) + 
    (2.10 <= eta && eta < 2.20) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 44.50) * (0.267886) + 
    (2.10 <= eta && eta < 2.20) * (44.50 <= pt * cosh(eta) && pt * cosh(eta) < 47.00) * (0.288721) + 
    (2.10 <= eta && eta < 2.20) * (47.00 <= pt * cosh(eta) && pt * cosh(eta) < 49.90) * (0.309707) + 
    (2.10 <= eta && eta < 2.20) * (49.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.319951) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.20 <= eta && eta < 2.30) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.000745) + 
    (2.20 <= eta && eta < 2.30) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.019403) + 
    (2.20 <= eta && eta < 2.30) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.70) * (0.039056) + 
    (2.20 <= eta && eta < 2.30) * (26.70 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.059148) + 
    (2.20 <= eta && eta < 2.30) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.60) * (0.079483) + 
    (2.20 <= eta && eta < 2.30) * (29.60 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.100128) + 
    (2.20 <= eta && eta < 2.30) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 32.40) * (0.121061) + 
    (2.20 <= eta && eta < 2.30) * (32.40 <= pt * cosh(eta) && pt * cosh(eta) < 33.80) * (0.141765) + 
    (2.20 <= eta && eta < 2.30) * (33.80 <= pt * cosh(eta) && pt * cosh(eta) < 35.30) * (0.162570) + 
    (2.20 <= eta && eta < 2.30) * (35.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.183803) + 
    (2.20 <= eta && eta < 2.30) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.205065) + 
    (2.20 <= eta && eta < 2.30) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 40.40) * (0.226004) + 
    (2.20 <= eta && eta < 2.30) * (40.40 <= pt * cosh(eta) && pt * cosh(eta) < 42.40) * (0.246846) + 
    (2.20 <= eta && eta < 2.30) * (42.40 <= pt * cosh(eta) && pt * cosh(eta) < 44.60) * (0.267693) + 
    (2.20 <= eta && eta < 2.30) * (44.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.10) * (0.288496) + 
    (2.20 <= eta && eta < 2.30) * (47.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.309457) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.30 <= eta && eta < 2.40) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.000768) + 
    (2.30 <= eta && eta < 2.40) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.019754) + 
    (2.30 <= eta && eta < 2.40) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.039501) + 
    (2.30 <= eta && eta < 2.40) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.059615) + 
    (2.30 <= eta && eta < 2.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.079935) + 
    (2.30 <= eta && eta < 2.40) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.100543) + 
    (2.30 <= eta && eta < 2.40) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.121426) + 
    (2.30 <= eta && eta < 2.40) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.90) * (0.142073) + 
    (2.30 <= eta && eta < 2.40) * (33.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.40) * (0.162818) + 
    (2.30 <= eta && eta < 2.40) * (35.40 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.183988) + 
    (2.30 <= eta && eta < 2.40) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 38.70) * (0.205189) + 
    (2.30 <= eta && eta < 2.40) * (38.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.50) * (0.226070) + 
    (2.30 <= eta && eta < 2.40) * (40.50 <= pt * cosh(eta) && pt * cosh(eta) < 42.50) * (0.246859) + 
    (2.30 <= eta && eta < 2.40) * (42.50 <= pt * cosh(eta) && pt * cosh(eta) < 44.70) * (0.267658) + 
    (2.30 <= eta && eta < 2.40) * (44.70 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.288418) + 
    (2.30 <= eta && eta < 2.40) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.308989) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.40 <= eta && eta < 2.50) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.000773) + 
    (2.40 <= eta && eta < 2.50) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.019842) + 
    (2.40 <= eta && eta < 2.50) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.039636) + 
    (2.40 <= eta && eta < 2.50) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.059781) + 
    (2.40 <= eta && eta < 2.50) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.080124) + 
    (2.40 <= eta && eta < 2.50) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.100746) + 
    (2.40 <= eta && eta < 2.50) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.121639) + 
    (2.40 <= eta && eta < 2.50) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.90) * (0.142291) + 
    (2.40 <= eta && eta < 2.50) * (33.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.40) * (0.163037) + 
    (2.40 <= eta && eta < 2.50) * (35.40 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.184205) + 
    (2.40 <= eta && eta < 2.50) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 38.70) * (0.205401) + 
    (2.40 <= eta && eta < 2.50) * (38.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.50) * (0.226275) + 
    (2.40 <= eta && eta < 2.50) * (40.50 <= pt * cosh(eta) && pt * cosh(eta) < 42.50) * (0.247054) + 
    (2.40 <= eta && eta < 2.50) * (42.50 <= pt * cosh(eta) && pt * cosh(eta) < 44.70) * (0.267842) + 
    (2.40 <= eta && eta < 2.50) * (44.70 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.288590) + 
    (2.40 <= eta && eta < 2.50) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.309148) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.50 <= eta && eta < 2.60) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.000746) + 
    (2.50 <= eta && eta < 2.60) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.019539) + 
    (2.50 <= eta && eta < 2.60) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.039556) + 
    (2.50 <= eta && eta < 2.60) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.060020) + 
    (2.50 <= eta && eta < 2.60) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.080706) + 
    (2.50 <= eta && eta < 2.60) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 30.70) * (0.101669) + 
    (2.50 <= eta && eta < 2.60) * (30.70 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.122884) + 
    (2.50 <= eta && eta < 2.60) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 33.50) * (0.143824) + 
    (2.50 <= eta && eta < 2.60) * (33.50 <= pt * cosh(eta) && pt * cosh(eta) < 35.00) * (0.164826) + 
    (2.50 <= eta && eta < 2.60) * (35.00 <= pt * cosh(eta) && pt * cosh(eta) < 36.60) * (0.186217) + 
    (2.50 <= eta && eta < 2.60) * (36.60 <= pt * cosh(eta) && pt * cosh(eta) < 38.30) * (0.207595) + 
    (2.50 <= eta && eta < 2.60) * (38.30 <= pt * cosh(eta) && pt * cosh(eta) < 40.10) * (0.228608) + 
    (2.50 <= eta && eta < 2.60) * (40.10 <= pt * cosh(eta) && pt * cosh(eta) < 42.10) * (0.249484) + 
    (2.50 <= eta && eta < 2.60) * (42.10 <= pt * cosh(eta) && pt * cosh(eta) < 44.30) * (0.270328) + 
    (2.50 <= eta && eta < 2.60) * (44.30 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.291089) + 
    (2.50 <= eta && eta < 2.60) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 49.70) * (0.311971) + 
    (2.50 <= eta && eta < 2.60) * (49.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.323152) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.60 <= eta && eta < 2.70) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.000749) + 
    (2.60 <= eta && eta < 2.70) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.019720) + 
    (2.60 <= eta && eta < 2.70) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.040133) + 
    (2.60 <= eta && eta < 2.70) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.60) * (0.060997) + 
    (2.60 <= eta && eta < 2.70) * (27.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.00) * (0.082052) + 
    (2.60 <= eta && eta < 2.70) * (29.00 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.103347) + 
    (2.60 <= eta && eta < 2.70) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.124852) + 
    (2.60 <= eta && eta < 2.70) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.146033) + 
    (2.60 <= eta && eta < 2.70) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.167231) + 
    (2.60 <= eta && eta < 2.70) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.20) * (0.188109) + 
    (2.60 <= eta && eta < 2.70) * (36.20 <= pt * cosh(eta) && pt * cosh(eta) < 37.90) * (0.209008) + 
    (2.60 <= eta && eta < 2.70) * (37.90 <= pt * cosh(eta) && pt * cosh(eta) < 39.70) * (0.230189) + 
    (2.60 <= eta && eta < 2.70) * (39.70 <= pt * cosh(eta) && pt * cosh(eta) < 41.70) * (0.251197) + 
    (2.60 <= eta && eta < 2.70) * (41.70 <= pt * cosh(eta) && pt * cosh(eta) < 43.90) * (0.272136) + 
    (2.60 <= eta && eta < 2.70) * (43.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.40) * (0.292956) + 
    (2.60 <= eta && eta < 2.70) * (46.40 <= pt * cosh(eta) && pt * cosh(eta) < 49.30) * (0.313858) + 
    (2.60 <= eta && eta < 2.70) * (49.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.326348) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.70 <= eta && eta < 2.80) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.000734) + 
    (2.70 <= eta && eta < 2.80) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.019136) + 
    (2.70 <= eta && eta < 2.80) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.039036) + 
    (2.70 <= eta && eta < 2.80) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.059994) + 
    (2.70 <= eta && eta < 2.80) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.60) * (0.081246) + 
    (2.70 <= eta && eta < 2.80) * (28.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.90) * (0.102010) + 
    (2.70 <= eta && eta < 2.80) * (29.90 <= pt * cosh(eta) && pt * cosh(eta) < 31.20) * (0.122240) + 
    (2.70 <= eta && eta < 2.80) * (31.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.60) * (0.142986) + 
    (2.70 <= eta && eta < 2.80) * (32.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.163845) + 
    (2.70 <= eta && eta < 2.80) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 35.50) * (0.184464) + 
    (2.70 <= eta && eta < 2.80) * (35.50 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.205221) + 
    (2.70 <= eta && eta < 2.80) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 38.80) * (0.225768) + 
    (2.70 <= eta && eta < 2.80) * (38.80 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.246347) + 
    (2.70 <= eta && eta < 2.80) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 42.80) * (0.267066) + 
    (2.70 <= eta && eta < 2.80) * (42.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.20) * (0.287866) + 
    (2.70 <= eta && eta < 2.80) * (45.20 <= pt * cosh(eta) && pt * cosh(eta) < 48.00) * (0.308948) + 
    (2.70 <= eta && eta < 2.80) * (48.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.326040) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.80 <= eta && eta < 2.90) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.000755) + 
    (2.80 <= eta && eta < 2.90) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.019582) + 
    (2.80 <= eta && eta < 2.90) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.039299) + 
    (2.80 <= eta && eta < 2.90) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.059144) + 
    (2.80 <= eta && eta < 2.90) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.079825) + 
    (2.80 <= eta && eta < 2.90) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.100784) + 
    (2.80 <= eta && eta < 2.90) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.121242) + 
    (2.80 <= eta && eta < 2.90) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.142240) + 
    (2.80 <= eta && eta < 2.90) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.163357) + 
    (2.80 <= eta && eta < 2.90) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.184227) + 
    (2.80 <= eta && eta < 2.90) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.205225) + 
    (2.80 <= eta && eta < 2.90) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.40) * (0.225993) + 
    (2.80 <= eta && eta < 2.90) * (38.40 <= pt * cosh(eta) && pt * cosh(eta) < 40.30) * (0.246773) + 
    (2.80 <= eta && eta < 2.90) * (40.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.40) * (0.267668) + 
    (2.80 <= eta && eta < 2.90) * (42.40 <= pt * cosh(eta) && pt * cosh(eta) < 44.80) * (0.288618) + 
    (2.80 <= eta && eta < 2.90) * (44.80 <= pt * cosh(eta) && pt * cosh(eta) < 47.50) * (0.309448) + 
    (2.80 <= eta && eta < 2.90) * (47.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.327951) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (2.90 <= eta && eta < 3.00) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.40) * (0.000763) + 
    (2.90 <= eta && eta < 3.00) * (21.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.60) * (0.019795) + 
    (2.90 <= eta && eta < 3.00) * (23.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.039828) + 
    (2.90 <= eta && eta < 3.00) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.059973) + 
    (2.90 <= eta && eta < 3.00) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.080154) + 
    (2.90 <= eta && eta < 3.00) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.100548) + 
    (2.90 <= eta && eta < 3.00) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.50) * (0.121212) + 
    (2.90 <= eta && eta < 3.00) * (30.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.141647) + 
    (2.90 <= eta && eta < 3.00) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.162245) + 
    (2.90 <= eta && eta < 3.00) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.70) * (0.183377) + 
    (2.90 <= eta && eta < 3.00) * (34.70 <= pt * cosh(eta) && pt * cosh(eta) < 36.30) * (0.204635) + 
    (2.90 <= eta && eta < 3.00) * (36.30 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.225649) + 
    (2.90 <= eta && eta < 3.00) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.90) * (0.246659) + 
    (2.90 <= eta && eta < 3.00) * (39.90 <= pt * cosh(eta) && pt * cosh(eta) < 42.00) * (0.267765) + 
    (2.90 <= eta && eta < 3.00) * (42.00 <= pt * cosh(eta) && pt * cosh(eta) < 44.40) * (0.288902) + 
    (2.90 <= eta && eta < 3.00) * (44.40 <= pt * cosh(eta) && pt * cosh(eta) < 47.10) * (0.309890) + 
    (2.90 <= eta && eta < 3.00) * (47.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.329801) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (3.00 <= eta && eta < 3.10) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.000757) + 
    (3.00 <= eta && eta < 3.10) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.40) * (0.019793) + 
    (3.00 <= eta && eta < 3.10) * (23.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.040032) + 
    (3.00 <= eta && eta < 3.10) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.060403) + 
    (3.00 <= eta && eta < 3.10) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.70) * (0.080804) + 
    (3.00 <= eta && eta < 3.10) * (27.70 <= pt * cosh(eta) && pt * cosh(eta) < 29.00) * (0.101403) + 
    (3.00 <= eta && eta < 3.10) * (29.00 <= pt * cosh(eta) && pt * cosh(eta) < 30.30) * (0.122252) + 
    (3.00 <= eta && eta < 3.10) * (30.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.142847) + 
    (3.00 <= eta && eta < 3.10) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 33.00) * (0.163582) + 
    (3.00 <= eta && eta < 3.10) * (33.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.184829) + 
    (3.00 <= eta && eta < 3.10) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 36.10) * (0.206176) + 
    (3.00 <= eta && eta < 3.10) * (36.10 <= pt * cosh(eta) && pt * cosh(eta) < 37.80) * (0.227253) + 
    (3.00 <= eta && eta < 3.10) * (37.80 <= pt * cosh(eta) && pt * cosh(eta) < 39.70) * (0.248298) + 
    (3.00 <= eta && eta < 3.10) * (39.70 <= pt * cosh(eta) && pt * cosh(eta) < 41.80) * (0.269416) + 
    (3.00 <= eta && eta < 3.10) * (41.80 <= pt * cosh(eta) && pt * cosh(eta) < 44.20) * (0.290537) + 
    (3.00 <= eta && eta < 3.10) * (44.20 <= pt * cosh(eta) && pt * cosh(eta) < 46.90) * (0.311485) + 
    (3.00 <= eta && eta < 3.10) * (46.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.331973) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (3.10 <= eta && eta < 3.20) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.000741) + 
    (3.10 <= eta && eta < 3.20) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.019121) + 
    (3.10 <= eta && eta < 3.20) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.038687) + 
    (3.10 <= eta && eta < 3.20) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.059011) + 
    (3.10 <= eta && eta < 3.20) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.079478) + 
    (3.10 <= eta && eta < 3.20) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.70) * (0.100208) + 
    (3.10 <= eta && eta < 3.20) * (28.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.121227) + 
    (3.10 <= eta && eta < 3.20) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.142009) + 
    (3.10 <= eta && eta < 3.20) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.162941) + 
    (3.10 <= eta && eta < 3.20) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.184390) + 
    (3.10 <= eta && eta < 3.20) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 35.80) * (0.205933) + 
    (3.10 <= eta && eta < 3.20) * (35.80 <= pt * cosh(eta) && pt * cosh(eta) < 37.50) * (0.227192) + 
    (3.10 <= eta && eta < 3.20) * (37.50 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.248406) + 
    (3.10 <= eta && eta < 3.20) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 41.50) * (0.269674) + 
    (3.10 <= eta && eta < 3.20) * (41.50 <= pt * cosh(eta) && pt * cosh(eta) < 43.80) * (0.290490) + 
    (3.10 <= eta && eta < 3.20) * (43.80 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.311207) + 
    (3.10 <= eta && eta < 3.20) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 49.70) * (0.332209) + 
    (3.10 <= eta && eta < 3.20) * (49.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.343307) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (3.20 <= eta && eta < 3.30) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.000768) + 
    (3.20 <= eta && eta < 3.30) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.019633) + 
    (3.20 <= eta && eta < 3.30) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.039580) + 
    (3.20 <= eta && eta < 3.30) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.060212) + 
    (3.20 <= eta && eta < 3.30) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.080924) + 
    (3.20 <= eta && eta < 3.30) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 28.60) * (0.101848) + 
    (3.20 <= eta && eta < 3.30) * (28.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.90) * (0.123017) + 
    (3.20 <= eta && eta < 3.30) * (29.90 <= pt * cosh(eta) && pt * cosh(eta) < 31.20) * (0.143910) + 
    (3.20 <= eta && eta < 3.30) * (31.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.60) * (0.164917) + 
    (3.20 <= eta && eta < 3.30) * (32.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.10) * (0.186410) + 
    (3.20 <= eta && eta < 3.30) * (34.10 <= pt * cosh(eta) && pt * cosh(eta) < 35.70) * (0.207967) + 
    (3.20 <= eta && eta < 3.30) * (35.70 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.229213) + 
    (3.20 <= eta && eta < 3.30) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 39.30) * (0.250388) + 
    (3.20 <= eta && eta < 3.30) * (39.30 <= pt * cosh(eta) && pt * cosh(eta) < 41.40) * (0.271594) + 
    (3.20 <= eta && eta < 3.30) * (41.40 <= pt * cosh(eta) && pt * cosh(eta) < 43.80) * (0.292761) + 
    (3.20 <= eta && eta < 3.30) * (43.80 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.313711) + 
    (3.20 <= eta && eta < 3.30) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 49.70) * (0.334474) + 
    (3.20 <= eta && eta < 3.30) * (49.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.345440) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (3.30 <= eta && eta < 3.40) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.000788) + 
    (3.30 <= eta && eta < 3.40) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.90) * (0.020005) + 
    (3.30 <= eta && eta < 3.40) * (22.90 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.040258) + 
    (3.30 <= eta && eta < 3.40) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.061145) + 
    (3.30 <= eta && eta < 3.40) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.082064) + 
    (3.30 <= eta && eta < 3.40) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.103157) + 
    (3.30 <= eta && eta < 3.40) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.80) * (0.124460) + 
    (3.30 <= eta && eta < 3.40) * (29.80 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.145453) + 
    (3.30 <= eta && eta < 3.40) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.166532) + 
    (3.30 <= eta && eta < 3.40) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.90) * (0.187357) + 
    (3.30 <= eta && eta < 3.40) * (33.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.40) * (0.207638) + 
    (3.30 <= eta && eta < 3.40) * (35.40 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.228413) + 
    (3.30 <= eta && eta < 3.40) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 38.90) * (0.249237) + 
    (3.30 <= eta && eta < 3.40) * (38.90 <= pt * cosh(eta) && pt * cosh(eta) < 40.90) * (0.269678) + 
    (3.30 <= eta && eta < 3.40) * (40.90 <= pt * cosh(eta) && pt * cosh(eta) < 43.20) * (0.290306) + 
    (3.30 <= eta && eta < 3.40) * (43.20 <= pt * cosh(eta) && pt * cosh(eta) < 45.90) * (0.311314) + 
    (3.30 <= eta && eta < 3.40) * (45.90 <= pt * cosh(eta) && pt * cosh(eta) < 49.00) * (0.332244) + 
    (3.30 <= eta && eta < 3.40) * (49.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.345265) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.500000) + 
    (3.40 <= eta && eta < 3.50) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.000743) + 
    (3.40 <= eta && eta < 3.50) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.019342) + 
    (3.40 <= eta && eta < 3.50) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.039448) + 
    (3.40 <= eta && eta < 3.50) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.060332) + 
    (3.40 <= eta && eta < 3.50) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.081319) + 
    (3.40 <= eta && eta < 3.50) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.102517) + 
    (3.40 <= eta && eta < 3.50) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.60) * (0.123948) + 
    (3.40 <= eta && eta < 3.50) * (29.60 <= pt * cosh(eta) && pt * cosh(eta) < 30.90) * (0.145075) + 
    (3.40 <= eta && eta < 3.50) * (30.90 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.166290) + 
    (3.40 <= eta && eta < 3.50) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.187247) + 
    (3.40 <= eta && eta < 3.50) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.207650) + 
    (3.40 <= eta && eta < 3.50) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.228542) + 
    (3.40 <= eta && eta < 3.50) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 38.70) * (0.249470) + 
    (3.40 <= eta && eta < 3.50) * (38.70 <= pt * cosh(eta) && pt * cosh(eta) < 40.70) * (0.270002) + 
    (3.40 <= eta && eta < 3.50) * (40.70 <= pt * cosh(eta) && pt * cosh(eta) < 43.00) * (0.290705) + 
    (3.40 <= eta && eta < 3.50) * (43.00 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.311775) + 
    (3.40 <= eta && eta < 3.50) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 48.80) * (0.332748) + 
    (3.40 <= eta && eta < 3.50) * (48.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.346376) } 

  add EfficiencyFormula {2212} {2212} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (1.00 <= eta && eta < 1.10) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.00) * (0.999348) + 
    (1.00 <= eta && eta < 1.10) * (18.00 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.980540) + 
    (1.00 <= eta && eta < 1.10) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.960532) + 
    (1.00 <= eta && eta < 1.10) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.60) * (0.939338) + 
    (1.00 <= eta && eta < 1.10) * (21.60 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.917587) + 
    (1.00 <= eta && eta < 1.10) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.40) * (0.895750) + 
    (1.00 <= eta && eta < 1.10) * (23.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.873577) + 
    (1.00 <= eta && eta < 1.10) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.851729) + 
    (1.00 <= eta && eta < 1.10) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.830689) + 
    (1.00 <= eta && eta < 1.10) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.809722) + 
    (1.00 <= eta && eta < 1.10) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.788275) + 
    (1.00 <= eta && eta < 1.10) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.40) * (0.766993) + 
    (1.00 <= eta && eta < 1.10) * (29.40 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.745712) + 
    (1.00 <= eta && eta < 1.10) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.40) * (0.724563) + 
    (1.00 <= eta && eta < 1.10) * (32.40 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.703482) + 
    (1.00 <= eta && eta < 1.10) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.682183) + 
    (1.00 <= eta && eta < 1.10) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 40.90) * (0.660314) + 
    (1.00 <= eta && eta < 1.10) * (40.90 <= pt * cosh(eta) && pt * cosh(eta) < 48.30) * (0.637091) + 
    (1.00 <= eta && eta < 1.10) * (48.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.625359) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (1.10 <= eta && eta < 1.20) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 18.50) * (0.999365) + 
    (1.10 <= eta && eta < 1.20) * (18.50 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.980488) + 
    (1.10 <= eta && eta < 1.20) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.20) * (0.960315) + 
    (1.10 <= eta && eta < 1.20) * (21.20 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.939740) + 
    (1.10 <= eta && eta < 1.20) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.918677) + 
    (1.10 <= eta && eta < 1.20) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.897519) + 
    (1.10 <= eta && eta < 1.20) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.875986) + 
    (1.10 <= eta && eta < 1.20) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.854702) + 
    (1.10 <= eta && eta < 1.20) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.833020) + 
    (1.10 <= eta && eta < 1.20) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.811444) + 
    (1.10 <= eta && eta < 1.20) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.90) * (0.790449) + 
    (1.10 <= eta && eta < 1.20) * (28.90 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.769533) + 
    (1.10 <= eta && eta < 1.20) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.748517) + 
    (1.10 <= eta && eta < 1.20) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.10) * (0.727517) + 
    (1.10 <= eta && eta < 1.20) * (33.10 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.706455) + 
    (1.10 <= eta && eta < 1.20) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 37.70) * (0.685021) + 
    (1.10 <= eta && eta < 1.20) * (37.70 <= pt * cosh(eta) && pt * cosh(eta) < 41.40) * (0.663362) + 
    (1.10 <= eta && eta < 1.20) * (41.40 <= pt * cosh(eta) && pt * cosh(eta) < 48.10) * (0.640661) + 
    (1.10 <= eta && eta < 1.20) * (48.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.628487) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (1.20 <= eta && eta < 1.30) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.00) * (0.999355) + 
    (1.20 <= eta && eta < 1.30) * (19.00 <= pt * cosh(eta) && pt * cosh(eta) < 20.60) * (0.980538) + 
    (1.20 <= eta && eta < 1.30) * (20.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.960090) + 
    (1.20 <= eta && eta < 1.30) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.939050) + 
    (1.20 <= eta && eta < 1.30) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.918492) + 
    (1.20 <= eta && eta < 1.30) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.897881) + 
    (1.20 <= eta && eta < 1.30) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.876905) + 
    (1.20 <= eta && eta < 1.30) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.856142) + 
    (1.20 <= eta && eta < 1.30) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.834944) + 
    (1.20 <= eta && eta < 1.30) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.40) * (0.813787) + 
    (1.20 <= eta && eta < 1.30) * (28.40 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.793129) + 
    (1.20 <= eta && eta < 1.30) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.70) * (0.772465) + 
    (1.20 <= eta && eta < 1.30) * (30.70 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.751611) + 
    (1.20 <= eta && eta < 1.30) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.730670) + 
    (1.20 <= eta && eta < 1.30) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 35.70) * (0.709548) + 
    (1.20 <= eta && eta < 1.30) * (35.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.30) * (0.687916) + 
    (1.20 <= eta && eta < 1.30) * (38.30 <= pt * cosh(eta) && pt * cosh(eta) < 41.90) * (0.666170) + 
    (1.20 <= eta && eta < 1.30) * (41.90 <= pt * cosh(eta) && pt * cosh(eta) < 48.10) * (0.643738) + 
    (1.20 <= eta && eta < 1.30) * (48.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.631380) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (1.30 <= eta && eta < 1.40) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.40) * (0.999373) + 
    (1.30 <= eta && eta < 1.40) * (19.40 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.981094) + 
    (1.30 <= eta && eta < 1.40) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.961438) + 
    (1.30 <= eta && eta < 1.40) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.20) * (0.941157) + 
    (1.30 <= eta && eta < 1.40) * (23.20 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.920158) + 
    (1.30 <= eta && eta < 1.40) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.898951) + 
    (1.30 <= eta && eta < 1.40) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.878449) + 
    (1.30 <= eta && eta < 1.40) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.857005) + 
    (1.30 <= eta && eta < 1.40) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.835164) + 
    (1.30 <= eta && eta < 1.40) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.813463) + 
    (1.30 <= eta && eta < 1.40) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.30) * (0.791437) + 
    (1.30 <= eta && eta < 1.40) * (30.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.60) * (0.769724) + 
    (1.30 <= eta && eta < 1.40) * (31.60 <= pt * cosh(eta) && pt * cosh(eta) < 33.10) * (0.748181) + 
    (1.30 <= eta && eta < 1.40) * (33.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.80) * (0.726923) + 
    (1.30 <= eta && eta < 1.40) * (34.80 <= pt * cosh(eta) && pt * cosh(eta) < 36.90) * (0.705865) + 
    (1.30 <= eta && eta < 1.40) * (36.90 <= pt * cosh(eta) && pt * cosh(eta) < 39.60) * (0.684675) + 
    (1.30 <= eta && eta < 1.40) * (39.60 <= pt * cosh(eta) && pt * cosh(eta) < 43.50) * (0.663164) + 
    (1.30 <= eta && eta < 1.40) * (43.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.641295) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (1.40 <= eta && eta < 1.50) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.999356) + 
    (1.40 <= eta && eta < 1.50) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.980314) + 
    (1.40 <= eta && eta < 1.50) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.959972) + 
    (1.40 <= eta && eta < 1.50) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.939834) + 
    (1.40 <= eta && eta < 1.50) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.919131) + 
    (1.40 <= eta && eta < 1.50) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.898295) + 
    (1.40 <= eta && eta < 1.50) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.50) * (0.878179) + 
    (1.40 <= eta && eta < 1.50) * (26.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.857142) + 
    (1.40 <= eta && eta < 1.50) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.835700) + 
    (1.40 <= eta && eta < 1.50) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.60) * (0.814364) + 
    (1.40 <= eta && eta < 1.50) * (29.60 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.792664) + 
    (1.40 <= eta && eta < 1.50) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.771217) + 
    (1.40 <= eta && eta < 1.50) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.749873) + 
    (1.40 <= eta && eta < 1.50) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.30) * (0.728739) + 
    (1.40 <= eta && eta < 1.50) * (35.30 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.707719) + 
    (1.40 <= eta && eta < 1.50) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 40.10) * (0.686472) + 
    (1.40 <= eta && eta < 1.50) * (40.10 <= pt * cosh(eta) && pt * cosh(eta) < 44.00) * (0.664783) + 
    (1.40 <= eta && eta < 1.50) * (44.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.643345) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (1.50 <= eta && eta < 1.60) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.10) * (0.999361) + 
    (1.50 <= eta && eta < 1.60) * (20.10 <= pt * cosh(eta) && pt * cosh(eta) < 21.80) * (0.980560) + 
    (1.50 <= eta && eta < 1.60) * (21.80 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.960696) + 
    (1.50 <= eta && eta < 1.60) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.940038) + 
    (1.50 <= eta && eta < 1.60) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.918647) + 
    (1.50 <= eta && eta < 1.60) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.897023) + 
    (1.50 <= eta && eta < 1.60) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.875024) + 
    (1.50 <= eta && eta < 1.60) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.853306) + 
    (1.50 <= eta && eta < 1.60) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.832354) + 
    (1.50 <= eta && eta < 1.60) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.20) * (0.811544) + 
    (1.50 <= eta && eta < 1.60) * (30.20 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.790394) + 
    (1.50 <= eta && eta < 1.60) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.769487) + 
    (1.50 <= eta && eta < 1.60) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.20) * (0.748660) + 
    (1.50 <= eta && eta < 1.60) * (34.20 <= pt * cosh(eta) && pt * cosh(eta) < 36.00) * (0.727424) + 
    (1.50 <= eta && eta < 1.60) * (36.00 <= pt * cosh(eta) && pt * cosh(eta) < 38.20) * (0.705955) + 
    (1.50 <= eta && eta < 1.60) * (38.20 <= pt * cosh(eta) && pt * cosh(eta) < 41.00) * (0.684629) + 
    (1.50 <= eta && eta < 1.60) * (41.00 <= pt * cosh(eta) && pt * cosh(eta) < 45.10) * (0.662965) + 
    (1.50 <= eta && eta < 1.60) * (45.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.643804) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (1.60 <= eta && eta < 1.70) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.40) * (0.999331) + 
    (1.60 <= eta && eta < 1.70) * (20.40 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.980092) + 
    (1.60 <= eta && eta < 1.70) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.30) * (0.960297) + 
    (1.60 <= eta && eta < 1.70) * (23.30 <= pt * cosh(eta) && pt * cosh(eta) < 24.40) * (0.939838) + 
    (1.60 <= eta && eta < 1.70) * (24.40 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.918702) + 
    (1.60 <= eta && eta < 1.70) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.40) * (0.897350) + 
    (1.60 <= eta && eta < 1.70) * (26.40 <= pt * cosh(eta) && pt * cosh(eta) < 27.40) * (0.875621) + 
    (1.60 <= eta && eta < 1.70) * (27.40 <= pt * cosh(eta) && pt * cosh(eta) < 28.40) * (0.854152) + 
    (1.60 <= eta && eta < 1.70) * (28.40 <= pt * cosh(eta) && pt * cosh(eta) < 29.40) * (0.833415) + 
    (1.60 <= eta && eta < 1.70) * (29.40 <= pt * cosh(eta) && pt * cosh(eta) < 30.50) * (0.812787) + 
    (1.60 <= eta && eta < 1.70) * (30.50 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.791786) + 
    (1.60 <= eta && eta < 1.70) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 33.00) * (0.770986) + 
    (1.60 <= eta && eta < 1.70) * (33.00 <= pt * cosh(eta) && pt * cosh(eta) < 34.50) * (0.750221) + 
    (1.60 <= eta && eta < 1.70) * (34.50 <= pt * cosh(eta) && pt * cosh(eta) < 36.30) * (0.728997) + 
    (1.60 <= eta && eta < 1.70) * (36.30 <= pt * cosh(eta) && pt * cosh(eta) < 38.50) * (0.707482) + 
    (1.60 <= eta && eta < 1.70) * (38.50 <= pt * cosh(eta) && pt * cosh(eta) < 41.30) * (0.686044) + 
    (1.60 <= eta && eta < 1.70) * (41.30 <= pt * cosh(eta) && pt * cosh(eta) < 45.30) * (0.664439) + 
    (1.60 <= eta && eta < 1.70) * (45.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.645456) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (1.70 <= eta && eta < 1.80) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.999379) + 
    (1.70 <= eta && eta < 1.80) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.981113) + 
    (1.70 <= eta && eta < 1.80) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.961171) + 
    (1.70 <= eta && eta < 1.80) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.60) * (0.940162) + 
    (1.70 <= eta && eta < 1.80) * (24.60 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.919254) + 
    (1.70 <= eta && eta < 1.80) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.898120) + 
    (1.70 <= eta && eta < 1.80) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.60) * (0.876590) + 
    (1.70 <= eta && eta < 1.80) * (27.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.60) * (0.855292) + 
    (1.70 <= eta && eta < 1.80) * (28.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.833694) + 
    (1.70 <= eta && eta < 1.80) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.812282) + 
    (1.70 <= eta && eta < 1.80) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.00) * (0.791504) + 
    (1.70 <= eta && eta < 1.80) * (32.00 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.770145) + 
    (1.70 <= eta && eta < 1.80) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.748958) + 
    (1.70 <= eta && eta < 1.80) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 36.70) * (0.728095) + 
    (1.70 <= eta && eta < 1.80) * (36.70 <= pt * cosh(eta) && pt * cosh(eta) < 38.90) * (0.706918) + 
    (1.70 <= eta && eta < 1.80) * (38.90 <= pt * cosh(eta) && pt * cosh(eta) < 41.80) * (0.685420) + 
    (1.70 <= eta && eta < 1.80) * (41.80 <= pt * cosh(eta) && pt * cosh(eta) < 45.90) * (0.663634) + 
    (1.70 <= eta && eta < 1.80) * (45.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.645858) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (1.80 <= eta && eta < 1.90) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.999339) + 
    (1.80 <= eta && eta < 1.90) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.40) * (0.980391) + 
    (1.80 <= eta && eta < 1.90) * (22.40 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.960255) + 
    (1.80 <= eta && eta < 1.90) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.80) * (0.939213) + 
    (1.80 <= eta && eta < 1.90) * (24.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.80) * (0.918359) + 
    (1.80 <= eta && eta < 1.90) * (25.80 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.897326) + 
    (1.80 <= eta && eta < 1.90) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.875930) + 
    (1.80 <= eta && eta < 1.90) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.80) * (0.854779) + 
    (1.80 <= eta && eta < 1.90) * (28.80 <= pt * cosh(eta) && pt * cosh(eta) < 29.90) * (0.833336) + 
    (1.80 <= eta && eta < 1.90) * (29.90 <= pt * cosh(eta) && pt * cosh(eta) < 31.00) * (0.812075) + 
    (1.80 <= eta && eta < 1.90) * (31.00 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.791439) + 
    (1.80 <= eta && eta < 1.90) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.770213) + 
    (1.80 <= eta && eta < 1.90) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.748478) + 
    (1.80 <= eta && eta < 1.90) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.727213) + 
    (1.80 <= eta && eta < 1.90) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.20) * (0.706313) + 
    (1.80 <= eta && eta < 1.90) * (39.20 <= pt * cosh(eta) && pt * cosh(eta) < 42.10) * (0.685068) + 
    (1.80 <= eta && eta < 1.90) * (42.10 <= pt * cosh(eta) && pt * cosh(eta) < 46.30) * (0.663255) + 
    (1.80 <= eta && eta < 1.90) * (46.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.646168) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (1.90 <= eta && eta < 2.00) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.999349) + 
    (1.90 <= eta && eta < 2.00) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.980635) + 
    (1.90 <= eta && eta < 2.00) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.960745) + 
    (1.90 <= eta && eta < 2.00) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.939924) + 
    (1.90 <= eta && eta < 2.00) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.919255) + 
    (1.90 <= eta && eta < 2.00) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.90) * (0.898376) + 
    (1.90 <= eta && eta < 2.00) * (26.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.877105) + 
    (1.90 <= eta && eta < 2.00) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 28.90) * (0.856046) + 
    (1.90 <= eta && eta < 2.00) * (28.90 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.834667) + 
    (1.90 <= eta && eta < 2.00) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.813442) + 
    (1.90 <= eta && eta < 2.00) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.792811) + 
    (1.90 <= eta && eta < 2.00) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.771564) + 
    (1.90 <= eta && eta < 2.00) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 35.30) * (0.749778) + 
    (1.90 <= eta && eta < 2.00) * (35.30 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.728433) + 
    (1.90 <= eta && eta < 2.00) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 39.30) * (0.707424) + 
    (1.90 <= eta && eta < 2.00) * (39.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.20) * (0.686038) + 
    (1.90 <= eta && eta < 2.00) * (42.20 <= pt * cosh(eta) && pt * cosh(eta) < 46.30) * (0.664289) + 
    (1.90 <= eta && eta < 2.00) * (46.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.647147) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (2.00 <= eta && eta < 2.10) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.999350) + 
    (2.00 <= eta && eta < 2.10) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 22.60) * (0.980688) + 
    (2.00 <= eta && eta < 2.10) * (22.60 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.960936) + 
    (2.00 <= eta && eta < 2.10) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.940260) + 
    (2.00 <= eta && eta < 2.10) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.919725) + 
    (2.00 <= eta && eta < 2.10) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.898967) + 
    (2.00 <= eta && eta < 2.10) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.877801) + 
    (2.00 <= eta && eta < 2.10) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.00) * (0.856829) + 
    (2.00 <= eta && eta < 2.10) * (29.00 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.835518) + 
    (2.00 <= eta && eta < 2.10) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.20) * (0.814341) + 
    (2.00 <= eta && eta < 2.10) * (31.20 <= pt * cosh(eta) && pt * cosh(eta) < 32.40) * (0.793738) + 
    (2.00 <= eta && eta < 2.10) * (32.40 <= pt * cosh(eta) && pt * cosh(eta) < 33.80) * (0.772498) + 
    (2.00 <= eta && eta < 2.10) * (33.80 <= pt * cosh(eta) && pt * cosh(eta) < 35.30) * (0.751363) + 
    (2.00 <= eta && eta < 2.10) * (35.30 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.730479) + 
    (2.00 <= eta && eta < 2.10) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 39.30) * (0.709203) + 
    (2.00 <= eta && eta < 2.10) * (39.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.10) * (0.687873) + 
    (2.00 <= eta && eta < 2.10) * (42.10 <= pt * cosh(eta) && pt * cosh(eta) < 46.10) * (0.666217) + 
    (2.00 <= eta && eta < 2.10) * (46.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.648389) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (2.10 <= eta && eta < 2.20) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.999331) + 
    (2.10 <= eta && eta < 2.20) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.980367) + 
    (2.10 <= eta && eta < 2.20) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.960541) + 
    (2.10 <= eta && eta < 2.20) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 25.10) * (0.939863) + 
    (2.10 <= eta && eta < 2.20) * (25.10 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.919362) + 
    (2.10 <= eta && eta < 2.20) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.898659) + 
    (2.10 <= eta && eta < 2.20) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.877561) + 
    (2.10 <= eta && eta < 2.20) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.10) * (0.856661) + 
    (2.10 <= eta && eta < 2.20) * (29.10 <= pt * cosh(eta) && pt * cosh(eta) < 30.20) * (0.835426) + 
    (2.10 <= eta && eta < 2.20) * (30.20 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.814321) + 
    (2.10 <= eta && eta < 2.20) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 32.50) * (0.793783) + 
    (2.10 <= eta && eta < 2.20) * (32.50 <= pt * cosh(eta) && pt * cosh(eta) < 33.90) * (0.772603) + 
    (2.10 <= eta && eta < 2.20) * (33.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.50) * (0.750854) + 
    (2.10 <= eta && eta < 2.20) * (35.50 <= pt * cosh(eta) && pt * cosh(eta) < 37.30) * (0.729512) + 
    (2.10 <= eta && eta < 2.20) * (37.30 <= pt * cosh(eta) && pt * cosh(eta) < 39.50) * (0.708467) + 
    (2.10 <= eta && eta < 2.20) * (39.50 <= pt * cosh(eta) && pt * cosh(eta) < 42.40) * (0.687002) + 
    (2.10 <= eta && eta < 2.20) * (42.40 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.665121) + 
    (2.10 <= eta && eta < 2.20) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.648182) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (2.20 <= eta && eta < 2.30) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.999363) + 
    (2.20 <= eta && eta < 2.30) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.980432) + 
    (2.20 <= eta && eta < 2.30) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.959984) + 
    (2.20 <= eta && eta < 2.30) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.939265) + 
    (2.20 <= eta && eta < 2.30) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.918772) + 
    (2.20 <= eta && eta < 2.30) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.898106) + 
    (2.20 <= eta && eta < 2.30) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.877065) + 
    (2.20 <= eta && eta < 2.30) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.856233) + 
    (2.20 <= eta && eta < 2.30) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.30) * (0.835074) + 
    (2.20 <= eta && eta < 2.30) * (30.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.814047) + 
    (2.20 <= eta && eta < 2.30) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 32.60) * (0.793585) + 
    (2.20 <= eta && eta < 2.30) * (32.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.772481) + 
    (2.20 <= eta && eta < 2.30) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 35.60) * (0.750803) + 
    (2.20 <= eta && eta < 2.30) * (35.60 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.729523) + 
    (2.20 <= eta && eta < 2.30) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 39.60) * (0.708529) + 
    (2.20 <= eta && eta < 2.30) * (39.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.50) * (0.687103) + 
    (2.20 <= eta && eta < 2.30) * (42.50 <= pt * cosh(eta) && pt * cosh(eta) < 46.60) * (0.665244) + 
    (2.20 <= eta && eta < 2.30) * (46.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.648481) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (2.30 <= eta && eta < 2.40) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.999336) + 
    (2.30 <= eta && eta < 2.40) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.90) * (0.979915) + 
    (2.30 <= eta && eta < 2.40) * (22.90 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.959285) + 
    (2.30 <= eta && eta < 2.40) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.30) * (0.938492) + 
    (2.30 <= eta && eta < 2.40) * (25.30 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.917985) + 
    (2.30 <= eta && eta < 2.40) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.897342) + 
    (2.30 <= eta && eta < 2.40) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.876349) + 
    (2.30 <= eta && eta < 2.40) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.30) * (0.855583) + 
    (2.30 <= eta && eta < 2.40) * (29.30 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.834500) + 
    (2.30 <= eta && eta < 2.40) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.813557) + 
    (2.30 <= eta && eta < 2.40) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.793180) + 
    (2.30 <= eta && eta < 2.40) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.10) * (0.772163) + 
    (2.30 <= eta && eta < 2.40) * (34.10 <= pt * cosh(eta) && pt * cosh(eta) < 35.70) * (0.750573) + 
    (2.30 <= eta && eta < 2.40) * (35.70 <= pt * cosh(eta) && pt * cosh(eta) < 37.50) * (0.729372) + 
    (2.30 <= eta && eta < 2.40) * (37.50 <= pt * cosh(eta) && pt * cosh(eta) < 39.70) * (0.708451) + 
    (2.30 <= eta && eta < 2.40) * (39.70 <= pt * cosh(eta) && pt * cosh(eta) < 42.60) * (0.687086) + 
    (2.30 <= eta && eta < 2.40) * (42.60 <= pt * cosh(eta) && pt * cosh(eta) < 46.70) * (0.665276) + 
    (2.30 <= eta && eta < 2.40) * (46.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.648713) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (2.40 <= eta && eta < 2.50) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.00) * (0.999385) + 
    (2.40 <= eta && eta < 2.50) * (21.00 <= pt * cosh(eta) && pt * cosh(eta) < 22.80) * (0.980880) + 
    (2.40 <= eta && eta < 2.50) * (22.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.960691) + 
    (2.40 <= eta && eta < 2.50) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.940147) + 
    (2.40 <= eta && eta < 2.50) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.919771) + 
    (2.40 <= eta && eta < 2.50) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.899182) + 
    (2.40 <= eta && eta < 2.50) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.878186) + 
    (2.40 <= eta && eta < 2.50) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.857370) + 
    (2.40 <= eta && eta < 2.50) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.30) * (0.836203) + 
    (2.40 <= eta && eta < 2.50) * (30.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.40) * (0.815148) + 
    (2.40 <= eta && eta < 2.50) * (31.40 <= pt * cosh(eta) && pt * cosh(eta) < 32.60) * (0.794641) + 
    (2.40 <= eta && eta < 2.50) * (32.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.00) * (0.773473) + 
    (2.40 <= eta && eta < 2.50) * (34.00 <= pt * cosh(eta) && pt * cosh(eta) < 35.60) * (0.751714) + 
    (2.40 <= eta && eta < 2.50) * (35.60 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.730339) + 
    (2.40 <= eta && eta < 2.50) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 39.60) * (0.709240) + 
    (2.40 <= eta && eta < 2.50) * (39.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.40) * (0.688051) + 
    (2.40 <= eta && eta < 2.50) * (42.40 <= pt * cosh(eta) && pt * cosh(eta) < 46.40) * (0.666489) + 
    (2.40 <= eta && eta < 2.50) * (46.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.649221) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (2.50 <= eta && eta < 2.60) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.999375) + 
    (2.50 <= eta && eta < 2.60) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.980538) + 
    (2.50 <= eta && eta < 2.60) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.959821) + 
    (2.50 <= eta && eta < 2.60) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.90) * (0.938769) + 
    (2.50 <= eta && eta < 2.60) * (24.90 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.917946) + 
    (2.50 <= eta && eta < 2.60) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.90) * (0.896967) + 
    (2.50 <= eta && eta < 2.60) * (26.90 <= pt * cosh(eta) && pt * cosh(eta) < 27.90) * (0.875638) + 
    (2.50 <= eta && eta < 2.60) * (27.90 <= pt * cosh(eta) && pt * cosh(eta) < 28.90) * (0.854560) + 
    (2.50 <= eta && eta < 2.60) * (28.90 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.833193) + 
    (2.50 <= eta && eta < 2.60) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.10) * (0.812008) + 
    (2.50 <= eta && eta < 2.60) * (31.10 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.791439) + 
    (2.50 <= eta && eta < 2.60) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 33.70) * (0.770277) + 
    (2.50 <= eta && eta < 2.60) * (33.70 <= pt * cosh(eta) && pt * cosh(eta) < 35.30) * (0.748599) + 
    (2.50 <= eta && eta < 2.60) * (35.30 <= pt * cosh(eta) && pt * cosh(eta) < 37.10) * (0.727378) + 
    (2.50 <= eta && eta < 2.60) * (37.10 <= pt * cosh(eta) && pt * cosh(eta) < 39.30) * (0.706508) + 
    (2.50 <= eta && eta < 2.60) * (39.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.20) * (0.685279) + 
    (2.50 <= eta && eta < 2.60) * (42.20 <= pt * cosh(eta) && pt * cosh(eta) < 46.40) * (0.663461) + 
    (2.50 <= eta && eta < 2.60) * (46.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.646521) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (2.60 <= eta && eta < 2.70) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.50) * (0.999329) + 
    (2.60 <= eta && eta < 2.70) * (20.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.20) * (0.980098) + 
    (2.60 <= eta && eta < 2.70) * (22.20 <= pt * cosh(eta) && pt * cosh(eta) < 23.40) * (0.960415) + 
    (2.60 <= eta && eta < 2.70) * (23.40 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.940083) + 
    (2.60 <= eta && eta < 2.70) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.50) * (0.919073) + 
    (2.60 <= eta && eta < 2.70) * (25.50 <= pt * cosh(eta) && pt * cosh(eta) < 26.50) * (0.897837) + 
    (2.60 <= eta && eta < 2.70) * (26.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.876213) + 
    (2.60 <= eta && eta < 2.70) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.854830) + 
    (2.60 <= eta && eta < 2.70) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.50) * (0.834161) + 
    (2.60 <= eta && eta < 2.70) * (29.50 <= pt * cosh(eta) && pt * cosh(eta) < 30.60) * (0.813583) + 
    (2.60 <= eta && eta < 2.70) * (30.60 <= pt * cosh(eta) && pt * cosh(eta) < 31.80) * (0.792615) + 
    (2.60 <= eta && eta < 2.70) * (31.80 <= pt * cosh(eta) && pt * cosh(eta) < 33.10) * (0.771829) + 
    (2.60 <= eta && eta < 2.70) * (33.10 <= pt * cosh(eta) && pt * cosh(eta) < 34.60) * (0.751057) + 
    (2.60 <= eta && eta < 2.70) * (34.60 <= pt * cosh(eta) && pt * cosh(eta) < 36.40) * (0.729805) + 
    (2.60 <= eta && eta < 2.70) * (36.40 <= pt * cosh(eta) && pt * cosh(eta) < 38.60) * (0.708237) + 
    (2.60 <= eta && eta < 2.70) * (38.60 <= pt * cosh(eta) && pt * cosh(eta) < 41.40) * (0.686720) + 
    (2.60 <= eta && eta < 2.70) * (41.40 <= pt * cosh(eta) && pt * cosh(eta) < 45.40) * (0.665007) + 
    (2.60 <= eta && eta < 2.70) * (45.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.646066) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (2.70 <= eta && eta < 2.80) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.20) * (0.999357) + 
    (2.70 <= eta && eta < 2.80) * (20.20 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.980532) + 
    (2.70 <= eta && eta < 2.80) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.960762) + 
    (2.70 <= eta && eta < 2.80) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.20) * (0.940220) + 
    (2.70 <= eta && eta < 2.80) * (24.20 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.918950) + 
    (2.70 <= eta && eta < 2.80) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.20) * (0.897439) + 
    (2.70 <= eta && eta < 2.80) * (26.20 <= pt * cosh(eta) && pt * cosh(eta) < 27.20) * (0.875543) + 
    (2.70 <= eta && eta < 2.80) * (27.20 <= pt * cosh(eta) && pt * cosh(eta) < 28.20) * (0.853913) + 
    (2.70 <= eta && eta < 2.80) * (28.20 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.833030) + 
    (2.70 <= eta && eta < 2.80) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.30) * (0.812272) + 
    (2.70 <= eta && eta < 2.80) * (30.30 <= pt * cosh(eta) && pt * cosh(eta) < 31.50) * (0.791158) + 
    (2.70 <= eta && eta < 2.80) * (31.50 <= pt * cosh(eta) && pt * cosh(eta) < 32.80) * (0.770268) + 
    (2.70 <= eta && eta < 2.80) * (32.80 <= pt * cosh(eta) && pt * cosh(eta) < 34.30) * (0.749439) + 
    (2.70 <= eta && eta < 2.80) * (34.30 <= pt * cosh(eta) && pt * cosh(eta) < 36.10) * (0.728180) + 
    (2.70 <= eta && eta < 2.80) * (36.10 <= pt * cosh(eta) && pt * cosh(eta) < 38.30) * (0.706664) + 
    (2.70 <= eta && eta < 2.80) * (38.30 <= pt * cosh(eta) && pt * cosh(eta) < 41.10) * (0.685265) + 
    (2.70 <= eta && eta < 2.80) * (41.10 <= pt * cosh(eta) && pt * cosh(eta) < 45.20) * (0.663498) + 
    (2.70 <= eta && eta < 2.80) * (45.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.644379) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (2.80 <= eta && eta < 2.90) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 20.00) * (0.999345) + 
    (2.80 <= eta && eta < 2.90) * (20.00 <= pt * cosh(eta) && pt * cosh(eta) < 21.70) * (0.980177) + 
    (2.80 <= eta && eta < 2.90) * (21.70 <= pt * cosh(eta) && pt * cosh(eta) < 22.90) * (0.959986) + 
    (2.80 <= eta && eta < 2.90) * (22.90 <= pt * cosh(eta) && pt * cosh(eta) < 24.00) * (0.939051) + 
    (2.80 <= eta && eta < 2.90) * (24.00 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.917427) + 
    (2.80 <= eta && eta < 2.90) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 26.00) * (0.895618) + 
    (2.80 <= eta && eta < 2.90) * (26.00 <= pt * cosh(eta) && pt * cosh(eta) < 27.00) * (0.873476) + 
    (2.80 <= eta && eta < 2.90) * (27.00 <= pt * cosh(eta) && pt * cosh(eta) < 28.00) * (0.851659) + 
    (2.80 <= eta && eta < 2.90) * (28.00 <= pt * cosh(eta) && pt * cosh(eta) < 29.00) * (0.830649) + 
    (2.80 <= eta && eta < 2.90) * (29.00 <= pt * cosh(eta) && pt * cosh(eta) < 30.10) * (0.809816) + 
    (2.80 <= eta && eta < 2.90) * (30.10 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.788676) + 
    (2.80 <= eta && eta < 2.90) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 32.60) * (0.767814) + 
    (2.80 <= eta && eta < 2.90) * (32.60 <= pt * cosh(eta) && pt * cosh(eta) < 34.10) * (0.747063) + 
    (2.80 <= eta && eta < 2.90) * (34.10 <= pt * cosh(eta) && pt * cosh(eta) < 35.90) * (0.725939) + 
    (2.80 <= eta && eta < 2.90) * (35.90 <= pt * cosh(eta) && pt * cosh(eta) < 38.10) * (0.704618) + 
    (2.80 <= eta && eta < 2.90) * (38.10 <= pt * cosh(eta) && pt * cosh(eta) < 41.00) * (0.683120) + 
    (2.80 <= eta && eta < 2.90) * (41.00 <= pt * cosh(eta) && pt * cosh(eta) < 45.20) * (0.661271) + 
    (2.80 <= eta && eta < 2.90) * (45.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.642594) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (2.90 <= eta && eta < 3.00) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.80) * (0.999346) + 
    (2.90 <= eta && eta < 3.00) * (19.80 <= pt * cosh(eta) && pt * cosh(eta) < 21.50) * (0.980111) + 
    (2.90 <= eta && eta < 3.00) * (21.50 <= pt * cosh(eta) && pt * cosh(eta) < 22.70) * (0.959654) + 
    (2.90 <= eta && eta < 3.00) * (22.70 <= pt * cosh(eta) && pt * cosh(eta) < 23.70) * (0.939441) + 
    (2.90 <= eta && eta < 3.00) * (23.70 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.918685) + 
    (2.90 <= eta && eta < 3.00) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.897815) + 
    (2.90 <= eta && eta < 3.00) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.50) * (0.877680) + 
    (2.90 <= eta && eta < 3.00) * (26.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.856635) + 
    (2.90 <= eta && eta < 3.00) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.835196) + 
    (2.90 <= eta && eta < 3.00) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.60) * (0.813874) + 
    (2.90 <= eta && eta < 3.00) * (29.60 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.792195) + 
    (2.90 <= eta && eta < 3.00) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.10) * (0.770777) + 
    (2.90 <= eta && eta < 3.00) * (32.10 <= pt * cosh(eta) && pt * cosh(eta) < 33.60) * (0.749470) + 
    (2.90 <= eta && eta < 3.00) * (33.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.30) * (0.728378) + 
    (2.90 <= eta && eta < 3.00) * (35.30 <= pt * cosh(eta) && pt * cosh(eta) < 37.40) * (0.707406) + 
    (2.90 <= eta && eta < 3.00) * (37.40 <= pt * cosh(eta) && pt * cosh(eta) < 40.10) * (0.686211) + 
    (2.90 <= eta && eta < 3.00) * (40.10 <= pt * cosh(eta) && pt * cosh(eta) < 44.00) * (0.664582) + 
    (2.90 <= eta && eta < 3.00) * (44.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.643207) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (3.00 <= eta && eta < 3.10) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.60) * (0.999361) + 
    (3.00 <= eta && eta < 3.10) * (19.60 <= pt * cosh(eta) && pt * cosh(eta) < 21.30) * (0.980308) + 
    (3.00 <= eta && eta < 3.10) * (21.30 <= pt * cosh(eta) && pt * cosh(eta) < 22.50) * (0.959730) + 
    (3.00 <= eta && eta < 3.10) * (22.50 <= pt * cosh(eta) && pt * cosh(eta) < 23.50) * (0.939334) + 
    (3.00 <= eta && eta < 3.10) * (23.50 <= pt * cosh(eta) && pt * cosh(eta) < 24.50) * (0.918374) + 
    (3.00 <= eta && eta < 3.10) * (24.50 <= pt * cosh(eta) && pt * cosh(eta) < 25.40) * (0.897301) + 
    (3.00 <= eta && eta < 3.10) * (25.40 <= pt * cosh(eta) && pt * cosh(eta) < 26.30) * (0.876983) + 
    (3.00 <= eta && eta < 3.10) * (26.30 <= pt * cosh(eta) && pt * cosh(eta) < 27.30) * (0.855766) + 
    (3.00 <= eta && eta < 3.10) * (27.30 <= pt * cosh(eta) && pt * cosh(eta) < 28.30) * (0.834175) + 
    (3.00 <= eta && eta < 3.10) * (28.30 <= pt * cosh(eta) && pt * cosh(eta) < 29.40) * (0.812730) + 
    (3.00 <= eta && eta < 3.10) * (29.40 <= pt * cosh(eta) && pt * cosh(eta) < 30.60) * (0.790958) + 
    (3.00 <= eta && eta < 3.10) * (30.60 <= pt * cosh(eta) && pt * cosh(eta) < 31.90) * (0.769483) + 
    (3.00 <= eta && eta < 3.10) * (31.90 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.748155) + 
    (3.00 <= eta && eta < 3.10) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 35.10) * (0.727081) + 
    (3.00 <= eta && eta < 3.10) * (35.10 <= pt * cosh(eta) && pt * cosh(eta) < 37.20) * (0.706169) + 
    (3.00 <= eta && eta < 3.10) * (37.20 <= pt * cosh(eta) && pt * cosh(eta) < 40.00) * (0.684714) + 
    (3.00 <= eta && eta < 3.10) * (40.00 <= pt * cosh(eta) && pt * cosh(eta) < 44.00) * (0.662815) + 
    (3.00 <= eta && eta < 3.10) * (44.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.641802) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (3.10 <= eta && eta < 3.20) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.50) * (0.999327) + 
    (3.10 <= eta && eta < 3.20) * (19.50 <= pt * cosh(eta) && pt * cosh(eta) < 21.10) * (0.980212) + 
    (3.10 <= eta && eta < 3.20) * (21.10 <= pt * cosh(eta) && pt * cosh(eta) < 22.30) * (0.960173) + 
    (3.10 <= eta && eta < 3.20) * (22.30 <= pt * cosh(eta) && pt * cosh(eta) < 23.30) * (0.939680) + 
    (3.10 <= eta && eta < 3.20) * (23.30 <= pt * cosh(eta) && pt * cosh(eta) < 24.30) * (0.918574) + 
    (3.10 <= eta && eta < 3.20) * (24.30 <= pt * cosh(eta) && pt * cosh(eta) < 25.20) * (0.897334) + 
    (3.10 <= eta && eta < 3.20) * (25.20 <= pt * cosh(eta) && pt * cosh(eta) < 26.10) * (0.876851) + 
    (3.10 <= eta && eta < 3.20) * (26.10 <= pt * cosh(eta) && pt * cosh(eta) < 27.10) * (0.855467) + 
    (3.10 <= eta && eta < 3.20) * (27.10 <= pt * cosh(eta) && pt * cosh(eta) < 28.10) * (0.833719) + 
    (3.10 <= eta && eta < 3.20) * (28.10 <= pt * cosh(eta) && pt * cosh(eta) < 29.20) * (0.812134) + 
    (3.10 <= eta && eta < 3.20) * (29.20 <= pt * cosh(eta) && pt * cosh(eta) < 30.40) * (0.790244) + 
    (3.10 <= eta && eta < 3.20) * (30.40 <= pt * cosh(eta) && pt * cosh(eta) < 31.70) * (0.768677) + 
    (3.10 <= eta && eta < 3.20) * (31.70 <= pt * cosh(eta) && pt * cosh(eta) < 33.20) * (0.747288) + 
    (3.10 <= eta && eta < 3.20) * (33.20 <= pt * cosh(eta) && pt * cosh(eta) < 34.90) * (0.726186) + 
    (3.10 <= eta && eta < 3.20) * (34.90 <= pt * cosh(eta) && pt * cosh(eta) < 37.00) * (0.705282) + 
    (3.10 <= eta && eta < 3.20) * (37.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.80) * (0.683878) + 
    (3.10 <= eta && eta < 3.20) * (39.80 <= pt * cosh(eta) && pt * cosh(eta) < 43.80) * (0.662083) + 
    (3.10 <= eta && eta < 3.20) * (43.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.640932) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (3.20 <= eta && eta < 3.30) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.30) * (0.999364) + 
    (3.20 <= eta && eta < 3.30) * (19.30 <= pt * cosh(eta) && pt * cosh(eta) < 20.90) * (0.980853) + 
    (3.20 <= eta && eta < 3.30) * (20.90 <= pt * cosh(eta) && pt * cosh(eta) < 22.10) * (0.960942) + 
    (3.20 <= eta && eta < 3.30) * (22.10 <= pt * cosh(eta) && pt * cosh(eta) < 23.10) * (0.940433) + 
    (3.20 <= eta && eta < 3.30) * (23.10 <= pt * cosh(eta) && pt * cosh(eta) < 24.10) * (0.919234) + 
    (3.20 <= eta && eta < 3.30) * (24.10 <= pt * cosh(eta) && pt * cosh(eta) < 25.00) * (0.897861) + 
    (3.20 <= eta && eta < 3.30) * (25.00 <= pt * cosh(eta) && pt * cosh(eta) < 25.90) * (0.877230) + 
    (3.20 <= eta && eta < 3.30) * (25.90 <= pt * cosh(eta) && pt * cosh(eta) < 26.80) * (0.856799) + 
    (3.20 <= eta && eta < 3.30) * (26.80 <= pt * cosh(eta) && pt * cosh(eta) < 27.80) * (0.835915) + 
    (3.20 <= eta && eta < 3.30) * (27.80 <= pt * cosh(eta) && pt * cosh(eta) < 28.90) * (0.814039) + 
    (3.20 <= eta && eta < 3.30) * (28.90 <= pt * cosh(eta) && pt * cosh(eta) < 30.00) * (0.792748) + 
    (3.20 <= eta && eta < 3.30) * (30.00 <= pt * cosh(eta) && pt * cosh(eta) < 31.30) * (0.771621) + 
    (3.20 <= eta && eta < 3.30) * (31.30 <= pt * cosh(eta) && pt * cosh(eta) < 32.70) * (0.750436) + 
    (3.20 <= eta && eta < 3.30) * (32.70 <= pt * cosh(eta) && pt * cosh(eta) < 34.40) * (0.729364) + 
    (3.20 <= eta && eta < 3.30) * (34.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.50) * (0.707748) + 
    (3.20 <= eta && eta < 3.30) * (36.50 <= pt * cosh(eta) && pt * cosh(eta) < 39.20) * (0.686021) + 
    (3.20 <= eta && eta < 3.30) * (39.20 <= pt * cosh(eta) && pt * cosh(eta) < 43.00) * (0.664270) + 
    (3.20 <= eta && eta < 3.30) * (43.00 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.641582) + 
    (3.20 <= eta && eta < 3.30) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.631269) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (3.30 <= eta && eta < 3.40) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.20) * (0.999349) + 
    (3.30 <= eta && eta < 3.40) * (19.20 <= pt * cosh(eta) && pt * cosh(eta) < 20.80) * (0.980513) + 
    (3.30 <= eta && eta < 3.40) * (20.80 <= pt * cosh(eta) && pt * cosh(eta) < 22.00) * (0.960291) + 
    (3.30 <= eta && eta < 3.40) * (22.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.00) * (0.939515) + 
    (3.30 <= eta && eta < 3.40) * (23.00 <= pt * cosh(eta) && pt * cosh(eta) < 23.90) * (0.919209) + 
    (3.30 <= eta && eta < 3.40) * (23.90 <= pt * cosh(eta) && pt * cosh(eta) < 24.80) * (0.898833) + 
    (3.30 <= eta && eta < 3.40) * (24.80 <= pt * cosh(eta) && pt * cosh(eta) < 25.70) * (0.878070) + 
    (3.30 <= eta && eta < 3.40) * (25.70 <= pt * cosh(eta) && pt * cosh(eta) < 26.60) * (0.857488) + 
    (3.30 <= eta && eta < 3.40) * (26.60 <= pt * cosh(eta) && pt * cosh(eta) < 27.60) * (0.836442) + 
    (3.30 <= eta && eta < 3.40) * (27.60 <= pt * cosh(eta) && pt * cosh(eta) < 28.60) * (0.815401) + 
    (3.30 <= eta && eta < 3.40) * (28.60 <= pt * cosh(eta) && pt * cosh(eta) < 29.70) * (0.794818) + 
    (3.30 <= eta && eta < 3.40) * (29.70 <= pt * cosh(eta) && pt * cosh(eta) < 30.90) * (0.774191) + 
    (3.30 <= eta && eta < 3.40) * (30.90 <= pt * cosh(eta) && pt * cosh(eta) < 32.30) * (0.753333) + 
    (3.30 <= eta && eta < 3.40) * (32.30 <= pt * cosh(eta) && pt * cosh(eta) < 33.90) * (0.732344) + 
    (3.30 <= eta && eta < 3.40) * (33.90 <= pt * cosh(eta) && pt * cosh(eta) < 35.90) * (0.711125) + 
    (3.30 <= eta && eta < 3.40) * (35.90 <= pt * cosh(eta) && pt * cosh(eta) < 38.40) * (0.689741) + 
    (3.30 <= eta && eta < 3.40) * (38.40 <= pt * cosh(eta) && pt * cosh(eta) < 41.90) * (0.668268) + 
    (3.30 <= eta && eta < 3.40) * (41.90 <= pt * cosh(eta) && pt * cosh(eta) < 47.80) * (0.645892) + 
    (3.30 <= eta && eta < 3.40) * (47.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.632957) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 2.60) * (0.601058) + 
    (3.40 <= eta && eta < 3.50) * (2.60 <= pt * cosh(eta) && pt * cosh(eta) < 19.10) * (0.999343) + 
    (3.40 <= eta && eta < 3.50) * (19.10 <= pt * cosh(eta) && pt * cosh(eta) < 20.70) * (0.980330) + 
    (3.40 <= eta && eta < 3.50) * (20.70 <= pt * cosh(eta) && pt * cosh(eta) < 21.90) * (0.959884) + 
    (3.40 <= eta && eta < 3.50) * (21.90 <= pt * cosh(eta) && pt * cosh(eta) < 22.90) * (0.938900) + 
    (3.40 <= eta && eta < 3.50) * (22.90 <= pt * cosh(eta) && pt * cosh(eta) < 23.80) * (0.918418) + 
    (3.40 <= eta && eta < 3.50) * (23.80 <= pt * cosh(eta) && pt * cosh(eta) < 24.70) * (0.897893) + 
    (3.40 <= eta && eta < 3.50) * (24.70 <= pt * cosh(eta) && pt * cosh(eta) < 25.60) * (0.877004) + 
    (3.40 <= eta && eta < 3.50) * (25.60 <= pt * cosh(eta) && pt * cosh(eta) < 26.50) * (0.856326) + 
    (3.40 <= eta && eta < 3.50) * (26.50 <= pt * cosh(eta) && pt * cosh(eta) < 27.50) * (0.835207) + 
    (3.40 <= eta && eta < 3.50) * (27.50 <= pt * cosh(eta) && pt * cosh(eta) < 28.50) * (0.814120) + 
    (3.40 <= eta && eta < 3.50) * (28.50 <= pt * cosh(eta) && pt * cosh(eta) < 29.60) * (0.793519) + 
    (3.40 <= eta && eta < 3.50) * (29.60 <= pt * cosh(eta) && pt * cosh(eta) < 30.80) * (0.772901) + 
    (3.40 <= eta && eta < 3.50) * (30.80 <= pt * cosh(eta) && pt * cosh(eta) < 32.20) * (0.752078) + 
    (3.40 <= eta && eta < 3.50) * (32.20 <= pt * cosh(eta) && pt * cosh(eta) < 33.80) * (0.731152) + 
    (3.40 <= eta && eta < 3.50) * (33.80 <= pt * cosh(eta) && pt * cosh(eta) < 35.80) * (0.710027) + 
    (3.40 <= eta && eta < 3.50) * (35.80 <= pt * cosh(eta) && pt * cosh(eta) < 38.40) * (0.688369) + 
    (3.40 <= eta && eta < 3.50) * (38.40 <= pt * cosh(eta) && pt * cosh(eta) < 42.00) * (0.666571) + 
    (3.40 <= eta && eta < 3.50) * (42.00 <= pt * cosh(eta) && pt * cosh(eta) < 48.20) * (0.644053) + 
    (3.40 <= eta && eta < 3.50) * (48.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.631748) } 

# Everything else is not ID'd at all.
  add EfficiencyFormula {0} {0} { 0.00 }

}


#My name is "dualRICH_c2f6" and I am described as follows: 
#   Detector Type =  1 [Barrel=0, Forward=1]
#    Eta coverage =  [1,4]
#           Radii =  [9.1608925814663991,212.7295320598304]
#               Z =  250

module IdentificationMap dualRICH_c2f6 { 
  set InputArray TrackSmearing/tracks
  set OutputArray tracks 

  # --- kaons ---

  add EfficiencyFormula {321} {321} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.00 <= eta && eta < 1.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.999386) + 
    (1.00 <= eta && eta < 1.10) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.10) * (0.980825) + 
    (1.00 <= eta && eta < 1.10) * (36.10 <= pt * cosh(eta) && pt * cosh(eta) < 38.10) * (0.961136) + 
    (1.00 <= eta && eta < 1.10) * (38.10 <= pt * cosh(eta) && pt * cosh(eta) < 39.80) * (0.940849) + 
    (1.00 <= eta && eta < 1.10) * (39.80 <= pt * cosh(eta) && pt * cosh(eta) < 41.40) * (0.920315) + 
    (1.00 <= eta && eta < 1.10) * (41.40 <= pt * cosh(eta) && pt * cosh(eta) < 42.90) * (0.899839) + 
    (1.00 <= eta && eta < 1.10) * (42.90 <= pt * cosh(eta) && pt * cosh(eta) < 44.40) * (0.879620) + 
    (1.00 <= eta && eta < 1.10) * (44.40 <= pt * cosh(eta) && pt * cosh(eta) < 46.00) * (0.858889) + 
    (1.00 <= eta && eta < 1.10) * (46.00 <= pt * cosh(eta) && pt * cosh(eta) < 47.60) * (0.838125) + 
    (1.00 <= eta && eta < 1.10) * (47.60 <= pt * cosh(eta) && pt * cosh(eta) < 49.30) * (0.817751) + 
    (1.00 <= eta && eta < 1.10) * (49.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.803676) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.10 <= eta && eta < 1.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.999404) + 
    (1.10 <= eta && eta < 1.20) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 38.10) * (0.981046) + 
    (1.10 <= eta && eta < 1.20) * (38.10 <= pt * cosh(eta) && pt * cosh(eta) < 40.20) * (0.961386) + 
    (1.10 <= eta && eta < 1.20) * (40.20 <= pt * cosh(eta) && pt * cosh(eta) < 42.00) * (0.941250) + 
    (1.10 <= eta && eta < 1.20) * (42.00 <= pt * cosh(eta) && pt * cosh(eta) < 43.70) * (0.920721) + 
    (1.10 <= eta && eta < 1.20) * (43.70 <= pt * cosh(eta) && pt * cosh(eta) < 45.30) * (0.900151) + 
    (1.10 <= eta && eta < 1.20) * (45.30 <= pt * cosh(eta) && pt * cosh(eta) < 46.90) * (0.879783) + 
    (1.10 <= eta && eta < 1.20) * (46.90 <= pt * cosh(eta) && pt * cosh(eta) < 48.60) * (0.858938) + 
    (1.10 <= eta && eta < 1.20) * (48.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.839889) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.20 <= eta && eta < 1.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 37.20) * (0.999403) + 
    (1.20 <= eta && eta < 1.30) * (37.20 <= pt * cosh(eta) && pt * cosh(eta) < 40.30) * (0.980937) + 
    (1.20 <= eta && eta < 1.30) * (40.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.50) * (0.961246) + 
    (1.20 <= eta && eta < 1.30) * (42.50 <= pt * cosh(eta) && pt * cosh(eta) < 44.40) * (0.941260) + 
    (1.20 <= eta && eta < 1.30) * (44.40 <= pt * cosh(eta) && pt * cosh(eta) < 46.20) * (0.920785) + 
    (1.20 <= eta && eta < 1.30) * (46.20 <= pt * cosh(eta) && pt * cosh(eta) < 47.90) * (0.900200) + 
    (1.20 <= eta && eta < 1.30) * (47.90 <= pt * cosh(eta) && pt * cosh(eta) < 49.60) * (0.879778) + 
    (1.20 <= eta && eta < 1.30) * (49.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.867191) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.30 <= eta && eta < 1.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.999386) + 
    (1.30 <= eta && eta < 1.40) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 42.60) * (0.980852) + 
    (1.30 <= eta && eta < 1.40) * (42.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.961166) + 
    (1.30 <= eta && eta < 1.40) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 47.00) * (0.940918) + 
    (1.30 <= eta && eta < 1.40) * (47.00 <= pt * cosh(eta) && pt * cosh(eta) < 48.90) * (0.920537) + 
    (1.30 <= eta && eta < 1.40) * (48.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.903953) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.40 <= eta && eta < 1.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 41.60) * (0.999406) + 
    (1.40 <= eta && eta < 1.50) * (41.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.10) * (0.980930) + 
    (1.40 <= eta && eta < 1.50) * (45.10 <= pt * cosh(eta) && pt * cosh(eta) < 47.60) * (0.961086) + 
    (1.40 <= eta && eta < 1.50) * (47.60 <= pt * cosh(eta) && pt * cosh(eta) < 49.70) * (0.941124) + 
    (1.40 <= eta && eta < 1.50) * (49.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.929574) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.50 <= eta && eta < 1.60) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 44.00) * (0.999402) + 
    (1.50 <= eta && eta < 1.60) * (44.00 <= pt * cosh(eta) && pt * cosh(eta) < 47.70) * (0.980852) + 
    (1.50 <= eta && eta < 1.60) * (47.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.962344) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.60 <= eta && eta < 1.70) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.999392) + 
    (1.60 <= eta && eta < 1.70) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.981712) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.70 <= eta && eta < 1.80) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 48.90) * (0.999398) + 
    (1.70 <= eta && eta < 1.80) * (48.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.987334) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.80 <= eta && eta < 1.90) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999578) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.90 <= eta && eta < 2.00) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999802) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.00 <= eta && eta < 2.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999889) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.10 <= eta && eta < 2.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999894) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.20 <= eta && eta < 2.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999898) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.30 <= eta && eta < 2.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999899) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.40 <= eta && eta < 2.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999886) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.50 <= eta && eta < 2.60) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999725) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.60 <= eta && eta < 2.70) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999448) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.70 <= eta && eta < 2.80) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 48.70) * (0.999390) + 
    (2.70 <= eta && eta < 2.80) * (48.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.986807) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.80 <= eta && eta < 2.90) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.30) * (0.999400) + 
    (2.80 <= eta && eta < 2.90) * (47.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.983941) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.90 <= eta && eta < 3.00) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 46.20) * (0.999386) + 
    (2.90 <= eta && eta < 3.00) * (46.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.980746) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (3.00 <= eta && eta < 3.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.20) * (0.999386) + 
    (3.00 <= eta && eta < 3.10) * (45.20 <= pt * cosh(eta) && pt * cosh(eta) < 48.90) * (0.980797) + 
    (3.00 <= eta && eta < 3.10) * (48.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.967085) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (3.10 <= eta && eta < 3.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 44.30) * (0.999395) + 
    (3.10 <= eta && eta < 3.20) * (44.30 <= pt * cosh(eta) && pt * cosh(eta) < 48.00) * (0.980787) + 
    (3.10 <= eta && eta < 3.20) * (48.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.963510) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (3.20 <= eta && eta < 3.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 43.60) * (0.999383) + 
    (3.20 <= eta && eta < 3.30) * (43.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.980645) + 
    (3.20 <= eta && eta < 3.30) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.960976) + 
    (3.20 <= eta && eta < 3.30) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.949897) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (3.30 <= eta && eta < 3.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.90) * (0.999396) + 
    (3.30 <= eta && eta < 3.40) * (42.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.980744) + 
    (3.30 <= eta && eta < 3.40) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 49.10) * (0.960764) + 
    (3.30 <= eta && eta < 3.40) * (49.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.946412) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (3.40 <= eta && eta < 3.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.999404) + 
    (3.40 <= eta && eta < 3.50) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.981060) + 
    (3.40 <= eta && eta < 3.50) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 48.30) * (0.961665) + 
    (3.40 <= eta && eta < 3.50) * (48.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.944030) } 

  add EfficiencyFormula {321} {-211} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.00 <= eta && eta < 1.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 36.00) * (0.000710) + 
    (1.00 <= eta && eta < 1.10) * (36.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.60) * (0.019224) + 
    (1.00 <= eta && eta < 1.10) * (39.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.039044) + 
    (1.00 <= eta && eta < 1.10) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 44.60) * (0.059154) + 
    (1.00 <= eta && eta < 1.10) * (44.60 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.079361) + 
    (1.00 <= eta && eta < 1.10) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 48.90) * (0.099640) + 
    (1.00 <= eta && eta < 1.10) * (48.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.114949) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.10 <= eta && eta < 1.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.000706) + 
    (1.10 <= eta && eta < 1.20) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 41.80) * (0.019135) + 
    (1.10 <= eta && eta < 1.20) * (41.80 <= pt * cosh(eta) && pt * cosh(eta) < 44.60) * (0.038662) + 
    (1.10 <= eta && eta < 1.20) * (44.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.10) * (0.058743) + 
    (1.10 <= eta && eta < 1.20) * (47.10 <= pt * cosh(eta) && pt * cosh(eta) < 49.40) * (0.079070) + 
    (1.10 <= eta && eta < 1.20) * (49.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.091902) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.20 <= eta && eta < 1.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 40.20) * (0.000714) + 
    (1.20 <= eta && eta < 1.30) * (40.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.20) * (0.019228) + 
    (1.20 <= eta && eta < 1.30) * (44.20 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.038821) + 
    (1.20 <= eta && eta < 1.30) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.058859) + 
    (1.20 <= eta && eta < 1.30) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.069819) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.30 <= eta && eta < 1.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.50) * (0.000708) + 
    (1.30 <= eta && eta < 1.40) * (42.50 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.019283) + 
    (1.30 <= eta && eta < 1.40) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.039127) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.40 <= eta && eta < 1.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.000717) + 
    (1.40 <= eta && eta < 1.50) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 49.50) * (0.019313) + 
    (1.40 <= eta && eta < 1.50) * (49.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.030744) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.50 <= eta && eta < 1.60) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.60) * (0.000720) + 
    (1.50 <= eta && eta < 1.60) * (47.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.014769) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.60 <= eta && eta < 1.70) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000663) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.70 <= eta && eta < 1.80) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000282) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.80 <= eta && eta < 1.90) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000127) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.90 <= eta && eta < 2.00) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000056) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.00 <= eta && eta < 2.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000030) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.10 <= eta && eta < 2.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000029) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.20 <= eta && eta < 2.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000028) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.30 <= eta && eta < 2.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000027) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.40 <= eta && eta < 2.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000031) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.50 <= eta && eta < 2.60) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000080) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.60 <= eta && eta < 2.70) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000169) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.70 <= eta && eta < 2.80) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000308) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.80 <= eta && eta < 2.90) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000499) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.90 <= eta && eta < 3.00) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 49.90) * (0.000716) + 
    (2.90 <= eta && eta < 3.00) * (49.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.010917) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (3.00 <= eta && eta < 3.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 48.80) * (0.000712) + 
    (3.00 <= eta && eta < 3.10) * (48.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.012607) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (3.10 <= eta && eta < 3.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.90) * (0.000720) + 
    (3.10 <= eta && eta < 3.20) * (47.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.014230) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (3.20 <= eta && eta < 3.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.00) * (0.000700) + 
    (3.20 <= eta && eta < 3.30) * (47.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.015613) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (3.30 <= eta && eta < 3.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 46.30) * (0.000701) + 
    (3.30 <= eta && eta < 3.40) * (46.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.017086) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (3.40 <= eta && eta < 3.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.000704) + 
    (3.40 <= eta && eta < 3.50) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.018507) } 

  add EfficiencyFormula {321} {2212} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.00 <= eta && eta < 1.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000012) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.10 <= eta && eta < 1.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000002) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.20 <= eta && eta < 1.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.30 <= eta && eta < 1.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.40 <= eta && eta < 1.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.50 <= eta && eta < 1.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.60 <= eta && eta < 1.70) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.70 <= eta && eta < 1.80) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.80 <= eta && eta < 1.90) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.90 <= eta && eta < 2.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.00 <= eta && eta < 2.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.10 <= eta && eta < 2.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.20 <= eta && eta < 2.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.30 <= eta && eta < 2.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.40 <= eta && eta < 2.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.50 <= eta && eta < 2.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.60 <= eta && eta < 2.70) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.70 <= eta && eta < 2.80) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.80 <= eta && eta < 2.90) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.90 <= eta && eta < 3.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (3.00 <= eta && eta < 3.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (3.10 <= eta && eta < 3.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (3.20 <= eta && eta < 3.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (3.30 <= eta && eta < 3.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (3.40 <= eta && eta < 3.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 


  # --- pions ---

  add EfficiencyFormula {-211} {321} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.00 <= eta && eta < 1.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 36.00) * (0.000710) + 
    (1.00 <= eta && eta < 1.10) * (36.00 <= pt * cosh(eta) && pt * cosh(eta) < 39.60) * (0.019224) + 
    (1.00 <= eta && eta < 1.10) * (39.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.039044) + 
    (1.00 <= eta && eta < 1.10) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 44.60) * (0.059154) + 
    (1.00 <= eta && eta < 1.10) * (44.60 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.079361) + 
    (1.00 <= eta && eta < 1.10) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 48.90) * (0.099640) + 
    (1.00 <= eta && eta < 1.10) * (48.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.114949) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.10 <= eta && eta < 1.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 38.00) * (0.000706) + 
    (1.10 <= eta && eta < 1.20) * (38.00 <= pt * cosh(eta) && pt * cosh(eta) < 41.80) * (0.019135) + 
    (1.10 <= eta && eta < 1.20) * (41.80 <= pt * cosh(eta) && pt * cosh(eta) < 44.60) * (0.038662) + 
    (1.10 <= eta && eta < 1.20) * (44.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.10) * (0.058743) + 
    (1.10 <= eta && eta < 1.20) * (47.10 <= pt * cosh(eta) && pt * cosh(eta) < 49.40) * (0.079070) + 
    (1.10 <= eta && eta < 1.20) * (49.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.091902) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.20 <= eta && eta < 1.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 40.20) * (0.000714) + 
    (1.20 <= eta && eta < 1.30) * (40.20 <= pt * cosh(eta) && pt * cosh(eta) < 44.20) * (0.019228) + 
    (1.20 <= eta && eta < 1.30) * (44.20 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.038821) + 
    (1.20 <= eta && eta < 1.30) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.058859) + 
    (1.20 <= eta && eta < 1.30) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.069819) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.30 <= eta && eta < 1.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.50) * (0.000708) + 
    (1.30 <= eta && eta < 1.40) * (42.50 <= pt * cosh(eta) && pt * cosh(eta) < 46.80) * (0.019283) + 
    (1.30 <= eta && eta < 1.40) * (46.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.039127) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.40 <= eta && eta < 1.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.000717) + 
    (1.40 <= eta && eta < 1.50) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 49.50) * (0.019313) + 
    (1.40 <= eta && eta < 1.50) * (49.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.030744) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.50 <= eta && eta < 1.60) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.60) * (0.000720) + 
    (1.50 <= eta && eta < 1.60) * (47.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.014769) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.60 <= eta && eta < 1.70) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000663) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.70 <= eta && eta < 1.80) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000282) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.80 <= eta && eta < 1.90) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000127) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (1.90 <= eta && eta < 2.00) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000056) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.00 <= eta && eta < 2.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000030) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.10 <= eta && eta < 2.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000029) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.20 <= eta && eta < 2.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000028) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.30 <= eta && eta < 2.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000027) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.40 <= eta && eta < 2.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000031) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.50 <= eta && eta < 2.60) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000080) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.60 <= eta && eta < 2.70) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000169) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.70 <= eta && eta < 2.80) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000308) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.80 <= eta && eta < 2.90) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000499) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (2.90 <= eta && eta < 3.00) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 49.90) * (0.000716) + 
    (2.90 <= eta && eta < 3.00) * (49.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.010917) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (3.00 <= eta && eta < 3.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 48.80) * (0.000712) + 
    (3.00 <= eta && eta < 3.10) * (48.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.012607) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (3.10 <= eta && eta < 3.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.90) * (0.000720) + 
    (3.10 <= eta && eta < 3.20) * (47.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.014230) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (3.20 <= eta && eta < 3.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.00) * (0.000700) + 
    (3.20 <= eta && eta < 3.30) * (47.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.015613) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (3.30 <= eta && eta < 3.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 46.30) * (0.000701) + 
    (3.30 <= eta && eta < 3.40) * (46.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.017086) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.500000) + 
    (3.40 <= eta && eta < 3.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.70) * (0.000704) + 
    (3.40 <= eta && eta < 3.50) * (45.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.018507) } 

  add EfficiencyFormula {-211} {-211} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.00 <= eta && eta < 1.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 33.40) * (0.999386) + 
    (1.00 <= eta && eta < 1.10) * (33.40 <= pt * cosh(eta) && pt * cosh(eta) < 36.10) * (0.980825) + 
    (1.00 <= eta && eta < 1.10) * (36.10 <= pt * cosh(eta) && pt * cosh(eta) < 38.10) * (0.961136) + 
    (1.00 <= eta && eta < 1.10) * (38.10 <= pt * cosh(eta) && pt * cosh(eta) < 39.80) * (0.940849) + 
    (1.00 <= eta && eta < 1.10) * (39.80 <= pt * cosh(eta) && pt * cosh(eta) < 41.40) * (0.920315) + 
    (1.00 <= eta && eta < 1.10) * (41.40 <= pt * cosh(eta) && pt * cosh(eta) < 42.90) * (0.899839) + 
    (1.00 <= eta && eta < 1.10) * (42.90 <= pt * cosh(eta) && pt * cosh(eta) < 44.40) * (0.879620) + 
    (1.00 <= eta && eta < 1.10) * (44.40 <= pt * cosh(eta) && pt * cosh(eta) < 46.00) * (0.858889) + 
    (1.00 <= eta && eta < 1.10) * (46.00 <= pt * cosh(eta) && pt * cosh(eta) < 47.60) * (0.838125) + 
    (1.00 <= eta && eta < 1.10) * (47.60 <= pt * cosh(eta) && pt * cosh(eta) < 49.30) * (0.817751) + 
    (1.00 <= eta && eta < 1.10) * (49.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.803676) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.10 <= eta && eta < 1.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 35.20) * (0.999404) + 
    (1.10 <= eta && eta < 1.20) * (35.20 <= pt * cosh(eta) && pt * cosh(eta) < 38.10) * (0.981046) + 
    (1.10 <= eta && eta < 1.20) * (38.10 <= pt * cosh(eta) && pt * cosh(eta) < 40.20) * (0.961386) + 
    (1.10 <= eta && eta < 1.20) * (40.20 <= pt * cosh(eta) && pt * cosh(eta) < 42.00) * (0.941250) + 
    (1.10 <= eta && eta < 1.20) * (42.00 <= pt * cosh(eta) && pt * cosh(eta) < 43.70) * (0.920721) + 
    (1.10 <= eta && eta < 1.20) * (43.70 <= pt * cosh(eta) && pt * cosh(eta) < 45.30) * (0.900151) + 
    (1.10 <= eta && eta < 1.20) * (45.30 <= pt * cosh(eta) && pt * cosh(eta) < 46.90) * (0.879783) + 
    (1.10 <= eta && eta < 1.20) * (46.90 <= pt * cosh(eta) && pt * cosh(eta) < 48.60) * (0.858938) + 
    (1.10 <= eta && eta < 1.20) * (48.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.839889) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.20 <= eta && eta < 1.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 37.20) * (0.999403) + 
    (1.20 <= eta && eta < 1.30) * (37.20 <= pt * cosh(eta) && pt * cosh(eta) < 40.30) * (0.980937) + 
    (1.20 <= eta && eta < 1.30) * (40.30 <= pt * cosh(eta) && pt * cosh(eta) < 42.50) * (0.961246) + 
    (1.20 <= eta && eta < 1.30) * (42.50 <= pt * cosh(eta) && pt * cosh(eta) < 44.40) * (0.941260) + 
    (1.20 <= eta && eta < 1.30) * (44.40 <= pt * cosh(eta) && pt * cosh(eta) < 46.20) * (0.920785) + 
    (1.20 <= eta && eta < 1.30) * (46.20 <= pt * cosh(eta) && pt * cosh(eta) < 47.90) * (0.900200) + 
    (1.20 <= eta && eta < 1.30) * (47.90 <= pt * cosh(eta) && pt * cosh(eta) < 49.60) * (0.879778) + 
    (1.20 <= eta && eta < 1.30) * (49.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.867191) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.30 <= eta && eta < 1.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 39.40) * (0.999386) + 
    (1.30 <= eta && eta < 1.40) * (39.40 <= pt * cosh(eta) && pt * cosh(eta) < 42.60) * (0.980852) + 
    (1.30 <= eta && eta < 1.40) * (42.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.00) * (0.961166) + 
    (1.30 <= eta && eta < 1.40) * (45.00 <= pt * cosh(eta) && pt * cosh(eta) < 47.00) * (0.940918) + 
    (1.30 <= eta && eta < 1.40) * (47.00 <= pt * cosh(eta) && pt * cosh(eta) < 48.90) * (0.920537) + 
    (1.30 <= eta && eta < 1.40) * (48.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.903953) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.40 <= eta && eta < 1.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 41.60) * (0.999406) + 
    (1.40 <= eta && eta < 1.50) * (41.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.10) * (0.980930) + 
    (1.40 <= eta && eta < 1.50) * (45.10 <= pt * cosh(eta) && pt * cosh(eta) < 47.60) * (0.961086) + 
    (1.40 <= eta && eta < 1.50) * (47.60 <= pt * cosh(eta) && pt * cosh(eta) < 49.70) * (0.941124) + 
    (1.40 <= eta && eta < 1.50) * (49.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.929574) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.50 <= eta && eta < 1.60) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 44.00) * (0.999402) + 
    (1.50 <= eta && eta < 1.60) * (44.00 <= pt * cosh(eta) && pt * cosh(eta) < 47.70) * (0.980852) + 
    (1.50 <= eta && eta < 1.60) * (47.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.962344) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.60 <= eta && eta < 1.70) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.999392) + 
    (1.60 <= eta && eta < 1.70) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.981712) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.70 <= eta && eta < 1.80) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 48.90) * (0.999398) + 
    (1.70 <= eta && eta < 1.80) * (48.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.987334) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.80 <= eta && eta < 1.90) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999578) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (1.90 <= eta && eta < 2.00) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999802) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.00 <= eta && eta < 2.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999889) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.10 <= eta && eta < 2.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999894) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.20 <= eta && eta < 2.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999898) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.30 <= eta && eta < 2.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999899) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.40 <= eta && eta < 2.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999886) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.50 <= eta && eta < 2.60) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999725) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.60 <= eta && eta < 2.70) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999448) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.70 <= eta && eta < 2.80) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 48.70) * (0.999390) + 
    (2.70 <= eta && eta < 2.80) * (48.70 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.986807) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.80 <= eta && eta < 2.90) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.30) * (0.999400) + 
    (2.80 <= eta && eta < 2.90) * (47.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.983941) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (2.90 <= eta && eta < 3.00) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 46.20) * (0.999386) + 
    (2.90 <= eta && eta < 3.00) * (46.20 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.980746) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (3.00 <= eta && eta < 3.10) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 45.20) * (0.999386) + 
    (3.00 <= eta && eta < 3.10) * (45.20 <= pt * cosh(eta) && pt * cosh(eta) < 48.90) * (0.980797) + 
    (3.00 <= eta && eta < 3.10) * (48.90 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.967085) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (3.10 <= eta && eta < 3.20) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 44.30) * (0.999395) + 
    (3.10 <= eta && eta < 3.20) * (44.30 <= pt * cosh(eta) && pt * cosh(eta) < 48.00) * (0.980787) + 
    (3.10 <= eta && eta < 3.20) * (48.00 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.963510) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (3.20 <= eta && eta < 3.30) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 43.60) * (0.999383) + 
    (3.20 <= eta && eta < 3.30) * (43.60 <= pt * cosh(eta) && pt * cosh(eta) < 47.20) * (0.980645) + 
    (3.20 <= eta && eta < 3.30) * (47.20 <= pt * cosh(eta) && pt * cosh(eta) < 49.80) * (0.960976) + 
    (3.20 <= eta && eta < 3.30) * (49.80 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.949897) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (3.30 <= eta && eta < 3.40) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.90) * (0.999396) + 
    (3.30 <= eta && eta < 3.40) * (42.90 <= pt * cosh(eta) && pt * cosh(eta) < 46.50) * (0.980744) + 
    (3.30 <= eta && eta < 3.40) * (46.50 <= pt * cosh(eta) && pt * cosh(eta) < 49.10) * (0.960764) + 
    (3.30 <= eta && eta < 3.40) * (49.10 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.946412) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 3.60) * (0.601058) + 
    (3.40 <= eta && eta < 3.50) * (3.60 <= pt * cosh(eta) && pt * cosh(eta) < 42.30) * (0.999404) + 
    (3.40 <= eta && eta < 3.50) * (42.30 <= pt * cosh(eta) && pt * cosh(eta) < 45.80) * (0.981060) + 
    (3.40 <= eta && eta < 3.50) * (45.80 <= pt * cosh(eta) && pt * cosh(eta) < 48.30) * (0.961665) + 
    (3.40 <= eta && eta < 3.50) * (48.30 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.944030) } 


  # --- protons ---

  add EfficiencyFormula {2212} {321} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.00 <= eta && eta < 1.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000012) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.10 <= eta && eta < 1.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000002) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.20 <= eta && eta < 1.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.30 <= eta && eta < 1.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.40 <= eta && eta < 1.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.50 <= eta && eta < 1.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.60 <= eta && eta < 1.70) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.70 <= eta && eta < 1.80) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.80 <= eta && eta < 1.90) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (1.90 <= eta && eta < 2.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.00 <= eta && eta < 2.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.10 <= eta && eta < 2.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.20 <= eta && eta < 2.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.30 <= eta && eta < 2.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.40 <= eta && eta < 2.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.50 <= eta && eta < 2.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.60 <= eta && eta < 2.70) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.70 <= eta && eta < 2.80) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.80 <= eta && eta < 2.90) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (2.90 <= eta && eta < 3.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (3.00 <= eta && eta < 3.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (3.10 <= eta && eta < 3.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (3.20 <= eta && eta < 3.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (3.30 <= eta && eta < 3.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.500000) + 
    (3.40 <= eta && eta < 3.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.000000) } 

  add EfficiencyFormula {2212} {2212} {
    (eta<1.00 || eta>=3.50 || pt * cosh(eta) < 0.10 || pt * cosh(eta) >= 50.00) * ( 0.00 ) + 
    (1.00 <= eta && eta < 1.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.00 <= eta && eta < 1.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (1.00 <= eta && eta < 1.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999953) + 
    (1.10 <= eta && eta < 1.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.10 <= eta && eta < 1.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (1.10 <= eta && eta < 1.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999991) + 
    (1.20 <= eta && eta < 1.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.20 <= eta && eta < 1.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (1.20 <= eta && eta < 1.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (0.999999) + 
    (1.30 <= eta && eta < 1.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.30 <= eta && eta < 1.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (1.30 <= eta && eta < 1.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.40 <= eta && eta < 1.50) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (1.40 <= eta && eta < 1.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.50 <= eta && eta < 1.60) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (1.50 <= eta && eta < 1.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.60 <= eta && eta < 1.70) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (1.60 <= eta && eta < 1.70) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.70 <= eta && eta < 1.80) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (1.70 <= eta && eta < 1.80) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.80 <= eta && eta < 1.90) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (1.80 <= eta && eta < 1.90) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (1.90 <= eta && eta < 2.00) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (1.90 <= eta && eta < 2.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.00 <= eta && eta < 2.10) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (2.00 <= eta && eta < 2.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.10 <= eta && eta < 2.20) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (2.10 <= eta && eta < 2.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.20 <= eta && eta < 2.30) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (2.20 <= eta && eta < 2.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.20) * (0.000000) + 
    (2.30 <= eta && eta < 2.40) * (0.20 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (2.30 <= eta && eta < 2.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.40 <= eta && eta < 2.50) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (2.40 <= eta && eta < 2.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.50 <= eta && eta < 2.60) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (2.50 <= eta && eta < 2.60) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.60 <= eta && eta < 2.70) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (2.60 <= eta && eta < 2.70) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.70 <= eta && eta < 2.80) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (2.70 <= eta && eta < 2.80) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.30) * (0.000000) + 
    (2.80 <= eta && eta < 2.90) * (0.30 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (2.80 <= eta && eta < 2.90) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (2.90 <= eta && eta < 3.00) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (2.90 <= eta && eta < 3.00) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.00 <= eta && eta < 3.10) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (3.00 <= eta && eta < 3.10) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.40) * (0.000000) + 
    (3.10 <= eta && eta < 3.20) * (0.40 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (3.10 <= eta && eta < 3.20) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.20 <= eta && eta < 3.30) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (3.20 <= eta && eta < 3.30) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.30 <= eta && eta < 3.40) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (3.30 <= eta && eta < 3.40) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.10 <= pt * cosh(eta) && pt * cosh(eta) < 0.50) * (0.000000) + 
    (3.40 <= eta && eta < 3.50) * (0.50 <= pt * cosh(eta) && pt * cosh(eta) < 12.40) * (0.601058) + 
    (3.40 <= eta && eta < 3.50) * (12.40 <= pt * cosh(eta) && pt * cosh(eta) < 50.00) * (1.000000) } 

# Everything else is not ID'd at all.
  add EfficiencyFormula {0} {0} { 0.00 }

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

  add Branch TrackSmearing/tracks Track Track
  add Branch Calorimeter/towers Tower Tower

  add Branch HCal/eflowTracks EFlowTrack Track
  add Branch ECal/eflowPhotons EFlowPhoton Tower
  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch mRICH/tracks mRICHTrack Track
  add Branch barrelDIRC/tracks barrelDIRCTrack Track
  add Branch dualRICH_aerogel/tracks dualRICHagTrack Track
  add Branch dualRICH_c2f6/tracks dualRICHcfTrack Track

  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenMissingET/momentum GenMissingET MissingET

  add Branch UniqueObjectFinder/jets Jet Jet
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/photons Photon Photon

  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}
