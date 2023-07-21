######################################################################################################################                                
# EIC detector model                                                                                                                                                                                     
# based on parameters from EIC detector matrix from EIC yellow report https://physdiv.jlab.org/DetectorMatrix/ (also in https://arxiv.org/abs/2103.05419). 
# as well as on assumptions on calorimeter granularity, and efficiency.
# email: miguel.arratia@ucr.edu, ssekula@smu.edu 
#######################################################################################################################


#######################################
# Load any external configurations
#######################################

#source customizations.tcl


#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency

 
  ChargedHadronMomentumSmearing
  ElectronMomentumSmearing
  

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

  PIDSystems

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
   
    # radius of the magnetic field coverage, in m
    set Radius 1.5
    # half-length of the magnetic field coverage, in m
    set HalfLength 1.20
    # magnetic field
    set Bz 3.0
}


####################################
# Common Tracking Efficiency Model
####################################
#Dummy efficiency (100%). Leaving structure to show how tracking dependent on pt and eta can be incorporated)
#

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
    (abs(eta) > 3.5)                                                  * (0.00)+
    0.0    
}

set CommonTrackingResolution {
    (abs(eta) <= 1.0)                        * sqrt((5e-3)^2 + (pt*cosh(eta))^2*(2e-4)^2) +
    (abs(eta) > 1.0 && abs(eta) <= 2.5)      * sqrt((1e-2)^2 + (pt*cosh(eta))^2*(2e-4)^2) +
    (abs(eta) > 2.5 && abs(eta) <= 3.5)      * sqrt((2e-2)^2 + (pt*cosh(eta))^2*(1e-3)^2) +
    (abs(eta) > 3.5)                                                  * (0.00)
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



########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set OutputArray chargedHadrons
  set ResolutionFormula  $CommonTrackingResolution
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
  set D0ResolutionFormula "( abs(eta) <= 2.5 && abs(eta)>1.0 ) * (sqrt( (0.010)^2 + (0.040/pt)^2)) +
                             ( abs(eta) <= 1.0 ) *                 (sqrt( (0.005)^2 + (0.030/pt)^2))"
  set DZResolutionFormula "( abs(eta) <= 2.5 && abs(eta)>1.0 ) *   (sqrt( (0.020)^2 + (0.100/pt)^2)) + 
                             ( abs(eta) <= 1.0 )               *   (sqrt( (0.005)^2 + (0.030/pt)^2)) "
    
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
# PID Efficiency Maps
##################

# These maps are built from the table in section 8.6 of the EIC Yellow Report.
# When the table says that two species are separated by "3 sigma" we assume the 
# following:
#
# 1. Each specie is Gaussian distributed in a variable x each with the same Gaussian width, s.
# 2. The separation is defined by S. For 3 sigma, S=3.
# 3. The separation takes into account the two widths of the individual Gaussians, S = sqrt(s^2 + s^2). 
#    Thus, s = S/sqrt(2) for each of the presumed Gaussians for each specie.
# 4. To compute the identification efficiency, the probability that specie A is idenfied as A, one
#    uses the integral of the Gaussian within width s, e.g. p = 1 - ROOT::Math::gaussian_pdf(s) in ROOT.
# 5. TO instead compute the misidentfication probability, the probability that B -> A, one instead uses
#    the one-sided cumulative distribution function, e.g. p = 1 - ROOT::Math::normal_cdf(s) in ROOT.
# EXAMPLE:
#
# Electrons need to be separable from pions at the level of 3 sigma for 1.0 <= eta <= 3.5. Thus S = 3.
# That means s = 3/sqrt(2) = 2.121. Thus p(e->e) = 1 - ROOT::Math::gaussian_pdf(2.121) = 0.958. 
# Conversely, p(pi -> e) = 1 - ROOT::Math::normal_cdf(2.121) = 0.017. Once can also compute from this
# p(e->pi) = 1 - ROOT::Math::normal_cdf(2.121) (assuming a 2-species model for PID) and p(pi->pi) = 0.958 (integral of gaussian within s)
#

module IdentificationMap PIDSystems { 
    set InputArray HCal/eflowTracks 
    set OutputArray tracks 
    
    
    # Electron/Pion identification, treated as a 2-species problem
    add EfficiencyFormula {-11} {-11} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.050) * (0.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.050) * (
							(-3.5 <= eta && eta < 1.0) * (1.00) +
							(1.0 <= eta && eta <= 3.5) * (0.95795179)) }
    
    add EfficiencyFormula {211} {-11} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.050) * (0.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.050) * (
							(-3.5 <= eta && eta < 1.0) * (1e-4) +
							(1.0 <= eta && eta <= 3.5) * (0.016947427)) }
    
    add EfficiencyFormula {-11} {211} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.050) * (1.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.050) * (
							(-3.5 <= eta && eta < 1.0) * (0.00) +
							(1.0 <= eta && eta <= 3.5) * (0.016947427)) }
    
    add EfficiencyFormula {211} {211} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.050) * (1.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.050) * (
							(-3.5 <= eta && eta < 1.0) * (0.99990000) +
							(1.0 <= eta && eta <= 3.5) * (0.95795179)) }

    # pi/K/proton identification, assuming 3 sigma separation between all species! (all pair-wise separations are 3 sigma)
    # Kaons must have p > 0.135 GeV/c
    # Pions must have p > 0.100 GeV/c

    add EfficiencyFormula {321} {321} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.135) * (0.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.135) * (
							(eta < -1.0 && pt * cosh(eta) <= 7) * (0.95795179) +
							(-1.0 <= eta && eta < 0.5 && pt * cosh(eta) <= 10) * (0.95795179) +
							(0.5 <= eta && eta < 1.0 && pt * cosh(eta) <= 15) * (0.95795179) +
							(1.0 <= eta && eta < 1.5 && pt * cosh(eta) <= 30) * (0.95795179) +
							(1.5 <= eta && eta < 2.5 && pt * cosh(eta) <= 50) * (0.95795179) +
							(2.5 <= eta && eta <= 3.5 && pt * cosh(eta) <= 45) * (0.95795179)) }

    add EfficiencyFormula {321} {-211} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.135) * (1.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.135) * (
							(eta < -1.0 && pt * cosh(eta) <= 7) * (0.016947427) +
							(-1.0 <= eta && eta < 0.5 && pt * cosh(eta) <= 10) * (0.016947427) +
							(0.5 <= eta && eta < 1.0 && pt * cosh(eta) <= 15) * (0.016947427) +
							(1.0 <= eta && eta < 1.5 && pt * cosh(eta) <= 30) * (0.016947427) +
							(1.5 <= eta && eta < 2.5 && pt * cosh(eta) <= 50) * (0.016947427) +
							(2.5 <= eta && eta <= 3.5 && pt * cosh(eta) <= 45) * (0.016947427)) }

    add EfficiencyFormula {321} {2212} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.135) * (1.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.135) * (
							(eta < -1.0 && pt * cosh(eta) <= 7) * (0.016947427) +
							(-1.0 <= eta && eta < 0.5 && pt * cosh(eta) <= 10) * (0.016947427) +
							(0.5 <= eta && eta < 1.0 && pt * cosh(eta) <= 15) * (0.016947427) +
							(1.0 <= eta && eta < 1.5 && pt * cosh(eta) <= 30) * (0.016947427) +
							(1.5 <= eta && eta < 2.5 && pt * cosh(eta) <= 50) * (0.016947427) +
							(2.5 <= eta && eta <= 3.5 && pt * cosh(eta) <= 45) * (0.016947427)) }

    
    add EfficiencyFormula {-211} {321} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.100) * (0.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.100) * (
							(eta < -1.0 && pt * cosh(eta) <= 7) * (0.016947427) +
							(-1.0 <= eta && eta < 0.5 && pt * cosh(eta) <= 10) * (0.016947427) +
							(0.5 <= eta && eta < 1.0 && pt * cosh(eta) <= 15) * (0.016947427) +
							(1.0 <= eta && eta < 1.5 && pt * cosh(eta) <= 30) * (0.016947427) +
							(1.5 <= eta && eta < 2.5 && pt * cosh(eta) <= 50) * (0.016947427) +
							(2.5 <= eta && eta <= 3.5 && pt * cosh(eta) <= 45) * (0.016947427)) }

    add EfficiencyFormula {-211} {2212} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.100) * (0.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.100) * (
							(eta < -1.0 && pt * cosh(eta) <= 7) * (0.016947427) +
							(-1.0 <= eta && eta < 0.5 && pt * cosh(eta) <= 10) * (0.016947427) +
							(0.5 <= eta && eta < 1.0 && pt * cosh(eta) <= 15) * (0.016947427) +
							(1.0 <= eta && eta < 1.5 && pt * cosh(eta) <= 30) * (0.016947427) +
							(1.5 <= eta && eta < 2.5 && pt * cosh(eta) <= 50) * (0.016947427) +
							(2.5 <= eta && eta <= 3.5 && pt * cosh(eta) <= 45) * (0.016947427)) }

    add EfficiencyFormula {211} {211} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.100) * (1.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.100) * (
							(eta < -1.0 && pt * cosh(eta) <= 7) * (0.95795179) +
							(-1.0 <= eta && eta < 0.5 && pt * cosh(eta) <= 10) * (0.95795179) +
							(0.5 <= eta && eta < 1.0 && pt * cosh(eta) <= 15) * (0.95795179) +
							(1.0 <= eta && eta < 1.5 && pt * cosh(eta) <= 30) * (0.95795179) +
							(1.5 <= eta && eta < 2.5 && pt * cosh(eta) <= 50) * (0.95795179) +
							(2.5 <= eta && eta <= 3.5 && pt * cosh(eta) <= 45) * (0.95795179)) }
    

    add EfficiencyFormula {2212} {2212} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.100) * (1.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.100) * (
							(eta < -1.0 && pt * cosh(eta) <= 7) * (0.95795179) +
							(-1.0 <= eta && eta < 0.5 && pt * cosh(eta) <= 10) * (0.95795179) +
							(0.5 <= eta && eta < 1.0 && pt * cosh(eta) <= 15) * (0.95795179) +
							(1.0 <= eta && eta < 1.5 && pt * cosh(eta) <= 30) * (0.95795179) +
							(1.5 <= eta && eta < 2.5 && pt * cosh(eta) <= 50) * (0.95795179) +
							(2.5 <= eta && eta <= 3.5 && pt * cosh(eta) <= 45) * (0.95795179)) }


    add EfficiencyFormula {2212} {321} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.100) * (0.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.100) * (
							(eta < -1.0 && pt * cosh(eta) <= 7) * (0.016947427) +
							(-1.0 <= eta && eta < 0.5 && pt * cosh(eta) <= 10) * (0.016947427) +
							(0.5 <= eta && eta < 1.0 && pt * cosh(eta) <= 15) * (0.016947427) +
							(1.0 <= eta && eta < 1.5 && pt * cosh(eta) <= 30) * (0.016947427) +
							(1.5 <= eta && eta < 2.5 && pt * cosh(eta) <= 50) * (0.016947427) +
							(2.5 <= eta && eta <= 3.5 && pt * cosh(eta) <= 45) * (0.016947427)) }


    add EfficiencyFormula {2212} {-221} { (abs(eta) > 3.5 || pt * cosh(eta) < 0.100) * (0.00) +
	(abs(eta) <= 3.5 && pt * cosh(eta) >= 0.100) * (
							(eta < -1.0 && pt * cosh(eta) <= 7) * (0.016947427) +
							(-1.0 <= eta && eta < 0.5 && pt * cosh(eta) <= 10) * (0.016947427) +
							(0.5 <= eta && eta < 1.0 && pt * cosh(eta) <= 15) * (0.016947427) +
							(1.0 <= eta && eta < 1.5 && pt * cosh(eta) <= 30) * (0.016947427) +
							(1.5 <= eta && eta < 2.5 && pt * cosh(eta) <= 50) * (0.016947427) +
							(2.5 <= eta && eta <= 3.5 && pt * cosh(eta) <= 45) * (0.016947427)) }

    # Everything else with no PID system coverage is 100% identified as pion
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

  add Branch PIDSystems/tracks PIDSystemsTrack Track

  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenMissingET/momentum GenMissingET MissingET
 
  add Branch UniqueObjectFinder/jets Jet Jet
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/photons Photon Photon
  
  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}


