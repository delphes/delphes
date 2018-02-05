#######################################
# CLICdet model
# based on CLICdp-Note-2017-001
# Ulrike Schnoor ulrike.schnoor@cern.ch
# 
# Jet finding with Valencia algorithm:
# use exclusive clustering with njets
# according to final state
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
    
    ECal
    HCal

    Calorimeter
    EFlowMerger

    
    PhotonEfficiency
    PhotonIsolation

    ElectronFilter
    ElectronEfficiency
    ElectronIsolation

    ChargedHadronFilter

    MuonEfficiency
    MuonIsolation

    EFlowFilter
    
    NeutrinoFilter
    GenJetFinder
    FastJetFinderKt
    FastJetFinderVLC_R05_N2
    FastJetFinderVLC_R05_N3
    FastJetFinderVLC_R05_N4
    FastJetFinderVLC_R05_N5
    FastJetFinderVLC_R05_N6
    FastJetFinderVLC_R07_N2
    FastJetFinderVLC_R07_N3
    FastJetFinderVLC_R07_N4
    FastJetFinderVLC_R07_N5
    FastJetFinderVLC_R07_N6
    FastJetFinderVLC_R10_N2
    FastJetFinderVLC_R10_N3
    FastJetFinderVLC_R10_N4
    FastJetFinderVLC_R10_N5
    FastJetFinderVLC_R10_N6
    FastJetFinderVLC_R12_N2
    FastJetFinderVLC_R12_N3
    FastJetFinderVLC_R12_N4
    FastJetFinderVLC_R12_N5
    FastJetFinderVLC_R12_N6
    FastJetFinderVLC_R15_N2
    FastJetFinderVLC_R15_N3
    FastJetFinderVLC_R15_N4
    FastJetFinderVLC_R15_N5
    FastJetFinderVLC_R15_N6



    MissingET
    GenMissingET

    

    JetFlavorAssociation_R05N2
	JetFlavorAssociation_R05N3
	JetFlavorAssociation_R05N4
	JetFlavorAssociation_R05N5
	JetFlavorAssociation_R05N6
	
	JetFlavorAssociation_R07N2
	JetFlavorAssociation_R07N3
	JetFlavorAssociation_R07N4
	JetFlavorAssociation_R07N5
	JetFlavorAssociation_R07N6

	JetFlavorAssociation_R10N2
	JetFlavorAssociation_R10N3
	JetFlavorAssociation_R10N4
	JetFlavorAssociation_R10N5
	JetFlavorAssociation_R10N6

	JetFlavorAssociation_R12N2
	JetFlavorAssociation_R12N3
	JetFlavorAssociation_R12N4
	JetFlavorAssociation_R12N5
	JetFlavorAssociation_R12N6

	JetFlavorAssociation_R15N2
	JetFlavorAssociation_R15N3
	JetFlavorAssociation_R15N4
	JetFlavorAssociation_R15N5
	JetFlavorAssociation_R15N6

	BTaggingWP50_R05N2
	BTaggingWP70_R05N2
	BTaggingWP90_R05N2
	BTaggingWP50_R05N3
	BTaggingWP70_R05N3
	BTaggingWP90_R05N3
	BTaggingWP50_R05N4
	BTaggingWP70_R05N4
	BTaggingWP90_R05N4
	BTaggingWP50_R05N5
	BTaggingWP70_R05N5
	BTaggingWP90_R05N5
	BTaggingWP50_R05N6
	BTaggingWP70_R05N6
	BTaggingWP90_R05N6
	BTaggingWP50_R07N2
	BTaggingWP70_R07N2
	BTaggingWP90_R07N2
	BTaggingWP50_R07N3
	BTaggingWP70_R07N3
	BTaggingWP90_R07N3
	BTaggingWP50_R07N4
	BTaggingWP70_R07N4
	BTaggingWP90_R07N4
	BTaggingWP50_R07N5
	BTaggingWP70_R07N5
	BTaggingWP90_R07N5
	BTaggingWP50_R07N6
	BTaggingWP70_R07N6
	BTaggingWP90_R07N6
	BTaggingWP50_R10N2
	BTaggingWP70_R10N2
	BTaggingWP90_R10N2
	BTaggingWP50_R10N3
	BTaggingWP70_R10N3
	BTaggingWP90_R10N3
	BTaggingWP50_R10N4
	BTaggingWP70_R10N4
	BTaggingWP90_R10N4
	BTaggingWP50_R10N5
	BTaggingWP70_R10N5
	BTaggingWP90_R10N5
	BTaggingWP50_R10N6
	BTaggingWP70_R10N6
	BTaggingWP90_R10N6
	BTaggingWP50_R12N2
	BTaggingWP70_R12N2
	BTaggingWP90_R12N2
	BTaggingWP50_R12N3
	BTaggingWP70_R12N3
	BTaggingWP90_R12N3
	BTaggingWP50_R12N4
	BTaggingWP70_R12N4
	BTaggingWP90_R12N4
	BTaggingWP50_R12N5
	BTaggingWP70_R12N5
	BTaggingWP90_R12N5
	BTaggingWP50_R12N6
	BTaggingWP70_R12N6
	BTaggingWP90_R12N6
	BTaggingWP50_R15N2
	BTaggingWP70_R15N2
	BTaggingWP90_R15N2
	BTaggingWP50_R15N3
	BTaggingWP70_R15N3
	BTaggingWP90_R15N3
	BTaggingWP50_R15N4
	BTaggingWP70_R15N4
	BTaggingWP90_R15N4
	BTaggingWP50_R15N5
	BTaggingWP70_R15N5
	BTaggingWP90_R15N5
	BTaggingWP50_R15N6
	BTaggingWP70_R15N6
	BTaggingWP90_R15N6

    
    TauTagging_R05N2

    ScalarHT


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

    # radius of the magnetic field coverage in the calorimeter, in m
    set Radius 1.5
    # half-length of the magnetic field coverage in the calorimeter, in m
    set HalfLength 2.31

    # magnetic field, in T
    set Bz 4.0
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
    set InputArray ParticlePropagator/chargedHadrons
    set OutputArray chargedHadrons
    # Current full simulation with CLICdet provides for pions:

    set EfficiencyFormula {
	(abs(eta) > 2.54) * (0.000) +
	(energy >= 100) * (abs(eta) < 2.54)  * (1.000) +
	(energy < 100 && energy >= 10) * (abs(eta) <=2.54 && abs(eta) > 2.34)  * (0.994) +
	(energy < 100 && energy >= 10) * (abs(eta) <= 2.34) * (1.000) +
	(energy < 10  && energy >= 1) * (abs(eta) <= 2.54 && abs(eta) > 0.55 ) * (0.000) +
	(energy < 10  && energy >= 1) * (abs(eta) <= 0.55 ) * (1.000) 
    }
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
    set InputArray ParticlePropagator/electrons
    set OutputArray electrons


    # Current full simulation with CLICdet provides for electrons:
    set EfficiencyFormula {
	(pt <= 1)                                                               * (0.000) +
	(abs(eta) > 2.54)                                                       * (0.000) +
	(energy > 100)                 * (abs(eta) <= 2.54 && abs(eta) > 2.44 ) * (0.993) +
	(energy > 100)                 * (abs(eta) <= 2.44 && abs(eta) > 2.34 ) * (0.997) +
	(energy > 100)                 * (abs(eta) <= 2.34  )                   * (1.000) +
	(energy <= 100 && energy > 10) * (abs(eta) <= 2.54 && abs(eta) > 2.17 ) * (0.998) +
	(energy <= 100 && energy > 10) * (abs(eta) <= 2.17)                     * (1.000) +
	(energy <= 10 && energy > 1)   * (abs(eta) > 2.34 )                     * (0.000) +
	(energy <= 10 && energy > 1)   * (abs(eta) <= 2.34 && abs(eta) > 0.76 ) * (0.997) +
	(energy <= 10 && energy > 1)   * (abs(eta) <= 0.76)                     * (0.999)
    }
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
    set InputArray ParticlePropagator/muons
    set OutputArray muons

    # Current full simulation with CLICdet provides for muons:
    set EfficiencyFormula {
	(pt < 1) * (0.000)+
	(abs(eta) > 2.54) * (0.0000) +
	(abs(eta) <= 2.54 && abs(eta) > 2.44 ) * (energy >= 100)                * (0.994) +
	(abs(eta) <= 2.54 && abs(eta) > 2.44 ) * (energy >= 10 && energy < 100) * (0.997) +
	(abs(eta) <= 2.54 && abs(eta) > 2.44 ) * (energy >= 1  && energy < 10)  * (0.996) +
	(abs(eta) <= 2.44 )                    * (energy >= 10)                 * (1.000) +
	(abs(eta) <= 2.44 && abs(eta) > 2.25 ) * (energy >= 1 && energy < 10)   * (0.999) +
	(abs(eta) <= 2.25 )                    * (energy >= 1)                  * (1.000) 
    }
}

########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
    set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
    set OutputArray chargedHadrons


    # Resolution given in dpT/pT.
    # CLICdet internal studies
    set ResolutionFormula {
	(abs(eta) < 2.66 && abs(eta) >= 2.03 ) * sqrt( 8.56036e-05^2 * pt^2 +0.0148987^2    ) +
	(abs(eta) < 2.03 && abs(eta) >= 1.01 ) * sqrt( 1.12382e-05^2 * pt^2 +0.00391722^2   ) +
	(abs(eta) < 1.01 && abs(eta) >= 0.55 ) * sqrt( 1.16768e-05^2 * pt^2 +0.00255204^2    ) +
	(abs(eta) < 0.55 && abs(eta) >= 0.18 ) * sqrt( 1.28327e-05^2 * pt^2 +0.00220587^2   ) +
	(abs(eta) < 0.18)                      * sqrt( 1.32845e-05^2 * pt^2 +0.00209325^2   )
	
    }
}

###################################
# Momentum resolution for electrons
###################################

module MomentumSmearing ElectronMomentumSmearing {
    set InputArray ElectronTrackingEfficiency/electrons
    set OutputArray electrons

    # Resolution given in dpT/pT.
    # CLICdet internal studies
    set ResolutionFormula {
	(abs(eta) < 2.66 && abs(eta) >= 2.03 ) * sqrt( 8.62283e-05^2 * pt^2  + 0.0177556^2   ) +
	(abs(eta) < 2.03 && abs(eta) >= 1.01 ) * sqrt( 1.0915e-05 ^2 * pt^2  + 0.00663766^2  ) +
	(abs(eta) < 1.01 && abs(eta) >= 0.55 ) * sqrt( 1.15518e-05^2 * pt^2  + 0.00398644^2  ) +
	(abs(eta) < 0.55 && abs(eta) >= 0.18 ) * sqrt( 1.3307e-05 ^2 * pt^2  + 0.00317807^2  ) +
	(abs(eta) < 0.18)                      * sqrt( 1.40722e-05^2 * pt^2  + 0.00292138^2  )

    }
}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
    set InputArray MuonTrackingEfficiency/muons
    set OutputArray muons

    # Resolution given in dpT/pT.

    # CLICdet internal studies
    set ResolutionFormula {

	(abs(eta) < 2.66 && abs(eta) >= 2.03 ) * sqrt(4.57439e-05^2 * pt^2*   + 0.0149328^2	   ) +
	(abs(eta) < 2.03 && abs(eta) >= 1.01 ) * sqrt(9.81626e-06^2 * pt^2*   + 0.00379895^2  ) +
	(abs(eta) < 1.01 && abs(eta) >= 0.55 ) * sqrt(1.1959e-05^2 * pt^2*   +  0.00242417^2 ) +
	(abs(eta) < 0.55 && abs(eta) >= 0.18 ) * sqrt(1.20149e-05^2 * pt^2  + 0.00219291^2  ) +
	(abs(eta) < 0.18)                      * sqrt(1.29686e-05^2 * pt^2  + 0.0020392^2      ) 

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

    #ECAL barrel: dphi = 0.2 degree, deta=0.003 towers up to |eta| <=1.2
    #ECAL endcaps: dphi = 0.8 degree, deta=0.02 towers up to |eta| <=2.5
    #ECAL plug: dphi = 1 degree, deta = 0.02 up to |eta| <=3
    #ECAL cell sizes always 5x5 mm^2

    #barrel:
    #dphi = 0.2 degree towers up to eta <=1.2
    set PhiBins {}
    for {set i -900} {$i <= 900} {incr i} {
	add PhiBins [expr {$i * $pi/900.0 }]
    }
    # 0.003 unit (5x5 mm^2) in eta up to eta <=1.2
    for {set i -400} {$i <=400} {incr i} {
	set eta [expr {$i * 0.003}]
	add EtaPhiBins $eta $PhiBins
    }

    #endcaps:
    #dphi = 0.8 degree towers for 1.2 < eta <=2.5
    set PhiBins {}
    for {set i -225} {$i <= 225} {incr i} {
	add PhiBins [expr {$i * $pi/225.}]
    }
    #deta=0.02 units for 1.2 < |eta| <=2.5
    #first, from -2.5 to -1.2, there will be (1.3/0.02=)65 segments
    for {set i 1} {$i <=66} {incr i} {
	set eta [expr {-2.52 + $i * 0.02}]
	add EtaPhiBins $eta $PhiBins
    }
    #same for 1.2 to 2.5
    for  {set i 1} {$i <=66} {incr i} {
	set eta [expr {1.18 + $i*0.02}]
	add EtaPhiBins $eta $PhiBins
    }
    
    #plug: 
    #dphi = 1 degree for 2.5 < eta <=3
    set PhiBins {}
    for {set i -180} {$i <= 180} {incr i} {
	add PhiBins [expr {$i * $pi/180.}]
    }
    # deta = 0.02 for 2.5 < |eta| <=3
    # from -3 to -2.5, there will be 25 segments
    for {set i 1} {$i <= 26} {incr i} {
	set eta [expr {-3.02 + $i * 0.02}]
	add EtaPhiBins $eta $PhiBins
    }
    #same for 2.5 to 3
    for {set i 1} {$i <= 26} {incr i} {
	set eta [expr {2.48 + $i*0.02}]
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
    add EnergyFraction {310} {0.3}
    add EnergyFraction {3122} {0.3}

    # set ECalResolutionFormula {resolution formula as a function of eta and energy}
    set ResolutionFormula { (abs(eta) <= 0.78)                   * sqrt(energy^2*0.01^2 + energy*0.156^2)+
	(abs(eta) > 0.78 && abs(eta) <=0.83 ) * sqrt( energy^0.01^2 + energy*0.175^2  ) +
	(abs(eta) <= 3 && abs(eta) > 0.83) * sqrt( energy^2*0.01^2 + energy*0.151^2  )}
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

    
    #HCAL barrel: dphi = 1 degree, deta= 0.02 towers up to |eta| <=0.8
    #HCAL ring: dphi = 1 degree, deta= 0.02 towers up to |eta| <=0.9
    #HCAL endcaps: dphi = 6 degree, deta = 0.1  up to |eta| <=3.5
    #HCAL cell sizes always 30x30 mm^2

    #barrel and ring:
    #dphi = 1 degree up to |eta| <=0.9
    set PhiBins {}
    for {set i -180} {$i <=180} {incr i} {
	add PhiBins [expr {$i * $pi/180.0}]
    }
    #deta= 0.02 towers up to |eta| <=0.9
    for {set i -45} {$i <=45} {incr i} {
	set eta [expr {$i * 0.02}]
	add EtaPhiBins $eta $PhiBins
    }

    #endcaps:
    # dphi = 6 degree
    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    # deta =0.1 for 0.9 < |eta| <=3.5
    #for -3.5 to -0.9, 26 segments
    for {set i 1} {$i <=27} {incr i} {
	set eta [expr {-3.6 + $i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }
    #same for 0.9 to 3.5
    for {set i 1} {$i <=27} {incr i} {
	set eta [expr {0.8 + $i * 0.1 }]
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

    # set HCalResolutionFormula {resolution formula as a function of eta and energy}
    #CLICdet internal studies
    set ResolutionFormula {
		(abs(eta)<= 0.78) * sqrt(1.48^2  + energy*0.251^2  + energy^2*0.0537^2) +
		(abs(eta)<=1.099 && abs(eta) > 0.78) * sqrt(  1.25^2 + energy*0.311^2 + energy^2*0.0515^2 ) +
		(abs(eta)<=3 && abs(eta)> 1.099) * sqrt( 1.09^2 + energy*0.321^2 + energy^2*0.0517^2  )
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


###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
    set InputArray ECal/eflowPhotons
    set OutputArray photons

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    # efficiency formula for photons
    # current full simulation of CLICdet yields:
    set EfficiencyFormula {
	(energy < 2.0 ) * (0.000) +
	(energy >= 2.0) * (abs(eta) < 0.7)*(0.94) +
	(energy >= 2.0) * (abs(eta) >=0.7 && abs(eta) <=3.0) * (0.9)	}

}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
    set CandidateInputArray PhotonEfficiency/photons
    set IsolationInputArray EFlowMerger/eflow

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
    set EfficiencyFormula {
	(energy < 5.0 ) * (0.000) +
	(energy >= 5.0) * (abs( eta)<0.7 )*(  0.89)+
	(energy >= 5.0) * (abs(eta) > 0.7 && abs(eta) < 3)*( 0.8)}

   
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
    set CandidateInputArray ElectronEfficiency/electrons
    set IsolationInputArray EFlowMerger/eflow

    set OutputArray electrons

    set DeltaRMax 0.5

    set PTMin 0.5

    set PTRatioMax 0.12
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
    set InputArray MuonMomentumSmearing/muons
    set OutputArray muons

    # set EfficiencyFormula {efficiency as a function of eta and pt}

    # efficiency formula for muons
    # current full simulation of CLICdet yields:
    set EfficiencyFormula {
	(energy < 2.0 ) * (0.000) +
	(energy  < 50 && energy >=2.0 )*(0.98) + (energy>=50) *(0.999)}
}

################
# Muon isolation
################

module Isolation MuonIsolation {
    set CandidateInputArray MuonEfficiency/muons
    set IsolationInputArray EFlowMerger/eflow

    set OutputArray muons

    set DeltaRMax 0.5

    set PTMin 0.5

    set PTRatioMax 0.25
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
######################
# EFlowFilter (UniqueObjectFinder)
######################
module UniqueObjectFinder EFlowFilter {
    add InputArray PhotonIsolation/photons photons
    add InputArray ElectronIsolation/electrons electrons
    add InputArray MuonIsolation/muons muons
    add InputArray EFlowMerger/eflow eflow
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

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt, 7 anti-kt with winner-take-all axis (for N-subjettiness), 8 N-jettiness, 9 Valencia
    set JetAlgorithm 9
    set ParameterR 0.5

    set JetPTMin 20.0
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

module FastJetFinder FastJetFinderKt {
    #  set InputArray Calorimeter/towers
    set InputArray EFlowMerger/eflow

    set OutputArray KTjets

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt, 7 anti-kt with winner-take-all axis (for N-subjettiness), 8 N-jettiness, 9 Valencia
    set JetAlgorithm 4
    set ParameterR 0.5

    set JetPTMin 20.0
}

################
# Jet finder VLC
################


source CLICdet/CLICdet_JetReco.tcl


###################
## Jet Energy Scale
###################
#
#module EnergyScale JetEnergyScale {
#    set InputArray FastJetFinderVLC_R10_N4/VLCjetsR10N4
#    set OutputArray jets
#
#    # scale formula for jets
#    set ScaleFormula {1.00}
#}


########################
# Jet Flavor Association
########################

source  CLICdet/CLICdet_JetFlavorAssociation.tcl

###########
# b-tagging
###########
# based on CLICdp-Note-2014-002    

source  CLICdet/CLICdet_BTagging.tcl


#############
# tau-tagging
#############


module TauTagging TauTagging_R05N2 {
    set ParticleInputArray Delphes/allParticles
    set PartonInputArray Delphes/partons
    set JetInputArray FastJetFinderVLC_R05_N2/VLCjetsR05N2

    set DeltaR 0.5

    set TauPTMin 1.0

    set TauEtaMax 4.0

    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
    # not from CLICdet
    # default efficiency formula (misidentification rate)
    add EfficiencyFormula {0} {0.001}
    # efficiency formula for tau-jets
    add EfficiencyFormula {15} {0.4}
}


##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
    # add Branch InputArray BranchName BranchClass
    add Branch Delphes/allParticles Particle GenParticle
    
    add Branch GenJetFinder/jets GenJet Jet

    add Branch FastJetFinderKt/KTjets KTjet Jet
    add Branch FastJetFinderVLC_R05_N2/VLCjetsR05N2 VLCjetR05N2 Jet
    add Branch FastJetFinderVLC_R05_N3/VLCjetsR05N3 VLCjetR05N3 Jet
    add Branch FastJetFinderVLC_R05_N4/VLCjetsR05N4 VLCjetR05N4 Jet
    add Branch FastJetFinderVLC_R05_N5/VLCjetsR05N5 VLCjetR05N5 Jet
    add Branch FastJetFinderVLC_R05_N6/VLCjetsR05N6 VLCjetR05N6 Jet

    add Branch FastJetFinderVLC_R07_N2/VLCjetsR07N2 VLCjetR07N2 Jet
    add Branch FastJetFinderVLC_R07_N3/VLCjetsR07N3 VLCjetR07N3 Jet
    add Branch FastJetFinderVLC_R07_N4/VLCjetsR07N4 VLCjetR07N4 Jet
    add Branch FastJetFinderVLC_R07_N5/VLCjetsR07N5 VLCjetR07N5 Jet
    add Branch FastJetFinderVLC_R07_N6/VLCjetsR07N6 VLCjetR07N6 Jet

    add Branch FastJetFinderVLC_R10_N2/VLCjetsR10N2 VLCjetR10N2 Jet
    add Branch FastJetFinderVLC_R10_N3/VLCjetsR10N3 VLCjetR10N3 Jet
    add Branch FastJetFinderVLC_R10_N4/VLCjetsR10N4 VLCjetR10N4 Jet
    add Branch FastJetFinderVLC_R10_N5/VLCjetsR10N5 VLCjetR10N5 Jet
    add Branch FastJetFinderVLC_R10_N6/VLCjetsR10N6 VLCjetR10N6 Jet

    add Branch FastJetFinderVLC_R12_N2/VLCjetsR12N2 VLCjetR12N2 Jet
    add Branch FastJetFinderVLC_R12_N3/VLCjetsR12N3 VLCjetR12N3 Jet
    add Branch FastJetFinderVLC_R12_N4/VLCjetsR12N4 VLCjetR12N4 Jet
    add Branch FastJetFinderVLC_R12_N5/VLCjetsR12N5 VLCjetR12N5 Jet
    add Branch FastJetFinderVLC_R12_N6/VLCjetsR12N6 VLCjetR12N6 Jet

    add Branch FastJetFinderVLC_R15_N2/VLCjetsR15N2 VLCjetR15N2 Jet
    add Branch FastJetFinderVLC_R15_N3/VLCjetsR15N3 VLCjetR15N3 Jet
    add Branch FastJetFinderVLC_R15_N4/VLCjetsR15N4 VLCjetR15N4 Jet
    add Branch FastJetFinderVLC_R15_N5/VLCjetsR15N5 VLCjetR15N5 Jet
    add Branch FastJetFinderVLC_R15_N6/VLCjetsR15N6 VLCjetR15N6 Jet

    add Branch GenMissingET/momentum GenMissingET MissingET

    add Branch TrackMerger/tracks Track Track
    add Branch Calorimeter/towers Tower Tower

    add Branch HCal/eflowTracks EFlowTrack Track
    add Branch ECal/eflowPhotons EFlowPhoton Tower
    add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower
    
    add Branch EFlowFilter/photons Photon Photon
    add Branch EFlowFilter/electrons Electron Electron
    add Branch EFlowFilter/muons Muon Muon
    
    add Branch MissingET/momentum MissingET MissingET
    add Branch ScalarHT/energy ScalarHT ScalarHT
}

