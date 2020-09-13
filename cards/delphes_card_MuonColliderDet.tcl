#######################################
# Muon Collider Detector TARGET model
#
# Michele Selvaggi michele.selvaggi@cern.ch
# Ulrike Schnoor ulrike.schnoor@cern.ch
#
#
# !!! DISCLAIMER !!!
#
# The parameterisation of the Muon Collider
# has to be intended as a target performance.
# This has not been validated by full simulation.
# Hybrid between FCC-hh and CLIC performance.
#
#
#######################################

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
    ParticlePropagator
    TrackMergerProp

    DenseProp
    DenseMergeTracks
    DenseTrackFilter

    ChargedHadronTrackingEfficiency
    ElectronTrackingEfficiency
    MuonTrackingEfficiency
    ForwardMuonEfficiency

    ChargedHadronMomentumSmearing
    ElectronMomentumSmearing
    MuonMomentumSmearing
    ForwardMuonMomentumSmearing

    TrackMerger

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

    MuonEfficiency
    MuonIsolation

    UniqueObjectFinder

    NeutrinoFilter
    GenJetFinder


    FastJetFinderKt
    FastJetFinderVLC_R02_N2
    FastJetFinderVLC_R02_N3
    FastJetFinderVLC_R02_N4
    FastJetFinderVLC_R02_N5
    FastJetFinderVLC_R02_N6
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

    FastJetFinderVLC_R05_inclusive
    FastJetFinderVLC_R02_inclusive
    FastJetFinderVLC_R07_inclusive
    FastJetFinderVLC_R10_inclusive
    FastJetFinderVLC_R12_inclusive
    FastJetFinderVLC_R15_inclusive

    MissingET
    GenMissingET

    JetMomentumSmearing_VLCR02N2
    JetMomentumSmearing_VLCR02N3
    JetMomentumSmearing_VLCR02N4
    JetMomentumSmearing_VLCR02N5
    JetMomentumSmearing_VLCR02N6
    JetMomentumSmearing_VLCR02_inclusive
    JetMomentumSmearing_VLCR05N2
    JetMomentumSmearing_VLCR05N3
    JetMomentumSmearing_VLCR05N4
    JetMomentumSmearing_VLCR05N5
    JetMomentumSmearing_VLCR05N6
    JetMomentumSmearing_VLCR05_inclusive
    JetMomentumSmearing_VLCR07N2
    JetMomentumSmearing_VLCR07N3
    JetMomentumSmearing_VLCR07N4
    JetMomentumSmearing_VLCR07N5
    JetMomentumSmearing_VLCR07N6
    JetMomentumSmearing_VLCR07_inclusive
    JetMomentumSmearing_VLCR10N2
    JetMomentumSmearing_VLCR10N3
    JetMomentumSmearing_VLCR10N4
    JetMomentumSmearing_VLCR10N5
    JetMomentumSmearing_VLCR10N6
    JetMomentumSmearing_VLCR10_inclusive
    JetMomentumSmearing_VLCR12N2
    JetMomentumSmearing_VLCR12N3
    JetMomentumSmearing_VLCR12N4
    JetMomentumSmearing_VLCR12N5
    JetMomentumSmearing_VLCR12N6
    JetMomentumSmearing_VLCR12_inclusive
    JetMomentumSmearing_VLCR15N2
    JetMomentumSmearing_VLCR15N3
    JetMomentumSmearing_VLCR15N4
    JetMomentumSmearing_VLCR15N5
    JetMomentumSmearing_VLCR15N6
    JetMomentumSmearing_VLCR15_inclusive


    JetFlavorAssociation_R02N2
    JetFlavorAssociation_R02N3
    JetFlavorAssociation_R02N4
    JetFlavorAssociation_R02N5
    JetFlavorAssociation_R02N6
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

    JetFlavorAssociation_R02_inclusive
    JetFlavorAssociation_R05_inclusive
    JetFlavorAssociation_R07_inclusive
    JetFlavorAssociation_R10_inclusive
    JetFlavorAssociation_R12_inclusive
    JetFlavorAssociation_R15_inclusive


    BTagging_WP50_R02N2
    BTagging_WP70_R02N2
    BTagging_WP90_R02N2
    BTagging_WP50_R02N3
    BTagging_WP70_R02N3
    BTagging_WP90_R02N3
    BTagging_WP50_R02N4
    BTagging_WP70_R02N4
    BTagging_WP90_R02N4
    BTagging_WP50_R02N5
    BTagging_WP70_R02N5
    BTagging_WP90_R02N5
    BTagging_WP50_R02N6
    BTagging_WP70_R02N6
    BTagging_WP90_R02N6
    BTagging_WP50_R05N2
    BTagging_WP70_R05N2
    BTagging_WP90_R05N2
    BTagging_WP50_R05N3
    BTagging_WP70_R05N3
    BTagging_WP90_R05N3
    BTagging_WP50_R05N4
    BTagging_WP70_R05N4
    BTagging_WP90_R05N4
    BTagging_WP50_R05N5
    BTagging_WP70_R05N5
    BTagging_WP90_R05N5
    BTagging_WP50_R05N6
    BTagging_WP70_R05N6
    BTagging_WP90_R05N6
    BTagging_WP50_R07N2
    BTagging_WP70_R07N2
    BTagging_WP90_R07N2
    BTagging_WP50_R07N3
    BTagging_WP70_R07N3
    BTagging_WP90_R07N3
    BTagging_WP50_R07N4
    BTagging_WP70_R07N4
    BTagging_WP90_R07N4
    BTagging_WP50_R07N5
    BTagging_WP70_R07N5
    BTagging_WP90_R07N5
    BTagging_WP50_R07N6
    BTagging_WP70_R07N6
    BTagging_WP90_R07N6
    BTagging_WP50_R10N2
    BTagging_WP70_R10N2
    BTagging_WP90_R10N2
    BTagging_WP50_R10N3
    BTagging_WP70_R10N3
    BTagging_WP90_R10N3
    BTagging_WP50_R10N4
    BTagging_WP70_R10N4
    BTagging_WP90_R10N4
    BTagging_WP50_R10N5
    BTagging_WP70_R10N5
    BTagging_WP90_R10N5
    BTagging_WP50_R10N6
    BTagging_WP70_R10N6
    BTagging_WP90_R10N6
    BTagging_WP50_R12N2
    BTagging_WP70_R12N2
    BTagging_WP90_R12N2
    BTagging_WP50_R12N3
    BTagging_WP70_R12N3
    BTagging_WP90_R12N3
    BTagging_WP50_R12N4
    BTagging_WP70_R12N4
    BTagging_WP90_R12N4
    BTagging_WP50_R12N5
    BTagging_WP70_R12N5
    BTagging_WP90_R12N5
    BTagging_WP50_R12N6
    BTagging_WP70_R12N6
    BTagging_WP90_R12N6
    BTagging_WP50_R15N2
    BTagging_WP70_R15N2
    BTagging_WP90_R15N2
    BTagging_WP50_R15N3
    BTagging_WP70_R15N3
    BTagging_WP90_R15N3
    BTagging_WP50_R15N4
    BTagging_WP70_R15N4
    BTagging_WP90_R15N4
    BTagging_WP50_R15N5
    BTagging_WP70_R15N5
    BTagging_WP90_R15N5
    BTagging_WP50_R15N6
    BTagging_WP70_R15N6
    BTagging_WP90_R15N6
    BTagging_WP50_R02_inclusive
    BTagging_WP70_R02_inclusive
    BTagging_WP90_R02_inclusive
    BTagging_WP50_R05_inclusive
    BTagging_WP70_R05_inclusive
    BTagging_WP90_R05_inclusive
    BTagging_WP50_R07_inclusive
    BTagging_WP70_R07_inclusive
    BTagging_WP90_R07_inclusive
    BTagging_WP50_R10_inclusive
    BTagging_WP70_R10_inclusive
    BTagging_WP90_R10_inclusive
    BTagging_WP50_R12_inclusive
    BTagging_WP70_R12_inclusive
    BTagging_WP90_R12_inclusive
    BTagging_WP50_R15_inclusive
    BTagging_WP70_R15_inclusive
    BTagging_WP90_R15_inclusive


    TauTagging_R02N2
    TauTagging_R02N3
    TauTagging_R02N4
    TauTagging_R02N5
    TauTagging_R02N6
    TauTagging_R05N2
    TauTagging_R05N3
    TauTagging_R05N4
    TauTagging_R05N5
    TauTagging_R05N6
    TauTagging_R07N2
    TauTagging_R07N3
    TauTagging_R07N4
    TauTagging_R07N5
    TauTagging_R07N6
    TauTagging_R10N2
    TauTagging_R10N3
    TauTagging_R10N4
    TauTagging_R10N5
    TauTagging_R10N6
    TauTagging_R12N2
    TauTagging_R12N3
    TauTagging_R12N4
    TauTagging_R12N5
    TauTagging_R12N6
    TauTagging_R15N2
    TauTagging_R15N3
    TauTagging_R15N4
    TauTagging_R15N5
    TauTagging_R15N6
    TauTagging_R02_inclusive
    TauTagging_R05_inclusive
    TauTagging_R07_inclusive
    TauTagging_R10_inclusive
    TauTagging_R12_inclusive
    TauTagging_R15_inclusive

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


##############
# Track merger
##############

module Merger TrackMergerProp {
# add InputArray InputArray
  add InputArray ParticlePropagator/chargedHadrons
  add InputArray ParticlePropagator/electrons
  add InputArray ParticlePropagator/muons
  set OutputArray tracks
}



####################################
# Track propagation to pseudo-pixel
####################################

module ParticlePropagator DenseProp {

  set InputArray TrackMergerProp/tracks

  # radius of the magnetic field coverage, in m
  set Radius 0.45
  set RadiusMax 1.5
  # half-length of the magnetic field coverage, in m
  set HalfLength 0.8
  set HalfLengthMax 2.31

  # magnetic field
  set Bz 4.0
}

#####################
# Dense Track merger
#####################

module Merger DenseMergeTracks {
# add InputArray InputArray
  add InputArray DenseProp/chargedHadrons
  add InputArray DenseProp/electrons
  add InputArray DenseProp/muons
  set OutputArray tracks
}


######################
#   Dense Track Filter
######################

module DenseTrackFilter DenseTrackFilter {

  set TrackInputArray DenseMergeTracks/tracks

  set TrackOutputArray tracks
  set ChargedHadronOutputArray chargedHadrons
  set ElectronOutputArray electrons
  set MuonOutputArray muons

  set EtaPhiRes 0.003
  set EtaMax 2.5

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
}




####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
    set InputArray DenseTrackFilter/chargedHadrons
    set OutputArray chargedHadrons
    # tracking efficiency formula for charged hadrons

    set EfficiencyFormula { (pt <= 0.5) * (0.00) +
    (abs(eta) <= 2.0) * (pt > 0.5 && pt <= 1) * (0.90) +
    (abs(eta) <= 2.0) * (pt > 1) * (0.95) +
    (abs(eta) > 2.0 && abs(eta) < 2.5) * (pt > 0.5 && pt <= 1) * (0.80) +
    (abs(eta) > 2.0 && abs(eta) < 2.5) * (pt > 1.0) * (0.85) +
    (abs(eta) > 2.5 ) * (0.00)
   }
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
    set InputArray DenseTrackFilter/electrons
    set OutputArray electrons

    set EfficiencyFormula { (pt <= 0.5) * (0.00) +
    (abs(eta) <= 2.0) * (pt > 0.5 && pt <= 1) * (0.90) +
    (abs(eta) <= 2.0) * (pt > 1) * (0.95) +
    (abs(eta) > 2.0 && abs(eta) < 2.5) * (pt > 0.5 && pt <= 1) * (0.80) +
    (abs(eta) > 2.0 && abs(eta) < 2.5) * (pt > 1.0) * (0.85) +
    (abs(eta) > 2.5 ) * (0.00)
   }
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
    set InputArray DenseTrackFilter/muons
    set OutputArray muons

    set EfficiencyFormula { (pt <= 0.5) * (0.00) +
    (abs(eta) <= 2.0) * (pt > 0.5 && pt <= 1) * (0.95) +
    (abs(eta) <= 2.0) * (pt > 1) * (0.99) +
    (abs(eta) > 2.0 && abs(eta) < 2.5) * (pt > 0.5 && pt <= 1) * (0.90) +
    (abs(eta) > 2.0 && abs(eta) < 2.5) * (pt > 1.0) * (0.95) +
    (abs(eta) > 2.5 ) * (0.00)
   }
}

##########################
# Forward Muon efficiency
##########################

## hypothetical forward muon spectrometer
module Efficiency ForwardMuonEfficiency {
    set InputArray ParticlePropagator/muons
    set OutputArray muons

    set EfficiencyFormula { (pt <= 0.5) * (0.00) +
    (abs(eta) > 2.5 && abs(eta) < 6.0) * (pt > 0.5 && pt <= 1) * (0.90) +
    (abs(eta) > 2.5 && abs(eta) < 6.0) * (pt > 1.0) * (0.95) +
    (abs(eta) > 6.0 ) * (0.00)
   }
}


########################################
# Momentum resolution for charged tracks
########################################

module MomentumSmearing ChargedHadronMomentumSmearing {
    set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
    set OutputArray chargedHadrons


    # Resolution given in dpT/pT (from FCC-hh)
    set ResolutionFormula {    (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 0.0000 && energy < 1.0000) * (0.00315864) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 1.0000 && energy < 2.0000) * (0.003159 + (energy-1.000000)* 0.000007) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 2.0000 && energy < 5.0000) * (0.003166 + (energy-2.000000)* 0.000011) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 5.0000 && energy < 10.0000) * (0.003198 + (energy-5.000000)* 0.000012) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 10.0000 && energy < 100.0000) * (0.003259 + (energy-10.000000)* 0.000010) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004173 + (energy-100.000000)* 0.000019) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.020916 + (energy-1000.000000)* 0.000021) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 10000.0000) * (0.205876*energy/10000.000000) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 0.0000 && energy < 1.0000) * (0.00316278) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 1.0000 && energy < 2.0000) * (0.003163 + (energy-1.000000)* 0.000006) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 2.0000 && energy < 5.0000) * (0.003169 + (energy-2.000000)* 0.000010) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 5.0000 && energy < 10.0000) * (0.003198 + (energy-5.000000)* 0.000011) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 10.0000 && energy < 100.0000) * (0.003255 + (energy-10.000000)* 0.000010) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004165 + (energy-100.000000)* 0.000019) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.020917 + (energy-1000.000000)* 0.000021) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 10000.0000) * (0.205952*energy/10000.000000) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 0.0000 && energy < 1.0000) * (0.00320482) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 1.0000 && energy < 2.0000) * (0.003205 + (energy-1.000000)* 0.000006) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 2.0000 && energy < 5.0000) * (0.003211 + (energy-2.000000)* 0.000009) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 5.0000 && energy < 10.0000) * (0.003238 + (energy-5.000000)* 0.000011) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 10.0000 && energy < 100.0000) * (0.003294 + (energy-10.000000)* 0.000010) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004176 + (energy-100.000000)* 0.000018) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.020586 + (energy-1000.000000)* 0.000020) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 10000.0000) * (0.202528*energy/10000.000000) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 0.0000 && energy < 1.0000) * (0.00325680) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 1.0000 && energy < 2.0000) * (0.003257 + (energy-1.000000)* 0.000001) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 2.0000 && energy < 5.0000) * (0.003257 + (energy-2.000000)* 0.000009) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 5.0000 && energy < 10.0000) * (0.003286 + (energy-5.000000)* 0.000011) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 10.0000 && energy < 100.0000) * (0.003342 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004195 + (energy-100.000000)* 0.000017) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.019873 + (energy-1000.000000)* 0.000019) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 10000.0000) * (0.195142*energy/10000.000000) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 0.0000 && energy < 1.0000) * (0.00354020) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 1.0000 && energy < 2.0000) * (0.003540 + (energy-1.000000)* -0.000201) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 2.0000 && energy < 5.0000) * (0.003340 + (energy-2.000000)* 0.000009) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 5.0000 && energy < 10.0000) * (0.003366 + (energy-5.000000)* 0.000011) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 10.0000 && energy < 100.0000) * (0.003422 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004265 + (energy-100.000000)* 0.000017) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.019240 + (energy-1000.000000)* 0.000019) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 10000.0000) * (0.188429*energy/10000.000000) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 0.0000 && energy < 1.0000) * (0.00362672) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 1.0000 && energy < 2.0000) * (0.003627 + (energy-1.000000)* -0.000223) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 2.0000 && energy < 5.0000) * (0.003403 + (energy-2.000000)* 0.000008) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 5.0000 && energy < 10.0000) * (0.003428 + (energy-5.000000)* 0.000010) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 10.0000 && energy < 100.0000) * (0.003479 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004256 + (energy-100.000000)* 0.000016) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.018471 + (energy-1000.000000)* 0.000018) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 10000.0000) * (0.180531*energy/10000.000000) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 0.0000 && energy < 1.0000) * (0.00386864) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 1.0000 && energy < 2.0000) * (0.003869 + (energy-1.000000)* -0.000354) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 2.0000 && energy < 5.0000) * (0.003515 + (energy-2.000000)* 0.000007) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 5.0000 && energy < 10.0000) * (0.003536 + (energy-5.000000)* 0.000009) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 10.0000 && energy < 100.0000) * (0.003583 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004349 + (energy-100.000000)* 0.000015) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.017442 + (energy-1000.000000)* 0.000017) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 10000.0000) * (0.169559*energy/10000.000000) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 0.0000 && energy < 1.0000) * (0.00418203) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 1.0000 && energy < 2.0000) * (0.004182 + (energy-1.000000)* -0.000556) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 2.0000 && energy < 5.0000) * (0.003626 + (energy-2.000000)* 0.000007) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 5.0000 && energy < 10.0000) * (0.003645 + (energy-5.000000)* 0.000008) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 10.0000 && energy < 100.0000) * (0.003687 + (energy-10.000000)* 0.000008) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004399 + (energy-100.000000)* 0.000013) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.016509 + (energy-1000.000000)* 0.000016) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 10000.0000) * (0.159676*energy/10000.000000) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 0.0000 && energy < 1.0000) * (0.00436103) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 1.0000 && energy < 2.0000) * (0.004361 + (energy-1.000000)* -0.000597) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 2.0000 && energy < 5.0000) * (0.003764 + (energy-2.000000)* 0.000006) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 5.0000 && energy < 10.0000) * (0.003781 + (energy-5.000000)* 0.000008) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 10.0000 && energy < 100.0000) * (0.003821 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004495 + (energy-100.000000)* 0.000012) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.015532 + (energy-1000.000000)* 0.000015) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 10000.0000) * (0.149090*energy/10000.000000) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 0.0000 && energy < 1.0000) * (0.00488279) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 1.0000 && energy < 2.0000) * (0.004883 + (energy-1.000000)* -0.000969) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 2.0000 && energy < 5.0000) * (0.003914 + (energy-2.000000)* 0.000006) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 5.0000 && energy < 10.0000) * (0.003930 + (energy-5.000000)* 0.000007) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 10.0000 && energy < 100.0000) * (0.003967 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004593 + (energy-100.000000)* 0.000011) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.014592 + (energy-1000.000000)* 0.000014) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 10000.0000) * (0.138764*energy/10000.000000) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 0.0000 && energy < 1.0000) * (0.00513716) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 1.0000 && energy < 2.0000) * (0.005137 + (energy-1.000000)* -0.001026) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 2.0000 && energy < 5.0000) * (0.004111 + (energy-2.000000)* 0.000005) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 5.0000 && energy < 10.0000) * (0.004125 + (energy-5.000000)* 0.000007) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 10.0000 && energy < 100.0000) * (0.004159 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004752 + (energy-100.000000)* 0.000010) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.013718 + (energy-1000.000000)* 0.000013) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 10000.0000) * (0.128750*energy/10000.000000) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 0.0000 && energy < 1.0000) * (0.00572019) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 1.0000 && energy < 2.0000) * (0.005720 + (energy-1.000000)* -0.001362) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 2.0000 && energy < 5.0000) * (0.004359 + (energy-2.000000)* 0.000005) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 5.0000 && energy < 10.0000) * (0.004372 + (energy-5.000000)* 0.000006) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 10.0000 && energy < 100.0000) * (0.004405 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005024 + (energy-100.000000)* 0.000010) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.013902 + (energy-1000.000000)* 0.000013) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 10000.0000) * (0.129437*energy/10000.000000) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 0.0000 && energy < 1.0000) * (0.00613558) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 1.0000 && energy < 2.0000) * (0.006136 + (energy-1.000000)* -0.001331) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 2.0000 && energy < 5.0000) * (0.004805 + (energy-2.000000)* 0.000005) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 5.0000 && energy < 10.0000) * (0.004818 + (energy-5.000000)* 0.000007) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 10.0000 && energy < 100.0000) * (0.004851 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005512 + (energy-100.000000)* 0.000012) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.016144 + (energy-1000.000000)* 0.000015) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 10000.0000) * (0.151739*energy/10000.000000) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 0.0000 && energy < 1.0000) * (0.00655464) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 1.0000 && energy < 2.0000) * (0.006555 + (energy-1.000000)* -0.001843) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 2.0000 && energy < 5.0000) * (0.004711 + (energy-2.000000)* 0.000004) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 5.0000 && energy < 10.0000) * (0.004724 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 10.0000 && energy < 100.0000) * (0.004748 + (energy-10.000000)* 0.000006) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005326 + (energy-100.000000)* 0.000009) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.013686 + (energy-1000.000000)* 0.000012) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 10000.0000) * (0.125361*energy/10000.000000) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 0.0000 && energy < 1.0000) * (0.00714442) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 1.0000 && energy < 2.0000) * (0.007144 + (energy-1.000000)* -0.002674) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 2.0000 && energy < 5.0000) * (0.004470 + (energy-2.000000)* -0.000070) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 5.0000 && energy < 10.0000) * (0.004259 + (energy-5.000000)* 0.000006) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 10.0000 && energy < 100.0000) * (0.004287 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005101 + (energy-100.000000)* 0.000008) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.012303 + (energy-1000.000000)* 0.000011) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 10000.0000) * (0.110091*energy/10000.000000) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 0.0000 && energy < 1.0000) * (0.00680449) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 1.0000 && energy < 2.0000) * (0.006804 + (energy-1.000000)* -0.002108) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 2.0000 && energy < 5.0000) * (0.004696 + (energy-2.000000)* -0.000116) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 5.0000 && energy < 10.0000) * (0.004348 + (energy-5.000000)* 0.000009) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 10.0000 && energy < 100.0000) * (0.004391 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005232 + (energy-100.000000)* 0.000007) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.011392 + (energy-1000.000000)* 0.000010) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 10000.0000) * (0.097988*energy/10000.000000) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 0.0000 && energy < 1.0000) * (0.00763793) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 1.0000 && energy < 2.0000) * (0.007638 + (energy-1.000000)* -0.003061) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 2.0000 && energy < 5.0000) * (0.004577 + (energy-2.000000)* -0.000168) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 5.0000 && energy < 10.0000) * (0.004074 + (energy-5.000000)* 0.000009) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 10.0000 && energy < 100.0000) * (0.004121 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004910 + (energy-100.000000)* 0.000006) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.010635 + (energy-1000.000000)* 0.000009) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 10000.0000) * (0.090704*energy/10000.000000) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 0.0000 && energy < 1.0000) * (0.00913948) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 1.0000 && energy < 2.0000) * (0.009139 + (energy-1.000000)* -0.004405) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 2.0000 && energy < 5.0000) * (0.004735 + (energy-2.000000)* -0.000231) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 5.0000 && energy < 10.0000) * (0.004043 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 10.0000 && energy < 100.0000) * (0.004066 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004739 + (energy-100.000000)* 0.000006) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.009891 + (energy-1000.000000)* 0.000008) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 10000.0000) * (0.083216*energy/10000.000000) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 0.0000 && energy < 1.0000) * (0.00956747) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 1.0000 && energy < 2.0000) * (0.009567 + (energy-1.000000)* -0.004497) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 2.0000 && energy < 5.0000) * (0.005070 + (energy-2.000000)* -0.000331) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 5.0000 && energy < 10.0000) * (0.004078 + (energy-5.000000)* 0.000004) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 10.0000 && energy < 100.0000) * (0.004100 + (energy-10.000000)* 0.000006) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004682 + (energy-100.000000)* 0.000005) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.009157 + (energy-1000.000000)* 0.000007) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 10000.0000) * (0.074702*energy/10000.000000) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 0.0000 && energy < 1.0000) * (0.00964334) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 1.0000 && energy < 2.0000) * (0.009643 + (energy-1.000000)* -0.003950) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 2.0000 && energy < 5.0000) * (0.005694 + (energy-2.000000)* -0.000479) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 5.0000 && energy < 10.0000) * (0.004256 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 10.0000 && energy < 100.0000) * (0.004281 + (energy-10.000000)* 0.000006) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004862 + (energy-100.000000)* 0.000004) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.008501 + (energy-1000.000000)* 0.000006) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 10000.0000) * (0.062525*energy/10000.000000) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 0.0000 && energy < 1.0000) * (0.01045039) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 1.0000 && energy < 2.0000) * (0.010450 + (energy-1.000000)* -0.005379) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 2.0000 && energy < 5.0000) * (0.005072 + (energy-2.000000)* -0.000321) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 5.0000 && energy < 10.0000) * (0.004109 + (energy-5.000000)* 0.000006) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 10.0000 && energy < 100.0000) * (0.004137 + (energy-10.000000)* 0.000006) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004666 + (energy-100.000000)* 0.000005) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.008951 + (energy-1000.000000)* 0.000007) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 10000.0000) * (0.073400*energy/10000.000000) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 0.0000 && energy < 1.0000) * (0.01046694) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 1.0000 && energy < 2.0000) * (0.010467 + (energy-1.000000)* -0.005023) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 2.0000 && energy < 5.0000) * (0.005444 + (energy-2.000000)* -0.000330) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 5.0000 && energy < 10.0000) * (0.004455 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 10.0000 && energy < 100.0000) * (0.004479 + (energy-10.000000)* 0.000004) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004883 + (energy-100.000000)* 0.000005) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.009382 + (energy-1000.000000)* 0.000008) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 10000.0000) * (0.078852*energy/10000.000000) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 0.0000 && energy < 1.0000) * (0.01090933) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 1.0000 && energy < 2.0000) * (0.010909 + (energy-1.000000)* -0.005299) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 2.0000 && energy < 5.0000) * (0.005610 + (energy-2.000000)* -0.000302) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 5.0000 && energy < 10.0000) * (0.004704 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 10.0000 && energy < 100.0000) * (0.004730 + (energy-10.000000)* 0.000005) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005146 + (energy-100.000000)* 0.000006) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.010323 + (energy-1000.000000)* 0.000009) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 10000.0000) * (0.088469*energy/10000.000000) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 0.0000 && energy < 1.0000) * (0.01271833) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 1.0000 && energy < 2.0000) * (0.012718 + (energy-1.000000)* -0.005764) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 2.0000 && energy < 5.0000) * (0.006954 + (energy-2.000000)* -0.000492) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 5.0000 && energy < 10.0000) * (0.005479 + (energy-5.000000)* 0.000003) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 10.0000 && energy < 100.0000) * (0.005495 + (energy-10.000000)* 0.000003) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005774 + (energy-100.000000)* 0.000006) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.011265 + (energy-1000.000000)* 0.000009) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 10000.0000) * (0.096592*energy/10000.000000) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 0.0000 && energy < 1.0000) * (0.01515272) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 1.0000 && energy < 2.0000) * (0.015153 + (energy-1.000000)* -0.007272) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 2.0000 && energy < 5.0000) * (0.007881 + (energy-2.000000)* -0.000660) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 5.0000 && energy < 10.0000) * (0.005900 + (energy-5.000000)* 0.000003) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 10.0000 && energy < 100.0000) * (0.005914 + (energy-10.000000)* 0.000003) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 100.0000 && energy < 1000.0000) * (0.006174 + (energy-100.000000)* 0.000007) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.012486 + (energy-1000.000000)* 0.000011) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 10000.0000) * (0.108659*energy/10000.000000)
    }
}

###################################
# Momentum resolution for electrons
###################################

module MomentumSmearing ElectronMomentumSmearing {
    set InputArray ElectronTrackingEfficiency/electrons
    set OutputArray electrons

    # Resolution given in dpT/pT (from FCC-hh)
    set ResolutionFormula {    (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 0.0000 && energy < 1.0000) * (0.00315864) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 1.0000 && energy < 2.0000) * (0.003159 + (energy-1.000000)* 0.000007) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 2.0000 && energy < 5.0000) * (0.003166 + (energy-2.000000)* 0.000011) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 5.0000 && energy < 10.0000) * (0.003198 + (energy-5.000000)* 0.000012) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 10.0000 && energy < 100.0000) * (0.003259 + (energy-10.000000)* 0.000010) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004173 + (energy-100.000000)* 0.000019) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.020916 + (energy-1000.000000)* 0.000021) + \
       (abs(eta) >= 0.0000 && abs(eta) < 0.1000) * (energy >= 10000.0000) * (0.205876*energy/10000.000000) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 0.0000 && energy < 1.0000) * (0.00316278) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 1.0000 && energy < 2.0000) * (0.003163 + (energy-1.000000)* 0.000006) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 2.0000 && energy < 5.0000) * (0.003169 + (energy-2.000000)* 0.000010) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 5.0000 && energy < 10.0000) * (0.003198 + (energy-5.000000)* 0.000011) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 10.0000 && energy < 100.0000) * (0.003255 + (energy-10.000000)* 0.000010) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004165 + (energy-100.000000)* 0.000019) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.020917 + (energy-1000.000000)* 0.000021) + \
       (abs(eta) >= 0.1000 && abs(eta) < 0.2000) * (energy >= 10000.0000) * (0.205952*energy/10000.000000) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 0.0000 && energy < 1.0000) * (0.00320482) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 1.0000 && energy < 2.0000) * (0.003205 + (energy-1.000000)* 0.000006) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 2.0000 && energy < 5.0000) * (0.003211 + (energy-2.000000)* 0.000009) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 5.0000 && energy < 10.0000) * (0.003238 + (energy-5.000000)* 0.000011) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 10.0000 && energy < 100.0000) * (0.003294 + (energy-10.000000)* 0.000010) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004176 + (energy-100.000000)* 0.000018) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.020586 + (energy-1000.000000)* 0.000020) + \
       (abs(eta) >= 0.2000 && abs(eta) < 0.3000) * (energy >= 10000.0000) * (0.202528*energy/10000.000000) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 0.0000 && energy < 1.0000) * (0.00325680) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 1.0000 && energy < 2.0000) * (0.003257 + (energy-1.000000)* 0.000001) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 2.0000 && energy < 5.0000) * (0.003257 + (energy-2.000000)* 0.000009) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 5.0000 && energy < 10.0000) * (0.003286 + (energy-5.000000)* 0.000011) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 10.0000 && energy < 100.0000) * (0.003342 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004195 + (energy-100.000000)* 0.000017) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.019873 + (energy-1000.000000)* 0.000019) + \
       (abs(eta) >= 0.3000 && abs(eta) < 0.4000) * (energy >= 10000.0000) * (0.195142*energy/10000.000000) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 0.0000 && energy < 1.0000) * (0.00354020) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 1.0000 && energy < 2.0000) * (0.003540 + (energy-1.000000)* -0.000201) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 2.0000 && energy < 5.0000) * (0.003340 + (energy-2.000000)* 0.000009) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 5.0000 && energy < 10.0000) * (0.003366 + (energy-5.000000)* 0.000011) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 10.0000 && energy < 100.0000) * (0.003422 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004265 + (energy-100.000000)* 0.000017) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.019240 + (energy-1000.000000)* 0.000019) + \
       (abs(eta) >= 0.4000 && abs(eta) < 0.5000) * (energy >= 10000.0000) * (0.188429*energy/10000.000000) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 0.0000 && energy < 1.0000) * (0.00362672) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 1.0000 && energy < 2.0000) * (0.003627 + (energy-1.000000)* -0.000223) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 2.0000 && energy < 5.0000) * (0.003403 + (energy-2.000000)* 0.000008) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 5.0000 && energy < 10.0000) * (0.003428 + (energy-5.000000)* 0.000010) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 10.0000 && energy < 100.0000) * (0.003479 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004256 + (energy-100.000000)* 0.000016) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.018471 + (energy-1000.000000)* 0.000018) + \
       (abs(eta) >= 0.5000 && abs(eta) < 0.6000) * (energy >= 10000.0000) * (0.180531*energy/10000.000000) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 0.0000 && energy < 1.0000) * (0.00386864) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 1.0000 && energy < 2.0000) * (0.003869 + (energy-1.000000)* -0.000354) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 2.0000 && energy < 5.0000) * (0.003515 + (energy-2.000000)* 0.000007) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 5.0000 && energy < 10.0000) * (0.003536 + (energy-5.000000)* 0.000009) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 10.0000 && energy < 100.0000) * (0.003583 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004349 + (energy-100.000000)* 0.000015) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.017442 + (energy-1000.000000)* 0.000017) + \
       (abs(eta) >= 0.6000 && abs(eta) < 0.7000) * (energy >= 10000.0000) * (0.169559*energy/10000.000000) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 0.0000 && energy < 1.0000) * (0.00418203) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 1.0000 && energy < 2.0000) * (0.004182 + (energy-1.000000)* -0.000556) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 2.0000 && energy < 5.0000) * (0.003626 + (energy-2.000000)* 0.000007) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 5.0000 && energy < 10.0000) * (0.003645 + (energy-5.000000)* 0.000008) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 10.0000 && energy < 100.0000) * (0.003687 + (energy-10.000000)* 0.000008) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004399 + (energy-100.000000)* 0.000013) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.016509 + (energy-1000.000000)* 0.000016) + \
       (abs(eta) >= 0.7000 && abs(eta) < 0.8000) * (energy >= 10000.0000) * (0.159676*energy/10000.000000) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 0.0000 && energy < 1.0000) * (0.00436103) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 1.0000 && energy < 2.0000) * (0.004361 + (energy-1.000000)* -0.000597) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 2.0000 && energy < 5.0000) * (0.003764 + (energy-2.000000)* 0.000006) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 5.0000 && energy < 10.0000) * (0.003781 + (energy-5.000000)* 0.000008) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 10.0000 && energy < 100.0000) * (0.003821 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004495 + (energy-100.000000)* 0.000012) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.015532 + (energy-1000.000000)* 0.000015) + \
       (abs(eta) >= 0.8000 && abs(eta) < 0.9000) * (energy >= 10000.0000) * (0.149090*energy/10000.000000) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 0.0000 && energy < 1.0000) * (0.00488279) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 1.0000 && energy < 2.0000) * (0.004883 + (energy-1.000000)* -0.000969) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 2.0000 && energy < 5.0000) * (0.003914 + (energy-2.000000)* 0.000006) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 5.0000 && energy < 10.0000) * (0.003930 + (energy-5.000000)* 0.000007) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 10.0000 && energy < 100.0000) * (0.003967 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004593 + (energy-100.000000)* 0.000011) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.014592 + (energy-1000.000000)* 0.000014) + \
       (abs(eta) >= 0.9000 && abs(eta) < 1.0000) * (energy >= 10000.0000) * (0.138764*energy/10000.000000) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 0.0000 && energy < 1.0000) * (0.00513716) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 1.0000 && energy < 2.0000) * (0.005137 + (energy-1.000000)* -0.001026) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 2.0000 && energy < 5.0000) * (0.004111 + (energy-2.000000)* 0.000005) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 5.0000 && energy < 10.0000) * (0.004125 + (energy-5.000000)* 0.000007) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 10.0000 && energy < 100.0000) * (0.004159 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004752 + (energy-100.000000)* 0.000010) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.013718 + (energy-1000.000000)* 0.000013) + \
       (abs(eta) >= 1.0000 && abs(eta) < 1.1000) * (energy >= 10000.0000) * (0.128750*energy/10000.000000) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 0.0000 && energy < 1.0000) * (0.00572019) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 1.0000 && energy < 2.0000) * (0.005720 + (energy-1.000000)* -0.001362) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 2.0000 && energy < 5.0000) * (0.004359 + (energy-2.000000)* 0.000005) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 5.0000 && energy < 10.0000) * (0.004372 + (energy-5.000000)* 0.000006) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 10.0000 && energy < 100.0000) * (0.004405 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005024 + (energy-100.000000)* 0.000010) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.013902 + (energy-1000.000000)* 0.000013) + \
       (abs(eta) >= 1.1000 && abs(eta) < 1.2000) * (energy >= 10000.0000) * (0.129437*energy/10000.000000) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 0.0000 && energy < 1.0000) * (0.00613558) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 1.0000 && energy < 2.0000) * (0.006136 + (energy-1.000000)* -0.001331) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 2.0000 && energy < 5.0000) * (0.004805 + (energy-2.000000)* 0.000005) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 5.0000 && energy < 10.0000) * (0.004818 + (energy-5.000000)* 0.000007) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 10.0000 && energy < 100.0000) * (0.004851 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005512 + (energy-100.000000)* 0.000012) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.016144 + (energy-1000.000000)* 0.000015) + \
       (abs(eta) >= 1.2000 && abs(eta) < 1.3000) * (energy >= 10000.0000) * (0.151739*energy/10000.000000) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 0.0000 && energy < 1.0000) * (0.00655464) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 1.0000 && energy < 2.0000) * (0.006555 + (energy-1.000000)* -0.001843) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 2.0000 && energy < 5.0000) * (0.004711 + (energy-2.000000)* 0.000004) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 5.0000 && energy < 10.0000) * (0.004724 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 10.0000 && energy < 100.0000) * (0.004748 + (energy-10.000000)* 0.000006) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005326 + (energy-100.000000)* 0.000009) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.013686 + (energy-1000.000000)* 0.000012) + \
       (abs(eta) >= 1.3000 && abs(eta) < 1.4000) * (energy >= 10000.0000) * (0.125361*energy/10000.000000) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 0.0000 && energy < 1.0000) * (0.00714442) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 1.0000 && energy < 2.0000) * (0.007144 + (energy-1.000000)* -0.002674) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 2.0000 && energy < 5.0000) * (0.004470 + (energy-2.000000)* -0.000070) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 5.0000 && energy < 10.0000) * (0.004259 + (energy-5.000000)* 0.000006) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 10.0000 && energy < 100.0000) * (0.004287 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005101 + (energy-100.000000)* 0.000008) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.012303 + (energy-1000.000000)* 0.000011) + \
       (abs(eta) >= 1.4000 && abs(eta) < 1.5000) * (energy >= 10000.0000) * (0.110091*energy/10000.000000) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 0.0000 && energy < 1.0000) * (0.00680449) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 1.0000 && energy < 2.0000) * (0.006804 + (energy-1.000000)* -0.002108) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 2.0000 && energy < 5.0000) * (0.004696 + (energy-2.000000)* -0.000116) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 5.0000 && energy < 10.0000) * (0.004348 + (energy-5.000000)* 0.000009) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 10.0000 && energy < 100.0000) * (0.004391 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005232 + (energy-100.000000)* 0.000007) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.011392 + (energy-1000.000000)* 0.000010) + \
       (abs(eta) >= 1.5000 && abs(eta) < 1.6000) * (energy >= 10000.0000) * (0.097988*energy/10000.000000) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 0.0000 && energy < 1.0000) * (0.00763793) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 1.0000 && energy < 2.0000) * (0.007638 + (energy-1.000000)* -0.003061) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 2.0000 && energy < 5.0000) * (0.004577 + (energy-2.000000)* -0.000168) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 5.0000 && energy < 10.0000) * (0.004074 + (energy-5.000000)* 0.000009) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 10.0000 && energy < 100.0000) * (0.004121 + (energy-10.000000)* 0.000009) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004910 + (energy-100.000000)* 0.000006) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.010635 + (energy-1000.000000)* 0.000009) + \
       (abs(eta) >= 1.6000 && abs(eta) < 1.7000) * (energy >= 10000.0000) * (0.090704*energy/10000.000000) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 0.0000 && energy < 1.0000) * (0.00913948) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 1.0000 && energy < 2.0000) * (0.009139 + (energy-1.000000)* -0.004405) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 2.0000 && energy < 5.0000) * (0.004735 + (energy-2.000000)* -0.000231) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 5.0000 && energy < 10.0000) * (0.004043 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 10.0000 && energy < 100.0000) * (0.004066 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004739 + (energy-100.000000)* 0.000006) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.009891 + (energy-1000.000000)* 0.000008) + \
       (abs(eta) >= 1.7000 && abs(eta) < 1.8000) * (energy >= 10000.0000) * (0.083216*energy/10000.000000) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 0.0000 && energy < 1.0000) * (0.00956747) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 1.0000 && energy < 2.0000) * (0.009567 + (energy-1.000000)* -0.004497) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 2.0000 && energy < 5.0000) * (0.005070 + (energy-2.000000)* -0.000331) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 5.0000 && energy < 10.0000) * (0.004078 + (energy-5.000000)* 0.000004) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 10.0000 && energy < 100.0000) * (0.004100 + (energy-10.000000)* 0.000006) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004682 + (energy-100.000000)* 0.000005) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.009157 + (energy-1000.000000)* 0.000007) + \
       (abs(eta) >= 1.8000 && abs(eta) < 1.9000) * (energy >= 10000.0000) * (0.074702*energy/10000.000000) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 0.0000 && energy < 1.0000) * (0.00964334) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 1.0000 && energy < 2.0000) * (0.009643 + (energy-1.000000)* -0.003950) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 2.0000 && energy < 5.0000) * (0.005694 + (energy-2.000000)* -0.000479) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 5.0000 && energy < 10.0000) * (0.004256 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 10.0000 && energy < 100.0000) * (0.004281 + (energy-10.000000)* 0.000006) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004862 + (energy-100.000000)* 0.000004) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.008501 + (energy-1000.000000)* 0.000006) + \
       (abs(eta) >= 1.9000 && abs(eta) < 2.0000) * (energy >= 10000.0000) * (0.062525*energy/10000.000000) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 0.0000 && energy < 1.0000) * (0.01045039) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 1.0000 && energy < 2.0000) * (0.010450 + (energy-1.000000)* -0.005379) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 2.0000 && energy < 5.0000) * (0.005072 + (energy-2.000000)* -0.000321) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 5.0000 && energy < 10.0000) * (0.004109 + (energy-5.000000)* 0.000006) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 10.0000 && energy < 100.0000) * (0.004137 + (energy-10.000000)* 0.000006) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004666 + (energy-100.000000)* 0.000005) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.008951 + (energy-1000.000000)* 0.000007) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 10000.0000) * (0.073400*energy/10000.000000) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 0.0000 && energy < 1.0000) * (0.01046694) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 1.0000 && energy < 2.0000) * (0.010467 + (energy-1.000000)* -0.005023) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 2.0000 && energy < 5.0000) * (0.005444 + (energy-2.000000)* -0.000330) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 5.0000 && energy < 10.0000) * (0.004455 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 10.0000 && energy < 100.0000) * (0.004479 + (energy-10.000000)* 0.000004) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004883 + (energy-100.000000)* 0.000005) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.009382 + (energy-1000.000000)* 0.000008) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 10000.0000) * (0.078852*energy/10000.000000) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 0.0000 && energy < 1.0000) * (0.01090933) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 1.0000 && energy < 2.0000) * (0.010909 + (energy-1.000000)* -0.005299) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 2.0000 && energy < 5.0000) * (0.005610 + (energy-2.000000)* -0.000302) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 5.0000 && energy < 10.0000) * (0.004704 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 10.0000 && energy < 100.0000) * (0.004730 + (energy-10.000000)* 0.000005) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005146 + (energy-100.000000)* 0.000006) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.010323 + (energy-1000.000000)* 0.000009) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 10000.0000) * (0.088469*energy/10000.000000) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 0.0000 && energy < 1.0000) * (0.01271833) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 1.0000 && energy < 2.0000) * (0.012718 + (energy-1.000000)* -0.005764) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 2.0000 && energy < 5.0000) * (0.006954 + (energy-2.000000)* -0.000492) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 5.0000 && energy < 10.0000) * (0.005479 + (energy-5.000000)* 0.000003) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 10.0000 && energy < 100.0000) * (0.005495 + (energy-10.000000)* 0.000003) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005774 + (energy-100.000000)* 0.000006) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.011265 + (energy-1000.000000)* 0.000009) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 10000.0000) * (0.096592*energy/10000.000000) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 0.0000 && energy < 1.0000) * (0.01515272) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 1.0000 && energy < 2.0000) * (0.015153 + (energy-1.000000)* -0.007272) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 2.0000 && energy < 5.0000) * (0.007881 + (energy-2.000000)* -0.000660) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 5.0000 && energy < 10.0000) * (0.005900 + (energy-5.000000)* 0.000003) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 10.0000 && energy < 100.0000) * (0.005914 + (energy-10.000000)* 0.000003) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 100.0000 && energy < 1000.0000) * (0.006174 + (energy-100.000000)* 0.000007) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.012486 + (energy-1000.000000)* 0.000011) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 10000.0000) * (0.108659*energy/10000.000000)
    }
}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
    set InputArray MuonTrackingEfficiency/muons
    set OutputArray muons

    # Resolution given in dpT/pT.

    set ResolutionFormula {

    ( abs(eta) < 1.35 ) *

    (sqrt(0.0000364164 + (
			  3*9.06262e-8 *pt^2* cosh(
         eta)^2 *(2.82074e-7/sin(2*atan(exp(-abs(eta))))^2 + (
          504.525 *(1/400000000 + (0.117945* 1/cosh(eta)^2)/(
             pt^2 *sin(2*atan(exp(-abs(eta))))^2)))/
          sin(2*atan(exp(-abs(eta))))^2) *sin(2*atan(exp(-abs(eta))))^2)/(
       0.00516429/sin(2*atan(exp(-abs(eta))))^2 + (
        96868.8 *(1/400000000 + 5*(0.117945 * 1/cosh(eta)^2)/(
           pt^2 *sin(2*atan(exp(-abs(eta))))^2)))/
        sin(2*atan(exp(-abs(eta))))^2))

    ) +

    ( abs(eta) > 1.35 && abs(eta) < 2.00) *

    ( 1.53846*(

     (
       (energy >= 0.0000 && energy < 1.0000) * (0.00953577) +
    (energy >= 1.0000 && energy < 2.0000) * (0.009536 + (energy-1.000000)* -0.003793) +
    (energy >= 2.0000 && energy < 5.0000) * (0.005742 + (energy-2.000000)* -0.000489) +
    (energy >= 5.0000 && energy < 10.0000) * (0.004277 + (energy-5.000000)* 0.000005) +
    (energy >= 10.0000 && energy < 20.0000) * (0.004302 + (energy-10.000000)* 0.000007) +
    (energy >= 20.0000 && energy < 50.0000) * (0.004368 + (energy-20.000000)* 0.000007) +
    (energy >= 50.0000 && energy < 100.0000) * (0.004581 + (energy-50.000000)* 0.000006) +
    (energy >= 100.0000 && energy < 200.0000) * (0.004875 + (energy-100.000000)* 0.000005) +
    (energy >= 200.0000 && energy < 500.0000) * (0.005344 + (energy-200.000000)* 0.000004) +
    (energy >= 500.0000 && energy < 1000.0000) * (0.006395 + (energy-500.000000)* 0.000004) +
    (energy >= 1000.0000 && energy < 2000.0000) * (0.008441 + (energy-1000.000000)* 0.000005) +
    (energy >= 2000.0000 && energy < 5000.0000) * (0.013700 + (energy-2000.000000)* 0.000006) +
    (energy >= 5000.0000 && energy < 10000.0000) * (0.031615 + (energy-5000.000000)* 0.000006) +
    (energy >= 10000.0000) * (0.062437*energy/10000.000000)
     )


     - sqrt(0.0000364164 + (9.06262e-8 *(1.19507e-6 + 2137.54 *(1/400000000 + 0.155982/pt^2)) *pt^2)/(0.0218797 + 410407. *(1/400000000 + 0.155982/pt^2))))*abs(eta) +

    3.07692*sqrt(0.0000364164 + (9.06262e-8 *(1.19507e-6 + 2137.54 *(1/400000000 + 0.155982/pt^2)) *pt^2)/(0.0218797 + 410407. *(1/400000000 + 0.155982/pt^2)))

     - 2.07692*
     (

      (energy >= 0.0000 && energy < 1.0000) * (0.00953577) +
    (energy >= 1.0000 && energy < 2.0000) * (0.009536 + (energy-1.000000)* -0.003793) +
    (energy >= 2.0000 && energy < 5.0000) * (0.005742 + (energy-2.000000)* -0.000489) +
    (energy >= 5.0000 && energy < 10.0000) * (0.004277 + (energy-5.000000)* 0.000005) +
    (energy >= 10.0000 && energy < 20.0000) * (0.004302 + (energy-10.000000)* 0.000007) +
    (energy >= 20.0000 && energy < 50.0000) * (0.004368 + (energy-20.000000)* 0.000007) +
    (energy >= 50.0000 && energy < 100.0000) * (0.004581 + (energy-50.000000)* 0.000006) +
    (energy >= 100.0000 && energy < 200.0000) * (0.004875 + (energy-100.000000)* 0.000005) +
    (energy >= 200.0000 && energy < 500.0000) * (0.005344 + (energy-200.000000)* 0.000004) +
    (energy >= 500.0000 && energy < 1000.0000) * (0.006395 + (energy-500.000000)* 0.000004) +
    (energy >= 1000.0000 && energy < 2000.0000) * (0.008441 + (energy-1000.000000)* 0.000005) +
    (energy >= 2000.0000 && energy < 5000.0000) * (0.013700 + (energy-2000.000000)* 0.000006) +
    (energy >= 5000.0000 && energy < 10000.0000) * (0.031615 + (energy-5000.000000)* 0.000006) +
    (energy >= 10000.0000) * (0.062437*energy/10000.000000)

     )


    ) +

       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 0.0000 && energy < 1.0000) * (0.01062416) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 1.0000 && energy < 2.0000) * (0.010624 + (energy-1.000000)* -0.005532) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 2.0000 && energy < 5.0000) * (0.005092 + (energy-2.000000)* -0.000326) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 5.0000 && energy < 10.0000) * (0.004115 + (energy-5.000000)* 0.000006) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 10.0000 && energy < 20.0000) * (0.004143 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 20.0000 && energy < 50.0000) * (0.004209 + (energy-20.000000)* 0.000007) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 50.0000 && energy < 100.0000) * (0.004413 + (energy-50.000000)* 0.000005) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 100.0000 && energy < 200.0000) * (0.004681 + (energy-100.000000)* 0.000004) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 200.0000 && energy < 500.0000) * (0.005076 + (energy-200.000000)* 0.000004) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 500.0000 && energy < 1000.0000) * (0.006270 + (energy-500.000000)* 0.000005) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 1000.0000 && energy < 2000.0000) * (0.008960 + (energy-1000.000000)* 0.000007) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 2000.0000 && energy < 5000.0000) * (0.015510 + (energy-2000.000000)* 0.000007) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 5000.0000 && energy < 10000.0000) * (0.036867 + (energy-5000.000000)* 0.000007) + \
       (abs(eta) >= 2.0000 && abs(eta) < 2.1000) * (energy >= 10000.0000) * (0.073168*energy/10000.000000) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 0.0000 && energy < 1.0000) * (0.01007098) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 1.0000 && energy < 2.0000) * (0.010071 + (energy-1.000000)* -0.004627) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 2.0000 && energy < 5.0000) * (0.005444 + (energy-2.000000)* -0.000322) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 5.0000 && energy < 10.0000) * (0.004478 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 10.0000 && energy < 20.0000) * (0.004501 + (energy-10.000000)* 0.000006) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 20.0000 && energy < 50.0000) * (0.004558 + (energy-20.000000)* 0.000005) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 50.0000 && energy < 100.0000) * (0.004701 + (energy-50.000000)* 0.000004) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 100.0000 && energy < 200.0000) * (0.004888 + (energy-100.000000)* 0.000003) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 200.0000 && energy < 500.0000) * (0.005213 + (energy-200.000000)* 0.000004) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 500.0000 && energy < 1000.0000) * (0.006429 + (energy-500.000000)* 0.000006) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 1000.0000 && energy < 2000.0000) * (0.009343 + (energy-1000.000000)* 0.000007) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 2000.0000 && energy < 5000.0000) * (0.016410 + (energy-2000.000000)* 0.000008) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 5000.0000 && energy < 10000.0000) * (0.039265 + (energy-5000.000000)* 0.000008) + \
       (abs(eta) >= 2.1000 && abs(eta) < 2.2000) * (energy >= 10000.0000) * (0.078014*energy/10000.000000) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 0.0000 && energy < 1.0000) * (0.01095892) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 1.0000 && energy < 2.0000) * (0.010959 + (energy-1.000000)* -0.005458) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 2.0000 && energy < 5.0000) * (0.005501 + (energy-2.000000)* -0.000281) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 5.0000 && energy < 10.0000) * (0.004660 + (energy-5.000000)* 0.000005) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 10.0000 && energy < 20.0000) * (0.004686 + (energy-10.000000)* 0.000007) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 20.0000 && energy < 50.0000) * (0.004757 + (energy-20.000000)* 0.000006) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 50.0000 && energy < 100.0000) * (0.004937 + (energy-50.000000)* 0.000004) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 100.0000 && energy < 200.0000) * (0.005143 + (energy-100.000000)* 0.000004) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 200.0000 && energy < 500.0000) * (0.005505 + (energy-200.000000)* 0.000005) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 500.0000 && energy < 1000.0000) * (0.006975 + (energy-500.000000)* 0.000007) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 1000.0000 && energy < 2000.0000) * (0.010462 + (energy-1000.000000)* 0.000008) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 2000.0000 && energy < 5000.0000) * (0.018731 + (energy-2000.000000)* 0.000009) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 5000.0000 && energy < 10000.0000) * (0.045159 + (energy-5000.000000)* 0.000009) + \
       (abs(eta) >= 2.2000 && abs(eta) < 2.3000) * (energy >= 10000.0000) * (0.089830*energy/10000.000000) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 0.0000 && energy < 1.0000) * (0.01279214) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 1.0000 && energy < 2.0000) * (0.012792 + (energy-1.000000)* -0.005763) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 2.0000 && energy < 5.0000) * (0.007029 + (energy-2.000000)* -0.000513) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 5.0000 && energy < 10.0000) * (0.005489 + (energy-5.000000)* 0.000003) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 10.0000 && energy < 20.0000) * (0.005503 + (energy-10.000000)* 0.000003) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 20.0000 && energy < 50.0000) * (0.005537 + (energy-20.000000)* 0.000003) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 50.0000 && energy < 100.0000) * (0.005636 + (energy-50.000000)* 0.000003) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 100.0000 && energy < 200.0000) * (0.005773 + (energy-100.000000)* 0.000003) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 200.0000 && energy < 500.0000) * (0.006060 + (energy-200.000000)* 0.000005) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 500.0000 && energy < 1000.0000) * (0.007489 + (energy-500.000000)* 0.000007) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 1000.0000 && energy < 2000.0000) * (0.011117 + (energy-1000.000000)* 0.000009) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 2000.0000 && energy < 5000.0000) * (0.019824 + (energy-2000.000000)* 0.000009) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 5000.0000 && energy < 10000.0000) * (0.047732 + (energy-5000.000000)* 0.000009) + \
       (abs(eta) >= 2.3000 && abs(eta) < 2.4000) * (energy >= 10000.0000) * (0.094931*energy/10000.000000) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 0.0000 && energy < 1.0000) * (0.01502671) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 1.0000 && energy < 2.0000) * (0.015027 + (energy-1.000000)* -0.007177) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 2.0000 && energy < 5.0000) * (0.007850 + (energy-2.000000)* -0.000651) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 5.0000 && energy < 10.0000) * (0.005898 + (energy-5.000000)* 0.000003) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 10.0000 && energy < 20.0000) * (0.005913 + (energy-10.000000)* 0.000003) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 20.0000 && energy < 50.0000) * (0.005947 + (energy-20.000000)* 0.000003) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 50.0000 && energy < 100.0000) * (0.006039 + (energy-50.000000)* 0.000003) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 100.0000 && energy < 200.0000) * (0.006170 + (energy-100.000000)* 0.000003) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 200.0000 && energy < 500.0000) * (0.006485 + (energy-200.000000)* 0.000006) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 500.0000 && energy < 1000.0000) * (0.008140 + (energy-500.000000)* 0.000008) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 1000.0000 && energy < 2000.0000) * (0.012304 + (energy-1000.000000)* 0.000010) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 2000.0000 && energy < 5000.0000) * (0.022168 + (energy-2000.000000)* 0.000010) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 5000.0000 && energy < 10000.0000) * (0.053585 + (energy-5000.000000)* 0.000011) + \
       (abs(eta) >= 2.4000 && abs(eta) < 2.5000) * (energy >= 10000.0000) * (0.106635*energy/10000.000000)
    }
}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing ForwardMuonMomentumSmearing {
    set InputArray ForwardMuonEfficiency/muons
    set OutputArray muons

    # Resolution given in dpT/pT (assume forward tracking spectrometer from FCC-hh)

    set ResolutionFormula {
    (abs(eta) >= 2.5000 && abs(eta) < 2.6000) * (energy >= 0.0000 && energy < 1.0000) * (0.01778926) + \
    (abs(eta) >= 2.5000 && abs(eta) < 2.6000) * (energy >= 1.0000 && energy < 2.0000) * (0.017789 + (energy-1.000000)* -0.008811) + \
    (abs(eta) >= 2.5000 && abs(eta) < 2.6000) * (energy >= 2.0000 && energy < 5.0000) * (0.008978 + (energy-2.000000)* -0.000844) + \
    (abs(eta) >= 2.5000 && abs(eta) < 2.6000) * (energy >= 5.0000 && energy < 10.0000) * (0.006446 + (energy-5.000000)* -0.000195) + \
    (abs(eta) >= 2.5000 && abs(eta) < 2.6000) * (energy >= 10.0000 && energy < 100.0000) * (0.005473 + (energy-10.000000)* 0.000003) + \
    (abs(eta) >= 2.5000 && abs(eta) < 2.6000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005782 + (energy-100.000000)* 0.000005) + \
    (abs(eta) >= 2.5000 && abs(eta) < 2.6000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.010655 + (energy-1000.000000)* 0.000009) + \
    (abs(eta) >= 2.5000 && abs(eta) < 2.6000) * (energy >= 10000.0000) * (0.087190*energy/10000.000000) + \
    (abs(eta) >= 2.6000 && abs(eta) < 2.7000) * (energy >= 0.0000 && energy < 1.0000) * (0.02061617) + \
    (abs(eta) >= 2.6000 && abs(eta) < 2.7000) * (energy >= 1.0000 && energy < 2.0000) * (0.020616 + (energy-1.000000)* -0.010732) + \
    (abs(eta) >= 2.6000 && abs(eta) < 2.7000) * (energy >= 2.0000 && energy < 5.0000) * (0.009884 + (energy-2.000000)* -0.001349) + \
    (abs(eta) >= 2.6000 && abs(eta) < 2.7000) * (energy >= 5.0000 && energy < 10.0000) * (0.005836 + (energy-5.000000)* -0.000484) + \
    (abs(eta) >= 2.6000 && abs(eta) < 2.7000) * (energy >= 10.0000 && energy < 100.0000) * (0.003417 + (energy-10.000000)* 0.000003) + \
    (abs(eta) >= 2.6000 && abs(eta) < 2.7000) * (energy >= 100.0000 && energy < 1000.0000) * (0.003655 + (energy-100.000000)* 0.000002) + \
    (abs(eta) >= 2.6000 && abs(eta) < 2.7000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.005164 + (energy-1000.000000)* 0.000003) + \
    (abs(eta) >= 2.6000 && abs(eta) < 2.7000) * (energy >= 10000.0000) * (0.029564*energy/10000.000000) + \
    (abs(eta) >= 2.7000 && abs(eta) < 2.8000) * (energy >= 0.0000 && energy < 1.0000) * (0.02328940) + \
    (abs(eta) >= 2.7000 && abs(eta) < 2.8000) * (energy >= 1.0000 && energy < 2.0000) * (0.023289 + (energy-1.000000)* -0.012327) + \
    (abs(eta) >= 2.7000 && abs(eta) < 2.8000) * (energy >= 2.0000 && energy < 5.0000) * (0.010962 + (energy-2.000000)* -0.001777) + \
    (abs(eta) >= 2.7000 && abs(eta) < 2.8000) * (energy >= 5.0000 && energy < 10.0000) * (0.005632 + (energy-5.000000)* -0.000436) + \
    (abs(eta) >= 2.7000 && abs(eta) < 2.8000) * (energy >= 10.0000 && energy < 100.0000) * (0.003452 + (energy-10.000000)* 0.000001) + \
    (abs(eta) >= 2.7000 && abs(eta) < 2.8000) * (energy >= 100.0000 && energy < 1000.0000) * (0.003539 + (energy-100.000000)* 0.000001) + \
    (abs(eta) >= 2.7000 && abs(eta) < 2.8000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.004876 + (energy-1000.000000)* 0.000002) + \
    (abs(eta) >= 2.7000 && abs(eta) < 2.8000) * (energy >= 10000.0000) * (0.026660*energy/10000.000000) + \
    (abs(eta) >= 2.8000 && abs(eta) < 2.9000) * (energy >= 0.0000 && energy < 1.0000) * (0.02303055) + \
    (abs(eta) >= 2.8000 && abs(eta) < 2.9000) * (energy >= 1.0000 && energy < 2.0000) * (0.023031 + (energy-1.000000)* -0.011212) + \
    (abs(eta) >= 2.8000 && abs(eta) < 2.9000) * (energy >= 2.0000 && energy < 5.0000) * (0.011819 + (energy-2.000000)* -0.002216) + \
    (abs(eta) >= 2.8000 && abs(eta) < 2.9000) * (energy >= 5.0000 && energy < 10.0000) * (0.005171 + (energy-5.000000)* -0.000352) + \
    (abs(eta) >= 2.8000 && abs(eta) < 2.9000) * (energy >= 10.0000 && energy < 100.0000) * (0.003412 + (energy-10.000000)* 0.000001) + \
    (abs(eta) >= 2.8000 && abs(eta) < 2.9000) * (energy >= 100.0000 && energy < 1000.0000) * (0.003516 + (energy-100.000000)* 0.000001) + \
    (abs(eta) >= 2.8000 && abs(eta) < 2.9000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.004565 + (energy-1000.000000)* 0.000002) + \
    (abs(eta) >= 2.8000 && abs(eta) < 2.9000) * (energy >= 10000.0000) * (0.025158*energy/10000.000000) + \
    (abs(eta) >= 2.9000 && abs(eta) < 3.0000) * (energy >= 0.0000 && energy < 1.0000) * (0.02611889) + \
    (abs(eta) >= 2.9000 && abs(eta) < 3.0000) * (energy >= 1.0000 && energy < 2.0000) * (0.026119 + (energy-1.000000)* -0.012835) + \
    (abs(eta) >= 2.9000 && abs(eta) < 3.0000) * (energy >= 2.0000 && energy < 5.0000) * (0.013284 + (energy-2.000000)* -0.002589) + \
    (abs(eta) >= 2.9000 && abs(eta) < 3.0000) * (energy >= 5.0000 && energy < 10.0000) * (0.005516 + (energy-5.000000)* -0.000388) + \
    (abs(eta) >= 2.9000 && abs(eta) < 3.0000) * (energy >= 10.0000 && energy < 100.0000) * (0.003574 + (energy-10.000000)* 0.000001) + \
    (abs(eta) >= 2.9000 && abs(eta) < 3.0000) * (energy >= 100.0000 && energy < 1000.0000) * (0.003690 + (energy-100.000000)* 0.000001) + \
    (abs(eta) >= 2.9000 && abs(eta) < 3.0000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.004493 + (energy-1000.000000)* 0.000002) + \
    (abs(eta) >= 2.9000 && abs(eta) < 3.0000) * (energy >= 10000.0000) * (0.022856*energy/10000.000000) + \
    (abs(eta) >= 3.0000 && abs(eta) < 3.1000) * (energy >= 0.0000 && energy < 1.0000) * (0.02991154) + \
    (abs(eta) >= 3.0000 && abs(eta) < 3.1000) * (energy >= 1.0000 && energy < 2.0000) * (0.029912 + (energy-1.000000)* -0.015226) + \
    (abs(eta) >= 3.0000 && abs(eta) < 3.1000) * (energy >= 2.0000 && energy < 5.0000) * (0.014686 + (energy-2.000000)* -0.002867) + \
    (abs(eta) >= 3.0000 && abs(eta) < 3.1000) * (energy >= 5.0000 && energy < 10.0000) * (0.006086 + (energy-5.000000)* -0.000457) + \
    (abs(eta) >= 3.0000 && abs(eta) < 3.1000) * (energy >= 10.0000 && energy < 100.0000) * (0.003803 + (energy-10.000000)* 0.000001) + \
    (abs(eta) >= 3.0000 && abs(eta) < 3.1000) * (energy >= 100.0000 && energy < 1000.0000) * (0.003932 + (energy-100.000000)* 0.000001) + \
    (abs(eta) >= 3.0000 && abs(eta) < 3.1000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.004619 + (energy-1000.000000)* 0.000002) + \
    (abs(eta) >= 3.0000 && abs(eta) < 3.1000) * (energy >= 10000.0000) * (0.020643*energy/10000.000000) + \
    (abs(eta) >= 3.1000 && abs(eta) < 3.2000) * (energy >= 0.0000 && energy < 1.0000) * (0.03556765) + \
    (abs(eta) >= 3.1000 && abs(eta) < 3.2000) * (energy >= 1.0000 && energy < 2.0000) * (0.035568 + (energy-1.000000)* -0.019079) + \
    (abs(eta) >= 3.1000 && abs(eta) < 3.2000) * (energy >= 2.0000 && energy < 5.0000) * (0.016488 + (energy-2.000000)* -0.003242) + \
    (abs(eta) >= 3.1000 && abs(eta) < 3.2000) * (energy >= 5.0000 && energy < 10.0000) * (0.006763 + (energy-5.000000)* -0.000526) + \
    (abs(eta) >= 3.1000 && abs(eta) < 3.2000) * (energy >= 10.0000 && energy < 100.0000) * (0.004133 + (energy-10.000000)* 0.000002) + \
    (abs(eta) >= 3.1000 && abs(eta) < 3.2000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004270 + (energy-100.000000)* 0.000001) + \
    (abs(eta) >= 3.1000 && abs(eta) < 3.2000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.004991 + (energy-1000.000000)* 0.000002) + \
    (abs(eta) >= 3.1000 && abs(eta) < 3.2000) * (energy >= 10000.0000) * (0.021331*energy/10000.000000) + \
    (abs(eta) >= 3.2000 && abs(eta) < 3.3000) * (energy >= 0.0000 && energy < 1.0000) * (0.03803690) + \
    (abs(eta) >= 3.2000 && abs(eta) < 3.3000) * (energy >= 1.0000 && energy < 2.0000) * (0.038037 + (energy-1.000000)* -0.020041) + \
    (abs(eta) >= 3.2000 && abs(eta) < 3.3000) * (energy >= 2.0000 && energy < 5.0000) * (0.017996 + (energy-2.000000)* -0.003460) + \
    (abs(eta) >= 3.2000 && abs(eta) < 3.3000) * (energy >= 5.0000 && energy < 10.0000) * (0.007617 + (energy-5.000000)* -0.000724) + \
    (abs(eta) >= 3.2000 && abs(eta) < 3.3000) * (energy >= 10.0000 && energy < 100.0000) * (0.003999 + (energy-10.000000)* 0.000002) + \
    (abs(eta) >= 3.2000 && abs(eta) < 3.3000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004159 + (energy-100.000000)* 0.000001) + \
    (abs(eta) >= 3.2000 && abs(eta) < 3.3000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.005300 + (energy-1000.000000)* 0.000002) + \
    (abs(eta) >= 3.2000 && abs(eta) < 3.3000) * (energy >= 10000.0000) * (0.023674*energy/10000.000000) + \
    (abs(eta) >= 3.3000 && abs(eta) < 3.4000) * (energy >= 0.0000 && energy < 1.0000) * (0.03760924) + \
    (abs(eta) >= 3.3000 && abs(eta) < 3.4000) * (energy >= 1.0000 && energy < 2.0000) * (0.037609 + (energy-1.000000)* -0.018328) + \
    (abs(eta) >= 3.3000 && abs(eta) < 3.4000) * (energy >= 2.0000 && energy < 5.0000) * (0.019281 + (energy-2.000000)* -0.003682) + \
    (abs(eta) >= 3.3000 && abs(eta) < 3.4000) * (energy >= 5.0000 && energy < 10.0000) * (0.008234 + (energy-5.000000)* -0.000772) + \
    (abs(eta) >= 3.3000 && abs(eta) < 3.4000) * (energy >= 10.0000 && energy < 100.0000) * (0.004372 + (energy-10.000000)* 0.000002) + \
    (abs(eta) >= 3.3000 && abs(eta) < 3.4000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004545 + (energy-100.000000)* 0.000001) + \
    (abs(eta) >= 3.3000 && abs(eta) < 3.4000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.005786 + (energy-1000.000000)* 0.000002) + \
    (abs(eta) >= 3.3000 && abs(eta) < 3.4000) * (energy >= 10000.0000) * (0.026924*energy/10000.000000) + \
    (abs(eta) >= 3.4000 && abs(eta) < 3.5000) * (energy >= 0.0000 && energy < 1.0000) * (0.03763696) + \
    (abs(eta) >= 3.4000 && abs(eta) < 3.5000) * (energy >= 1.0000 && energy < 2.0000) * (0.037637 + (energy-1.000000)* -0.016962) + \
    (abs(eta) >= 3.4000 && abs(eta) < 3.5000) * (energy >= 2.0000 && energy < 5.0000) * (0.020675 + (energy-2.000000)* -0.003893) + \
    (abs(eta) >= 3.4000 && abs(eta) < 3.5000) * (energy >= 5.0000 && energy < 10.0000) * (0.008995 + (energy-5.000000)* -0.000848) + \
    (abs(eta) >= 3.4000 && abs(eta) < 3.5000) * (energy >= 10.0000 && energy < 100.0000) * (0.004755 + (energy-10.000000)* 0.000002) + \
    (abs(eta) >= 3.4000 && abs(eta) < 3.5000) * (energy >= 100.0000 && energy < 1000.0000) * (0.004949 + (energy-100.000000)* 0.000002) + \
    (abs(eta) >= 3.4000 && abs(eta) < 3.5000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.006336 + (energy-1000.000000)* 0.000003) + \
    (abs(eta) >= 3.4000 && abs(eta) < 3.5000) * (energy >= 10000.0000) * (0.030130*energy/10000.000000) + \
    (abs(eta) >= 3.5000 && abs(eta) < 3.6000) * (energy >= 0.0000 && energy < 1.0000) * (0.04368466) + \
    (abs(eta) >= 3.5000 && abs(eta) < 3.6000) * (energy >= 1.0000 && energy < 2.0000) * (0.043685 + (energy-1.000000)* -0.020430) + \
    (abs(eta) >= 3.5000 && abs(eta) < 3.6000) * (energy >= 2.0000 && energy < 5.0000) * (0.023254 + (energy-2.000000)* -0.004385) + \
    (abs(eta) >= 3.5000 && abs(eta) < 3.6000) * (energy >= 5.0000 && energy < 10.0000) * (0.010100 + (energy-5.000000)* -0.001016) + \
    (abs(eta) >= 3.5000 && abs(eta) < 3.6000) * (energy >= 10.0000 && energy < 100.0000) * (0.005020 + (energy-10.000000)* 0.000002) + \
    (abs(eta) >= 3.5000 && abs(eta) < 3.6000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005206 + (energy-100.000000)* 0.000002) + \
    (abs(eta) >= 3.5000 && abs(eta) < 3.6000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.006679 + (energy-1000.000000)* 0.000003) + \
    (abs(eta) >= 3.5000 && abs(eta) < 3.6000) * (energy >= 10000.0000) * (0.033540*energy/10000.000000) + \
    (abs(eta) >= 3.6000 && abs(eta) < 3.7000) * (energy >= 0.0000 && energy < 1.0000) * (0.05055352) + \
    (abs(eta) >= 3.6000 && abs(eta) < 3.7000) * (energy >= 1.0000 && energy < 2.0000) * (0.050554 + (energy-1.000000)* -0.024554) + \
    (abs(eta) >= 3.6000 && abs(eta) < 3.7000) * (energy >= 2.0000 && energy < 5.0000) * (0.025999 + (energy-2.000000)* -0.004964) + \
    (abs(eta) >= 3.6000 && abs(eta) < 3.7000) * (energy >= 5.0000 && energy < 10.0000) * (0.011106 + (energy-5.000000)* -0.001134) + \
    (abs(eta) >= 3.6000 && abs(eta) < 3.7000) * (energy >= 10.0000 && energy < 100.0000) * (0.005436 + (energy-10.000000)* 0.000002) + \
    (abs(eta) >= 3.6000 && abs(eta) < 3.7000) * (energy >= 100.0000 && energy < 1000.0000) * (0.005582 + (energy-100.000000)* 0.000002) + \
    (abs(eta) >= 3.6000 && abs(eta) < 3.7000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.007047 + (energy-1000.000000)* 0.000003) + \
    (abs(eta) >= 3.6000 && abs(eta) < 3.7000) * (energy >= 10000.0000) * (0.036533*energy/10000.000000) + \
    (abs(eta) >= 3.7000 && abs(eta) < 3.8000) * (energy >= 0.0000 && energy < 1.0000) * (0.05790950) + \
    (abs(eta) >= 3.7000 && abs(eta) < 3.8000) * (energy >= 1.0000 && energy < 2.0000) * (0.057910 + (energy-1.000000)* -0.028845) + \
    (abs(eta) >= 3.7000 && abs(eta) < 3.8000) * (energy >= 2.0000 && energy < 5.0000) * (0.029065 + (energy-2.000000)* -0.005581) + \
    (abs(eta) >= 3.7000 && abs(eta) < 3.8000) * (energy >= 5.0000 && energy < 10.0000) * (0.012322 + (energy-5.000000)* -0.001294) + \
    (abs(eta) >= 3.7000 && abs(eta) < 3.8000) * (energy >= 10.0000 && energy < 100.0000) * (0.005853 + (energy-10.000000)* 0.000002) + \
    (abs(eta) >= 3.7000 && abs(eta) < 3.8000) * (energy >= 100.0000 && energy < 1000.0000) * (0.006005 + (energy-100.000000)* 0.000002) + \
    (abs(eta) >= 3.7000 && abs(eta) < 3.8000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.007594 + (energy-1000.000000)* 0.000004) + \
    (abs(eta) >= 3.7000 && abs(eta) < 3.8000) * (energy >= 10000.0000) * (0.040788*energy/10000.000000) + \
    (abs(eta) >= 3.8000 && abs(eta) < 3.9000) * (energy >= 0.0000 && energy < 1.0000) * (0.06787795) + \
    (abs(eta) >= 3.8000 && abs(eta) < 3.9000) * (energy >= 1.0000 && energy < 2.0000) * (0.067878 + (energy-1.000000)* -0.035271) + \
    (abs(eta) >= 3.8000 && abs(eta) < 3.9000) * (energy >= 2.0000 && energy < 5.0000) * (0.032607 + (energy-2.000000)* -0.006346) + \
    (abs(eta) >= 3.8000 && abs(eta) < 3.9000) * (energy >= 5.0000 && energy < 10.0000) * (0.013570 + (energy-5.000000)* -0.001453) + \
    (abs(eta) >= 3.8000 && abs(eta) < 3.9000) * (energy >= 10.0000 && energy < 100.0000) * (0.006304 + (energy-10.000000)* 0.000002) + \
    (abs(eta) >= 3.8000 && abs(eta) < 3.9000) * (energy >= 100.0000 && energy < 1000.0000) * (0.006489 + (energy-100.000000)* 0.000002) + \
    (abs(eta) >= 3.8000 && abs(eta) < 3.9000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.008302 + (energy-1000.000000)* 0.000004) + \
    (abs(eta) >= 3.8000 && abs(eta) < 3.9000) * (energy >= 10000.0000) * (0.045331*energy/10000.000000) + \
    (abs(eta) >= 3.9000 && abs(eta) < 4.0000) * (energy >= 0.0000 && energy < 1.0000) * (0.08769401) + \
    (abs(eta) >= 3.9000 && abs(eta) < 4.0000) * (energy >= 1.0000 && energy < 2.0000) * (0.087694 + (energy-1.000000)* -0.050535) + \
    (abs(eta) >= 3.9000 && abs(eta) < 4.0000) * (energy >= 2.0000 && energy < 5.0000) * (0.037159 + (energy-2.000000)* -0.007293) + \
    (abs(eta) >= 3.9000 && abs(eta) < 4.0000) * (energy >= 5.0000 && energy < 10.0000) * (0.015278 + (energy-5.000000)* -0.001678) + \
    (abs(eta) >= 3.9000 && abs(eta) < 4.0000) * (energy >= 10.0000 && energy < 100.0000) * (0.006886 + (energy-10.000000)* 0.000003) + \
    (abs(eta) >= 3.9000 && abs(eta) < 4.0000) * (energy >= 100.0000 && energy < 1000.0000) * (0.007116 + (energy-100.000000)* 0.000002) + \
    (abs(eta) >= 3.9000 && abs(eta) < 4.0000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.009229 + (energy-1000.000000)* 0.000005) + \
    (abs(eta) >= 3.9000 && abs(eta) < 4.0000) * (energy >= 10000.0000) * (0.051460*energy/10000.000000) + \
    (abs(eta) >= 4.0000 && abs(eta) < 4.1000) * (energy >= 0.0000 && energy < 1.0000) * (0.12378717) + \
    (abs(eta) >= 4.0000 && abs(eta) < 4.1000) * (energy >= 1.0000 && energy < 2.0000) * (0.123787 + (energy-1.000000)* -0.081534) + \
    (abs(eta) >= 4.0000 && abs(eta) < 4.1000) * (energy >= 2.0000 && energy < 5.0000) * (0.042253 + (energy-2.000000)* -0.008504) + \
    (abs(eta) >= 4.0000 && abs(eta) < 4.1000) * (energy >= 5.0000 && energy < 10.0000) * (0.016742 + (energy-5.000000)* -0.001830) + \
    (abs(eta) >= 4.0000 && abs(eta) < 4.1000) * (energy >= 10.0000 && energy < 100.0000) * (0.007595 + (energy-10.000000)* 0.000003) + \
    (abs(eta) >= 4.0000 && abs(eta) < 4.1000) * (energy >= 100.0000 && energy < 1000.0000) * (0.007850 + (energy-100.000000)* 0.000003) + \
    (abs(eta) >= 4.0000 && abs(eta) < 4.1000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.010277 + (energy-1000.000000)* 0.000005) + \
    (abs(eta) >= 4.0000 && abs(eta) < 4.1000) * (energy >= 10000.0000) * (0.059714*energy/10000.000000) + \
    (abs(eta) >= 4.1000 && abs(eta) < 4.2000) * (energy >= 0.0000 && energy < 1.0000) * (0.14293546) + \
    (abs(eta) >= 4.1000 && abs(eta) < 4.2000) * (energy >= 1.0000 && energy < 2.0000) * (0.142935 + (energy-1.000000)* -0.095632) + \
    (abs(eta) >= 4.1000 && abs(eta) < 4.2000) * (energy >= 2.0000 && energy < 5.0000) * (0.047303 + (energy-2.000000)* -0.009550) + \
    (abs(eta) >= 4.1000 && abs(eta) < 4.2000) * (energy >= 5.0000 && energy < 10.0000) * (0.018654 + (energy-5.000000)* -0.002035) + \
    (abs(eta) >= 4.1000 && abs(eta) < 4.2000) * (energy >= 10.0000 && energy < 100.0000) * (0.008480 + (energy-10.000000)* 0.000003) + \
    (abs(eta) >= 4.1000 && abs(eta) < 4.2000) * (energy >= 100.0000 && energy < 1000.0000) * (0.008767 + (energy-100.000000)* 0.000003) + \
    (abs(eta) >= 4.1000 && abs(eta) < 4.2000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.011459 + (energy-1000.000000)* 0.000006) + \
    (abs(eta) >= 4.1000 && abs(eta) < 4.2000) * (energy >= 10000.0000) * (0.066014*energy/10000.000000) + \
    (abs(eta) >= 4.2000 && abs(eta) < 4.3000) * (energy >= 0.0000 && energy < 1.0000) * (0.16175578) + \
    (abs(eta) >= 4.2000 && abs(eta) < 4.3000) * (energy >= 1.0000 && energy < 2.0000) * (0.161756 + (energy-1.000000)* -0.108957) + \
    (abs(eta) >= 4.2000 && abs(eta) < 4.3000) * (energy >= 2.0000 && energy < 5.0000) * (0.052799 + (energy-2.000000)* -0.010753) + \
    (abs(eta) >= 4.2000 && abs(eta) < 4.3000) * (energy >= 5.0000 && energy < 10.0000) * (0.020541 + (energy-5.000000)* -0.002235) + \
    (abs(eta) >= 4.2000 && abs(eta) < 4.3000) * (energy >= 10.0000 && energy < 100.0000) * (0.009366 + (energy-10.000000)* 0.000004) + \
    (abs(eta) >= 4.2000 && abs(eta) < 4.3000) * (energy >= 100.0000 && energy < 1000.0000) * (0.009686 + (energy-100.000000)* 0.000003) + \
    (abs(eta) >= 4.2000 && abs(eta) < 4.3000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.012672 + (energy-1000.000000)* 0.000007) + \
    (abs(eta) >= 4.2000 && abs(eta) < 4.3000) * (energy >= 10000.0000) * (0.072840*energy/10000.000000) + \
    (abs(eta) >= 4.3000 && abs(eta) < 4.4000) * (energy >= 0.0000 && energy < 1.0000) * (0.17768991) + \
    (abs(eta) >= 4.3000 && abs(eta) < 4.4000) * (energy >= 1.0000 && energy < 2.0000) * (0.177690 + (energy-1.000000)* -0.119322) + \
    (abs(eta) >= 4.3000 && abs(eta) < 4.4000) * (energy >= 2.0000 && energy < 5.0000) * (0.058368 + (energy-2.000000)* -0.011801) + \
    (abs(eta) >= 4.3000 && abs(eta) < 4.4000) * (energy >= 5.0000 && energy < 10.0000) * (0.022966 + (energy-5.000000)* -0.002532) + \
    (abs(eta) >= 4.3000 && abs(eta) < 4.4000) * (energy >= 10.0000 && energy < 100.0000) * (0.010308 + (energy-10.000000)* 0.000004) + \
    (abs(eta) >= 4.3000 && abs(eta) < 4.4000) * (energy >= 100.0000 && energy < 1000.0000) * (0.010654 + (energy-100.000000)* 0.000004) + \
    (abs(eta) >= 4.3000 && abs(eta) < 4.4000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.013982 + (energy-1000.000000)* 0.000007) + \
    (abs(eta) >= 4.3000 && abs(eta) < 4.4000) * (energy >= 10000.0000) * (0.081143*energy/10000.000000) + \
    (abs(eta) >= 4.4000 && abs(eta) < 4.5000) * (energy >= 0.0000 && energy < 1.0000) * (0.19156432) + \
    (abs(eta) >= 4.4000 && abs(eta) < 4.5000) * (energy >= 1.0000 && energy < 2.0000) * (0.191564 + (energy-1.000000)* -0.125357) + \
    (abs(eta) >= 4.4000 && abs(eta) < 4.5000) * (energy >= 2.0000 && energy < 5.0000) * (0.066207 + (energy-2.000000)* -0.013680) + \
    (abs(eta) >= 4.4000 && abs(eta) < 4.5000) * (energy >= 5.0000 && energy < 10.0000) * (0.025168 + (energy-5.000000)* -0.002755) + \
    (abs(eta) >= 4.4000 && abs(eta) < 4.5000) * (energy >= 10.0000 && energy < 100.0000) * (0.011394 + (energy-10.000000)* 0.000004) + \
    (abs(eta) >= 4.4000 && abs(eta) < 4.5000) * (energy >= 100.0000 && energy < 1000.0000) * (0.011796 + (energy-100.000000)* 0.000004) + \
    (abs(eta) >= 4.4000 && abs(eta) < 4.5000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.015576 + (energy-1000.000000)* 0.000009) + \
    (abs(eta) >= 4.4000 && abs(eta) < 4.5000) * (energy >= 10000.0000) * (0.092193*energy/10000.000000) + \
    (abs(eta) >= 4.5000 && abs(eta) < 4.6000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 4.5000 && abs(eta) < 4.6000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.075652) + \
    (abs(eta) >= 4.5000 && abs(eta) < 4.6000) * (energy >= 2.0000 && energy < 5.0000) * (0.075652 + (energy-2.000000)* -0.015621) + \
    (abs(eta) >= 4.5000 && abs(eta) < 4.6000) * (energy >= 5.0000 && energy < 10.0000) * (0.028790 + (energy-5.000000)* -0.003236) + \
    (abs(eta) >= 4.5000 && abs(eta) < 4.6000) * (energy >= 10.0000 && energy < 100.0000) * (0.012611 + (energy-10.000000)* 0.000005) + \
    (abs(eta) >= 4.5000 && abs(eta) < 4.6000) * (energy >= 100.0000 && energy < 1000.0000) * (0.013037 + (energy-100.000000)* 0.000005) + \
    (abs(eta) >= 4.5000 && abs(eta) < 4.6000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.017269 + (energy-1000.000000)* 0.000010) + \
    (abs(eta) >= 4.5000 && abs(eta) < 4.6000) * (energy >= 10000.0000) * (0.103252*energy/10000.000000) + \
    (abs(eta) >= 4.6000 && abs(eta) < 4.7000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 4.6000 && abs(eta) < 4.7000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.088768) + \
    (abs(eta) >= 4.6000 && abs(eta) < 4.7000) * (energy >= 2.0000 && energy < 5.0000) * (0.088768 + (energy-2.000000)* -0.018880) + \
    (abs(eta) >= 4.6000 && abs(eta) < 4.7000) * (energy >= 5.0000 && energy < 10.0000) * (0.032126 + (energy-5.000000)* -0.003627) + \
    (abs(eta) >= 4.6000 && abs(eta) < 4.7000) * (energy >= 10.0000 && energy < 100.0000) * (0.013990 + (energy-10.000000)* 0.000005) + \
    (abs(eta) >= 4.6000 && abs(eta) < 4.7000) * (energy >= 100.0000 && energy < 1000.0000) * (0.014468 + (energy-100.000000)* 0.000005) + \
    (abs(eta) >= 4.6000 && abs(eta) < 4.7000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.019295 + (energy-1000.000000)* 0.000011) + \
    (abs(eta) >= 4.6000 && abs(eta) < 4.7000) * (energy >= 10000.0000) * (0.117161*energy/10000.000000) + \
    (abs(eta) >= 4.7000 && abs(eta) < 4.8000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 4.7000 && abs(eta) < 4.8000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.099752) + \
    (abs(eta) >= 4.7000 && abs(eta) < 4.8000) * (energy >= 2.0000 && energy < 5.0000) * (0.099752 + (energy-2.000000)* -0.021799) + \
    (abs(eta) >= 4.7000 && abs(eta) < 4.8000) * (energy >= 5.0000 && energy < 10.0000) * (0.034356 + (energy-5.000000)* -0.003793) + \
    (abs(eta) >= 4.7000 && abs(eta) < 4.8000) * (energy >= 10.0000 && energy < 100.0000) * (0.015389 + (energy-10.000000)* 0.000006) + \
    (abs(eta) >= 4.7000 && abs(eta) < 4.8000) * (energy >= 100.0000 && energy < 1000.0000) * (0.015906 + (energy-100.000000)* 0.000006) + \
    (abs(eta) >= 4.7000 && abs(eta) < 4.8000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.021315 + (energy-1000.000000)* 0.000012) + \
    (abs(eta) >= 4.7000 && abs(eta) < 4.8000) * (energy >= 10000.0000) * (0.133463*energy/10000.000000) + \
    (abs(eta) >= 4.8000 && abs(eta) < 4.9000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 4.8000 && abs(eta) < 4.9000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.129119) + \
    (abs(eta) >= 4.8000 && abs(eta) < 4.9000) * (energy >= 2.0000 && energy < 5.0000) * (0.129119 + (energy-2.000000)* -0.029895) + \
    (abs(eta) >= 4.8000 && abs(eta) < 4.9000) * (energy >= 5.0000 && energy < 10.0000) * (0.039433 + (energy-5.000000)* -0.004396) + \
    (abs(eta) >= 4.8000 && abs(eta) < 4.9000) * (energy >= 10.0000 && energy < 100.0000) * (0.017452 + (energy-10.000000)* 0.000006) + \
    (abs(eta) >= 4.8000 && abs(eta) < 4.9000) * (energy >= 100.0000 && energy < 1000.0000) * (0.018028 + (energy-100.000000)* 0.000007) + \
    (abs(eta) >= 4.8000 && abs(eta) < 4.9000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.024225 + (energy-1000.000000)* 0.000014) + \
    (abs(eta) >= 4.8000 && abs(eta) < 4.9000) * (energy >= 10000.0000) * (0.151204*energy/10000.000000) + \
    (abs(eta) >= 4.9000 && abs(eta) < 5.0000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 4.9000 && abs(eta) < 5.0000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.148587) + \
    (abs(eta) >= 4.9000 && abs(eta) < 5.0000) * (energy >= 2.0000 && energy < 5.0000) * (0.148587 + (energy-2.000000)* -0.035133) + \
    (abs(eta) >= 4.9000 && abs(eta) < 5.0000) * (energy >= 5.0000 && energy < 10.0000) * (0.043189 + (energy-5.000000)* -0.004829) + \
    (abs(eta) >= 4.9000 && abs(eta) < 5.0000) * (energy >= 10.0000 && energy < 100.0000) * (0.019045 + (energy-10.000000)* 0.000006) + \
    (abs(eta) >= 4.9000 && abs(eta) < 5.0000) * (energy >= 100.0000 && energy < 1000.0000) * (0.019625 + (energy-100.000000)* 0.000007) + \
    (abs(eta) >= 4.9000 && abs(eta) < 5.0000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.026335 + (energy-1000.000000)* 0.000015) + \
    (abs(eta) >= 4.9000 && abs(eta) < 5.0000) * (energy >= 10000.0000) * (0.165614*energy/10000.000000) + \
    (abs(eta) >= 5.0000 && abs(eta) < 5.1000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 5.0000 && abs(eta) < 5.1000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.230194) + \
    (abs(eta) >= 5.0000 && abs(eta) < 5.1000) * (energy >= 2.0000 && energy < 5.0000) * (0.230194 + (energy-2.000000)* -0.059949) + \
    (abs(eta) >= 5.0000 && abs(eta) < 5.1000) * (energy >= 5.0000 && energy < 10.0000) * (0.050347 + (energy-5.000000)* -0.005793) + \
    (abs(eta) >= 5.0000 && abs(eta) < 5.1000) * (energy >= 10.0000 && energy < 100.0000) * (0.021380 + (energy-10.000000)* 0.000007) + \
    (abs(eta) >= 5.0000 && abs(eta) < 5.1000) * (energy >= 100.0000 && energy < 1000.0000) * (0.022044 + (energy-100.000000)* 0.000009) + \
    (abs(eta) >= 5.0000 && abs(eta) < 5.1000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.029961 + (energy-1000.000000)* 0.000018) + \
    (abs(eta) >= 5.0000 && abs(eta) < 5.1000) * (energy >= 10000.0000) * (0.188507*energy/10000.000000) + \
    (abs(eta) >= 5.1000 && abs(eta) < 5.2000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 5.1000 && abs(eta) < 5.2000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.256226) + \
    (abs(eta) >= 5.1000 && abs(eta) < 5.2000) * (energy >= 2.0000 && energy < 5.0000) * (0.256226 + (energy-2.000000)* -0.067581) + \
    (abs(eta) >= 5.1000 && abs(eta) < 5.2000) * (energy >= 5.0000 && energy < 10.0000) * (0.053485 + (energy-5.000000)* -0.006077) + \
    (abs(eta) >= 5.1000 && abs(eta) < 5.2000) * (energy >= 10.0000 && energy < 100.0000) * (0.023102 + (energy-10.000000)* 0.000008) + \
    (abs(eta) >= 5.1000 && abs(eta) < 5.2000) * (energy >= 100.0000 && energy < 1000.0000) * (0.023822 + (energy-100.000000)* 0.000010) + \
    (abs(eta) >= 5.1000 && abs(eta) < 5.2000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.032699 + (energy-1000.000000)* 0.000020) + \
    (abs(eta) >= 5.1000 && abs(eta) < 5.2000) * (energy >= 10000.0000) * (0.212432*energy/10000.000000) + \
    (abs(eta) >= 5.2000 && abs(eta) < 5.3000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 5.2000 && abs(eta) < 5.3000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.428061) + \
    (abs(eta) >= 5.2000 && abs(eta) < 5.3000) * (energy >= 2.0000 && energy < 5.0000) * (0.428061 + (energy-2.000000)* -0.122293) + \
    (abs(eta) >= 5.2000 && abs(eta) < 5.3000) * (energy >= 5.0000 && energy < 10.0000) * (0.061183 + (energy-5.000000)* -0.007029) + \
    (abs(eta) >= 5.2000 && abs(eta) < 5.3000) * (energy >= 10.0000 && energy < 100.0000) * (0.026038 + (energy-10.000000)* 0.000009) + \
    (abs(eta) >= 5.2000 && abs(eta) < 5.3000) * (energy >= 100.0000 && energy < 1000.0000) * (0.026834 + (energy-100.000000)* 0.000011) + \
    (abs(eta) >= 5.2000 && abs(eta) < 5.3000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.037168 + (energy-1000.000000)* 0.000023) + \
    (abs(eta) >= 5.2000 && abs(eta) < 5.3000) * (energy >= 10000.0000) * (0.240865*energy/10000.000000) + \
    (abs(eta) >= 5.3000 && abs(eta) < 5.4000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 5.3000 && abs(eta) < 5.4000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.302517) + \
    (abs(eta) >= 5.3000 && abs(eta) < 5.4000) * (energy >= 2.0000 && energy < 5.0000) * (0.302517 + (energy-2.000000)* -0.078082) + \
    (abs(eta) >= 5.3000 && abs(eta) < 5.4000) * (energy >= 5.0000 && energy < 10.0000) * (0.068271 + (energy-5.000000)* -0.007915) + \
    (abs(eta) >= 5.3000 && abs(eta) < 5.4000) * (energy >= 10.0000 && energy < 100.0000) * (0.028698 + (energy-10.000000)* 0.000010) + \
    (abs(eta) >= 5.3000 && abs(eta) < 5.4000) * (energy >= 100.0000 && energy < 1000.0000) * (0.029589 + (energy-100.000000)* 0.000014) + \
    (abs(eta) >= 5.3000 && abs(eta) < 5.4000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.042319 + (energy-1000.000000)* 0.000027) + \
    (abs(eta) >= 5.3000 && abs(eta) < 5.4000) * (energy >= 10000.0000) * (0.284210*energy/10000.000000) + \
    (abs(eta) >= 5.4000 && abs(eta) < 5.5000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 5.4000 && abs(eta) < 5.5000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.000000) + \
    (abs(eta) >= 5.4000 && abs(eta) < 5.5000) * (energy >= 2.0000 && energy < 5.0000) * (0.000000 + (energy-2.000000)* 0.025424) + \
    (abs(eta) >= 5.4000 && abs(eta) < 5.5000) * (energy >= 5.0000 && energy < 10.0000) * (0.076273 + (energy-5.000000)* -0.008946) + \
    (abs(eta) >= 5.4000 && abs(eta) < 5.5000) * (energy >= 10.0000 && energy < 100.0000) * (0.031543 + (energy-10.000000)* 0.000010) + \
    (abs(eta) >= 5.4000 && abs(eta) < 5.5000) * (energy >= 100.0000 && energy < 1000.0000) * (0.032442 + (energy-100.000000)* 0.000017) + \
    (abs(eta) >= 5.4000 && abs(eta) < 5.5000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.047712 + (energy-1000.000000)* 0.000031) + \
    (abs(eta) >= 5.4000 && abs(eta) < 5.5000) * (energy >= 10000.0000) * (0.326399*energy/10000.000000) + \
    (abs(eta) >= 5.5000 && abs(eta) < 5.6000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 5.5000 && abs(eta) < 5.6000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.000000) + \
    (abs(eta) >= 5.5000 && abs(eta) < 5.6000) * (energy >= 2.0000 && energy < 5.0000) * (0.000000 + (energy-2.000000)* 0.029080) + \
    (abs(eta) >= 5.5000 && abs(eta) < 5.6000) * (energy >= 5.0000 && energy < 10.0000) * (0.087241 + (energy-5.000000)* -0.010397) + \
    (abs(eta) >= 5.5000 && abs(eta) < 5.6000) * (energy >= 10.0000 && energy < 100.0000) * (0.035256 + (energy-10.000000)* 0.000011) + \
    (abs(eta) >= 5.5000 && abs(eta) < 5.6000) * (energy >= 100.0000 && energy < 1000.0000) * (0.036252 + (energy-100.000000)* 0.000019) + \
    (abs(eta) >= 5.5000 && abs(eta) < 5.6000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.053456 + (energy-1000.000000)* 0.000034) + \
    (abs(eta) >= 5.5000 && abs(eta) < 5.6000) * (energy >= 10000.0000) * (0.357698*energy/10000.000000) + \
    (abs(eta) >= 5.6000 && abs(eta) < 5.7000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 5.6000 && abs(eta) < 5.7000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.000000) + \
    (abs(eta) >= 5.6000 && abs(eta) < 5.7000) * (energy >= 2.0000 && energy < 5.0000) * (0.000000 + (energy-2.000000)* 0.033344) + \
    (abs(eta) >= 5.6000 && abs(eta) < 5.7000) * (energy >= 5.0000 && energy < 10.0000) * (0.100031 + (energy-5.000000)* -0.012242) + \
    (abs(eta) >= 5.6000 && abs(eta) < 5.7000) * (energy >= 10.0000 && energy < 100.0000) * (0.038819 + (energy-10.000000)* 0.000015) + \
    (abs(eta) >= 5.6000 && abs(eta) < 5.7000) * (energy >= 100.0000 && energy < 1000.0000) * (0.040162 + (energy-100.000000)* 0.000023) + \
    (abs(eta) >= 5.6000 && abs(eta) < 5.7000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.061269 + (energy-1000.000000)* 0.000040) + \
    (abs(eta) >= 5.6000 && abs(eta) < 5.7000) * (energy >= 10000.0000) * (0.420216*energy/10000.000000) + \
    (abs(eta) >= 5.7000 && abs(eta) < 5.8000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 5.7000 && abs(eta) < 5.8000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.000000) + \
    (abs(eta) >= 5.7000 && abs(eta) < 5.8000) * (energy >= 2.0000 && energy < 5.0000) * (0.000000 + (energy-2.000000)* 0.039257) + \
    (abs(eta) >= 5.7000 && abs(eta) < 5.8000) * (energy >= 5.0000 && energy < 10.0000) * (0.117770 + (energy-5.000000)* -0.014993) + \
    (abs(eta) >= 5.7000 && abs(eta) < 5.8000) * (energy >= 10.0000 && energy < 100.0000) * (0.042808 + (energy-10.000000)* 0.000025) + \
    (abs(eta) >= 5.7000 && abs(eta) < 5.8000) * (energy >= 100.0000 && energy < 1000.0000) * (0.045061 + (energy-100.000000)* 0.000029) + \
    (abs(eta) >= 5.7000 && abs(eta) < 5.8000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.071538 + (energy-1000.000000)* 0.000049) + \
    (abs(eta) >= 5.7000 && abs(eta) < 5.8000) * (energy >= 10000.0000) * (0.510223*energy/10000.000000) + \
    (abs(eta) >= 5.8000 && abs(eta) < 5.9000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 5.8000 && abs(eta) < 5.9000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.000000) + \
    (abs(eta) >= 5.8000 && abs(eta) < 5.9000) * (energy >= 2.0000 && energy < 5.0000) * (0.000000 + (energy-2.000000)* 0.041967) + \
    (abs(eta) >= 5.8000 && abs(eta) < 5.9000) * (energy >= 5.0000 && energy < 10.0000) * (0.125902 + (energy-5.000000)* -0.016273) + \
    (abs(eta) >= 5.8000 && abs(eta) < 5.9000) * (energy >= 10.0000 && energy < 100.0000) * (0.044536 + (energy-10.000000)* 0.000023) + \
    (abs(eta) >= 5.8000 && abs(eta) < 5.9000) * (energy >= 100.0000 && energy < 1000.0000) * (0.046571 + (energy-100.000000)* 0.000043) + \
    (abs(eta) >= 5.8000 && abs(eta) < 5.9000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.085014 + (energy-1000.000000)* 0.000066) + \
    (abs(eta) >= 5.8000 && abs(eta) < 5.9000) * (energy >= 10000.0000) * (0.680639*energy/10000.000000) + \
    (abs(eta) >= 5.9000 && abs(eta) < 6.0000) * (energy >= 0.0000 && energy < 1.0000) * (0.00000000) + \
    (abs(eta) >= 5.9000 && abs(eta) < 6.0000) * (energy >= 1.0000 && energy < 2.0000) * (0.000000 + (energy-1.000000)* 0.000000) + \
    (abs(eta) >= 5.9000 && abs(eta) < 6.0000) * (energy >= 2.0000 && energy < 5.0000) * (0.000000 + (energy-2.000000)* 0.000000) + \
    (abs(eta) >= 5.9000 && abs(eta) < 6.0000) * (energy >= 5.0000 && energy < 10.0000) * (0.000000 + (energy-5.000000)* 0.013098) + \
    (abs(eta) >= 5.9000 && abs(eta) < 6.0000) * (energy >= 10.0000 && energy < 100.0000) * (0.065492 + (energy-10.000000)* 0.000107) + \
    (abs(eta) >= 5.9000 && abs(eta) < 6.0000) * (energy >= 100.0000 && energy < 1000.0000) * (0.075097 + (energy-100.000000)* 0.000243) + \
    (abs(eta) >= 5.9000 && abs(eta) < 6.0000) * (energy >= 1000.0000 && energy < 10000.0000) * (0.294198 + (energy-1000.000000)* 0.000280) + \
    (abs(eta) >= 5.9000 && abs(eta) < 6.0000) * (energy >= 10000.0000) * (2.814894*energy/10000.000000)
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
    set ResolutionFormula {
	(abs(eta) <= 0.78 )  * sqrt(energy^2*0.01^2 + energy*0.156^2)+
	(abs(eta) > 0.78 && abs(eta) <=0.83 ) * sqrt( energy^0.01^2 + energy*0.175^2  ) +
	(abs(eta) <= 2.5 && abs(eta) > 0.83) * sqrt( energy^2*0.01^2 + energy*0.151^2  )}
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
    #HCAL endcaps: dphi = 6 degree, deta = 0.1  up to |eta| <=2.5
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
    # deta =0.1 for 0.9 < |eta| <=2.5
    #for -2.5 to -0.9, 21 segments
    for {set i 1} {$i <=17} {incr i} {
	set eta [expr {-2.5 + $i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }
    #same for 0.9 to 2.5
    for {set i 1} {$i <=17} {incr i} {
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
	(abs(eta)<= 0.3) * sqrt(1.38^2  + energy*0.308^2  + energy^2*0.050^2) +
	 (abs(eta)<= 0.78 && abs(eta) > 0.3) * sqrt(1.25^2  + energy*0.322^2  + energy^2*0.048^2) +
	 (abs(eta)<=1.099 && abs(eta) > 0.78) * sqrt(  1.159^2 + energy*0.341^2 + energy^2*0.049^2 ) +
	 (abs(eta)<=2.5 && abs(eta)> 1.099) * sqrt( 1.09^2 + energy*0.319^2 + energy^2*0.052^2  )
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
    add InputArray MuonMomentumSmearing/muons
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
    set EfficiencyFormula {
	(energy < 2.0 ) * (0.000) +
	(energy >= 2.0) * (abs(eta) < 0.7)*(0.94) +
	(energy >= 2.0) * (abs(eta) >=0.7 && abs(eta) <=2.5) * (0.9)	}

}


##################
# Photon isolation
##################

module Isolation PhotonIsolation {
    set CandidateInputArray PhotonEfficiency/photons
    set IsolationInputArray EFlowMerger/eflow

    set OutputArray photons

    set DeltaRMax 0.1

    set PTMin 0.5

    set PTRatioMax 0.2
}


#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
    set InputArray ElectronFilter/electrons
    set OutputArray electrons

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    set EfficiencyFormula {
	(energy < 3.0 ) * ( 0.00 ) +
  (abs(eta) > 2.50) * ( 0.00 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (0.58 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * ( 0.7 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * ( 0.6 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.7 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.8 ) +
	( energy >=3 && energy < 8  ) * (abs(eta) <= 0.69)                    * (0.84 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (  0.6 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * ( 0.76 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * ( 0.67 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.78 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.86 ) +
	( energy >=8 && energy < 13  ) * (abs(eta) <= 0.69)                    * ( 0.88 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) > 1.95 && abs(eta) < 2.50) * (  0.6 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (  0.8 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * ( 0.68 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.84 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.88 ) +
	( energy >=13 && energy < 18  ) * (abs(eta) <= 0.69)                    * (  0.9 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (0.64 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (0.82 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * ( 0.7 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * (0.84 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.9 ) +
	( energy >=18 && energy < 23  ) * (abs(eta) <= 0.69)                    * (0.92 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) > 1.95 && abs(eta) < 2.50) * (0.64 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (0.86 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (0.74 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * (0.87 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * (0.91 ) +
	( energy >= 23 && energy < 28  ) * (abs(eta) <= 0.69)                    * (0.94 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)   * (0.67 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (0.88 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (0.78 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.9 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * (0.94 ) +
	( energy >=28 && energy < 35  ) * (abs(eta) <= 0.69)                    * (0.94 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (0.68 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * ( 0.9 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (0.86 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * (0.92 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * (0.94 ) +
	( energy >=35 && energy < 45  ) * (abs(eta) <= 0.69)                    * (0.96 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (  0.7 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * ( 0.92 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (  0.8 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.94 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.96 ) +
	( energy >=45 && energy < 80  ) * (abs(eta) <= 0.69)                    * ( 0.97 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) > 1.95 && abs(eta) < 2.50) * (0.68 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (0.96 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (0.84 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * (0.94 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * (0.98 ) +
	( energy >=80 && energy < 200  ) * (abs(eta) <= 0.69)                    * (0.98 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) > 1.95 && abs(eta) < 2.50) * ( 0.68 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * ( 0.97 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * ( 0.86 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * ( 0.96 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * ( 0.98 ) +
	( energy >=200 && energy < 400  ) * (abs(eta) <= 0.69)                    * ( 0.98 ) +
	( energy >=400  ) * (abs(eta) > 1.95 && abs(eta) < 2.50)  * (0.68 ) +
	( energy >=400  ) * (abs(eta) <= 1.95 && abs(eta) > 1.22) * (0.96 ) +
	( energy >=400  ) * (abs(eta) <= 1.22 && abs(eta) > 1.1 ) * (0.82 ) +
	( energy >=400  ) * (abs(eta) <= 1.1 && abs(eta) > 0.91 ) * (0.96 ) +
	( energy >=400  ) * (abs(eta) <= 0.91 && abs(eta) > 0.69) * (0.98 ) +
	( energy >=400  ) * (abs(eta) <= 0.69)                    * (0.98 )
    }
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
    set CandidateInputArray ElectronEfficiency/electrons
    set IsolationInputArray EFlowMerger/eflow

    set OutputArray electrons

    set DeltaRMax 0.1

    set PTMin 0.5

    set PTRatioMax 0.2
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
    set InputArray MuonMomentumSmearing/muons
    set OutputArray muons

    # set EfficiencyFormula {efficiency as a function of eta and pt}

    set EfficiencyFormula {
	(energy < 2.5 )     * (0.00) +
	(energy>=2.5  )     * (0.999)
    }
}

################
# Muon isolation
################

module Isolation MuonIsolation {
    set CandidateInputArray MuonEfficiency/muons
    set IsolationInputArray EFlowMerger/eflow

    set OutputArray muons

    set DeltaRMax 0.1

    set PTMin 0.5

    set PTRatioMax 0.2
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
module UniqueObjectFinder UniqueObjectFinder {
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


source MuonCollider/MuonColliderDet_JetReco.tcl


#########################################
# Jet Momentum Smearing to mimick overlay
#########################################


source MuonCollider/MuonColliderDet_JetSmearing.tcl



########################
# Jet Flavor Association
########################

source  MuonCollider/MuonColliderDet_JetFlavorAssociation.tcl

###########
# b-tagging
###########
# based on CLICdp-Note-2014-002

source  MuonCollider/MuonColliderDet_BTagging.tcl


#############
# tau-tagging
#############

source MuonCollider/MuonColliderDet_TauTagging.tcl

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
    # add Branch InputArray BranchName BranchClass
    add Branch Delphes/allParticles Particle GenParticle

    add Branch GenJetFinder/jets GenJet Jet

    add Branch FastJetFinderKt/KTjets KTjet Jet

    add Branch JetMomentumSmearing_VLCR05N2/JER_VLCjetsR05N2 VLCjetR05N2 Jet
    add Branch JetMomentumSmearing_VLCR05N3/JER_VLCjetsR05N3 VLCjetR05N3 Jet
    add Branch JetMomentumSmearing_VLCR05N4/JER_VLCjetsR05N4 VLCjetR05N4 Jet
    add Branch JetMomentumSmearing_VLCR05N5/JER_VLCjetsR05N5 VLCjetR05N5 Jet
    add Branch JetMomentumSmearing_VLCR05N6/JER_VLCjetsR05N6 VLCjetR05N6 Jet

    add Branch JetMomentumSmearing_VLCR07N2/JER_VLCjetsR07N2 VLCjetR07N2 Jet
    add Branch JetMomentumSmearing_VLCR07N3/JER_VLCjetsR07N3 VLCjetR07N3 Jet
    add Branch JetMomentumSmearing_VLCR07N4/JER_VLCjetsR07N4 VLCjetR07N4 Jet
    add Branch JetMomentumSmearing_VLCR07N5/JER_VLCjetsR07N5 VLCjetR07N5 Jet
    add Branch JetMomentumSmearing_VLCR07N6/JER_VLCjetsR07N6 VLCjetR07N6 Jet

    add Branch JetMomentumSmearing_VLCR10N2/JER_VLCjetsR10N2 VLCjetR10N2 Jet
    add Branch JetMomentumSmearing_VLCR10N3/JER_VLCjetsR10N3 VLCjetR10N3 Jet
    add Branch JetMomentumSmearing_VLCR10N4/JER_VLCjetsR10N4 VLCjetR10N4 Jet
    add Branch JetMomentumSmearing_VLCR10N5/JER_VLCjetsR10N5 VLCjetR10N5 Jet
    add Branch JetMomentumSmearing_VLCR10N6/JER_VLCjetsR10N6 VLCjetR10N6 Jet

    add Branch JetMomentumSmearing_VLCR12N2/JER_VLCjetsR12N2 VLCjetR12N2 Jet
    add Branch JetMomentumSmearing_VLCR12N3/JER_VLCjetsR12N3 VLCjetR12N3 Jet
    add Branch JetMomentumSmearing_VLCR12N4/JER_VLCjetsR12N4 VLCjetR12N4 Jet
    add Branch JetMomentumSmearing_VLCR12N5/JER_VLCjetsR12N5 VLCjetR12N5 Jet
    add Branch JetMomentumSmearing_VLCR12N6/JER_VLCjetsR12N6 VLCjetR12N6 Jet

    add Branch JetMomentumSmearing_VLCR15N2/JER_VLCjetsR15N2 VLCjetR15N2 Jet
    add Branch JetMomentumSmearing_VLCR15N3/JER_VLCjetsR15N3 VLCjetR15N3 Jet
    add Branch JetMomentumSmearing_VLCR15N4/JER_VLCjetsR15N4 VLCjetR15N4 Jet
    add Branch JetMomentumSmearing_VLCR15N5/JER_VLCjetsR15N5 VLCjetR15N5 Jet
    add Branch JetMomentumSmearing_VLCR15N6/JER_VLCjetsR15N6 VLCjetR15N6 Jet

    add Branch JetMomentumSmearing_VLCR02_inclusive/JER_VLCjetsR02_inclusive VLCjetR02_inclusive Jet
    add Branch JetMomentumSmearing_VLCR05_inclusive/JER_VLCjetsR05_inclusive VLCjetR05_inclusive Jet
    add Branch JetMomentumSmearing_VLCR07_inclusive/JER_VLCjetsR07_inclusive VLCjetR07_inclusive Jet
    add Branch JetMomentumSmearing_VLCR10_inclusive/JER_VLCjetsR10_inclusive VLCjetR10_inclusive Jet
    add Branch JetMomentumSmearing_VLCR12_inclusive/JER_VLCjetsR12_inclusive VLCjetR12_inclusive Jet
    add Branch JetMomentumSmearing_VLCR15_inclusive/JER_VLCjetsR15_inclusive VLCjetR15_inclusive Jet


    ####

    add Branch GenMissingET/momentum GenMissingET MissingET

    add Branch TrackMerger/tracks Track Track
    add Branch Calorimeter/towers Tower Tower

    add Branch HCal/eflowTracks EFlowTrack Track
    add Branch ECal/eflowPhotons EFlowPhoton Tower
    add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

    add Branch UniqueObjectFinder/photons Photon Photon
    add Branch UniqueObjectFinder/electrons Electron Electron
    add Branch UniqueObjectFinder/muons Muon Muon
    add Branch ForwardMuonMomentumSmearing/muons ForwardMuon Muon

    add Branch MissingET/momentum MissingET MissingET
    add Branch ScalarHT/energy ScalarHT ScalarHT
}
