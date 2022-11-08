
###########################################
#
#  Main authors: Michele Selvaggi (CERN)
#
#  Released on: May 2018
#
#  Based on CMS_PhaseII_200PU_v03
#
#  Test card for CMS PhaseII detector without pile-up (use for testing only)
#
############################################


#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronMomentumSmearing
  ElectronEnergySmearing
  MuonMomentumSmearing

  TrackMerger

  ECal
  HCal

  PhotonEnergySmearing
  ElectronFilter

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

  MissingET

  NeutrinoFilter
  GenJetFinder
  GenMissingET

  FastJetFinder
  FatJetFinder

  JetEnergyScale

  JetFlavorAssociation

  BTagging
  TauTagging

  UniqueObjectFinder
  GenParticleFilter

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

  # radius of the magnetic field coverage, in m
  set Radius 1.29
  # half-length of the magnetic field coverage, in m
  set HalfLength 3.00

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
          (abs(eta) <= 1.2) * (pt > 0.2 && pt <= 1.0) * (pt * 1.00) + \
          (abs(eta) <= 1.2) * (pt > 1.0) * (1.00) + \
          (abs(eta) > 1.2 && abs(eta) <= 2.8) * (pt > 0.2 && pt <= 1.0) * (pt*1.00) + \
          (abs(eta) > 1.2 && abs(eta) <= 2.8) * (pt > 1.0) * (1.00) + \
          (abs(eta) > 2.8 && abs(eta) <= 4.0) * (pt > 0.2 && pt <= 1.0) * (pt*0.95) + \
          (abs(eta) > 2.8 && abs(eta) <= 4.0) * (pt > 1.0) * (0.95) + \
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
  # resolution formula for charged hadrons ,

  #
  # Automatically generated tracker resolution formula for layout: OT612IT4025
  #
  #  By Unknown author on: 2017-06-30.17:03:00
  #
  set ResolutionFormula {    (abs(eta) >= 0.0000 && abs(eta) < 0.2000) * (pt >= 0.0000 && pt < 1.0000) * (0.00457888) + \
     (abs(eta) >= 0.0000 && abs(eta) < 0.2000) * (pt >= 1.0000 && pt < 10.0000) * (0.004579 + (pt-1.000000)* 0.000045) + \
     (abs(eta) >= 0.0000 && abs(eta) < 0.2000) * (pt >= 10.0000 && pt < 100.0000) * (0.004983 + (pt-10.000000)* 0.000047) + \
     (abs(eta) >= 0.0000 && abs(eta) < 0.2000) * (pt >= 100.0000) * (0.009244*pt/100.000000) + \
     (abs(eta) >= 0.2000 && abs(eta) < 0.4000) * (pt >= 0.0000 && pt < 1.0000) * (0.00505011) + \
     (abs(eta) >= 0.2000 && abs(eta) < 0.4000) * (pt >= 1.0000 && pt < 10.0000) * (0.005050 + (pt-1.000000)* 0.000033) + \
     (abs(eta) >= 0.2000 && abs(eta) < 0.4000) * (pt >= 10.0000 && pt < 100.0000) * (0.005343 + (pt-10.000000)* 0.000043) + \
     (abs(eta) >= 0.2000 && abs(eta) < 0.4000) * (pt >= 100.0000) * (0.009172*pt/100.000000) + \
     (abs(eta) >= 0.4000 && abs(eta) < 0.6000) * (pt >= 0.0000 && pt < 1.0000) * (0.00510573) + \
     (abs(eta) >= 0.4000 && abs(eta) < 0.6000) * (pt >= 1.0000 && pt < 10.0000) * (0.005106 + (pt-1.000000)* 0.000023) + \
     (abs(eta) >= 0.4000 && abs(eta) < 0.6000) * (pt >= 10.0000 && pt < 100.0000) * (0.005317 + (pt-10.000000)* 0.000042) + \
     (abs(eta) >= 0.4000 && abs(eta) < 0.6000) * (pt >= 100.0000) * (0.009077*pt/100.000000) + \
     (abs(eta) >= 0.6000 && abs(eta) < 0.8000) * (pt >= 0.0000 && pt < 1.0000) * (0.00578020) + \
     (abs(eta) >= 0.6000 && abs(eta) < 0.8000) * (pt >= 1.0000 && pt < 10.0000) * (0.005780 + (pt-1.000000)* -0.000000) + \
     (abs(eta) >= 0.6000 && abs(eta) < 0.8000) * (pt >= 10.0000 && pt < 100.0000) * (0.005779 + (pt-10.000000)* 0.000038) + \
     (abs(eta) >= 0.6000 && abs(eta) < 0.8000) * (pt >= 100.0000) * (0.009177*pt/100.000000) + \
     (abs(eta) >= 0.8000 && abs(eta) < 1.0000) * (pt >= 0.0000 && pt < 1.0000) * (0.00728723) + \
     (abs(eta) >= 0.8000 && abs(eta) < 1.0000) * (pt >= 1.0000 && pt < 10.0000) * (0.007287 + (pt-1.000000)* -0.000031) + \
     (abs(eta) >= 0.8000 && abs(eta) < 1.0000) * (pt >= 10.0000 && pt < 100.0000) * (0.007011 + (pt-10.000000)* 0.000038) + \
     (abs(eta) >= 0.8000 && abs(eta) < 1.0000) * (pt >= 100.0000) * (0.010429*pt/100.000000) + \
     (abs(eta) >= 1.0000 && abs(eta) < 1.2000) * (pt >= 0.0000 && pt < 1.0000) * (0.01045117) + \
     (abs(eta) >= 1.0000 && abs(eta) < 1.2000) * (pt >= 1.0000 && pt < 10.0000) * (0.010451 + (pt-1.000000)* -0.000051) + \
     (abs(eta) >= 1.0000 && abs(eta) < 1.2000) * (pt >= 10.0000 && pt < 100.0000) * (0.009989 + (pt-10.000000)* 0.000043) + \
     (abs(eta) >= 1.0000 && abs(eta) < 1.2000) * (pt >= 100.0000) * (0.013867*pt/100.000000) + \
     (abs(eta) >= 1.2000 && abs(eta) < 1.4000) * (pt >= 0.0000 && pt < 1.0000) * (0.01477199) + \
     (abs(eta) >= 1.2000 && abs(eta) < 1.4000) * (pt >= 1.0000 && pt < 10.0000) * (0.014772 + (pt-1.000000)* -0.000128) + \
     (abs(eta) >= 1.2000 && abs(eta) < 1.4000) * (pt >= 10.0000 && pt < 100.0000) * (0.013616 + (pt-10.000000)* 0.000035) + \
     (abs(eta) >= 1.2000 && abs(eta) < 1.4000) * (pt >= 100.0000) * (0.016800*pt/100.000000) + \
     (abs(eta) >= 1.4000 && abs(eta) < 1.6000) * (pt >= 0.0000 && pt < 1.0000) * (0.01731474) + \
     (abs(eta) >= 1.4000 && abs(eta) < 1.6000) * (pt >= 1.0000 && pt < 10.0000) * (0.017315 + (pt-1.000000)* -0.000208) + \
     (abs(eta) >= 1.4000 && abs(eta) < 1.6000) * (pt >= 10.0000 && pt < 100.0000) * (0.015439 + (pt-10.000000)* 0.000030) + \
     (abs(eta) >= 1.4000 && abs(eta) < 1.6000) * (pt >= 100.0000) * (0.018161*pt/100.000000) + \
     (abs(eta) >= 1.6000 && abs(eta) < 1.8000) * (pt >= 0.0000 && pt < 1.0000) * (0.01942025) + \
     (abs(eta) >= 1.6000 && abs(eta) < 1.8000) * (pt >= 1.0000 && pt < 10.0000) * (0.019420 + (pt-1.000000)* -0.000417) + \
     (abs(eta) >= 1.6000 && abs(eta) < 1.8000) * (pt >= 10.0000 && pt < 100.0000) * (0.015669 + (pt-10.000000)* 0.000026) + \
     (abs(eta) >= 1.6000 && abs(eta) < 1.8000) * (pt >= 100.0000) * (0.018039*pt/100.000000) + \
     (abs(eta) >= 1.8000 && abs(eta) < 2.0000) * (pt >= 0.0000 && pt < 1.0000) * (0.02201432) + \
     (abs(eta) >= 1.8000 && abs(eta) < 2.0000) * (pt >= 1.0000 && pt < 10.0000) * (0.022014 + (pt-1.000000)* -0.000667) + \
     (abs(eta) >= 1.8000 && abs(eta) < 2.0000) * (pt >= 10.0000 && pt < 100.0000) * (0.016012 + (pt-10.000000)* 0.000045) + \
     (abs(eta) >= 1.8000 && abs(eta) < 2.0000) * (pt >= 100.0000) * (0.020098*pt/100.000000) + \
     (abs(eta) >= 2.0000 && abs(eta) < 2.2000) * (pt >= 0.0000 && pt < 1.0000) * (0.02574300) + \
     (abs(eta) >= 2.0000 && abs(eta) < 2.2000) * (pt >= 1.0000 && pt < 10.0000) * (0.025743 + (pt-1.000000)* -0.001118) + \
     (abs(eta) >= 2.0000 && abs(eta) < 2.2000) * (pt >= 10.0000 && pt < 100.0000) * (0.015681 + (pt-10.000000)* 0.000051) + \
     (abs(eta) >= 2.0000 && abs(eta) < 2.2000) * (pt >= 100.0000) * (0.020289*pt/100.000000) + \
     (abs(eta) >= 2.2000 && abs(eta) < 2.4000) * (pt >= 0.0000 && pt < 1.0000) * (0.02885821) + \
     (abs(eta) >= 2.2000 && abs(eta) < 2.4000) * (pt >= 1.0000 && pt < 10.0000) * (0.028858 + (pt-1.000000)* -0.001345) + \
     (abs(eta) >= 2.2000 && abs(eta) < 2.4000) * (pt >= 10.0000 && pt < 100.0000) * (0.016753 + (pt-10.000000)* 0.000053) + \
     (abs(eta) >= 2.2000 && abs(eta) < 2.4000) * (pt >= 100.0000) * (0.021524*pt/100.000000) + \
     (abs(eta) >= 2.4000 && abs(eta) < 2.6000) * (pt >= 0.0000 && pt < 1.0000) * (0.03204812) + \
     (abs(eta) >= 2.4000 && abs(eta) < 2.6000) * (pt >= 1.0000 && pt < 10.0000) * (0.032048 + (pt-1.000000)* -0.001212) + \
     (abs(eta) >= 2.4000 && abs(eta) < 2.6000) * (pt >= 10.0000 && pt < 100.0000) * (0.021138 + (pt-10.000000)* 0.000037) + \
     (abs(eta) >= 2.4000 && abs(eta) < 2.6000) * (pt >= 100.0000) * (0.024477*pt/100.000000) + \
     (abs(eta) >= 2.6000 && abs(eta) < 2.8000) * (pt >= 0.0000 && pt < 1.0000) * (0.03950405) + \
     (abs(eta) >= 2.6000 && abs(eta) < 2.8000) * (pt >= 1.0000 && pt < 10.0000) * (0.039504 + (pt-1.000000)* -0.001386) + \
     (abs(eta) >= 2.6000 && abs(eta) < 2.8000) * (pt >= 10.0000 && pt < 100.0000) * (0.027026 + (pt-10.000000)* 0.000037) + \
     (abs(eta) >= 2.6000 && abs(eta) < 2.8000) * (pt >= 100.0000) * (0.030392*pt/100.000000) + \
     (abs(eta) >= 2.8000 && abs(eta) < 3.0000) * (pt >= 0.0000 && pt < 1.0000) * (0.04084751) + \
     (abs(eta) >= 2.8000 && abs(eta) < 3.0000) * (pt >= 1.0000 && pt < 10.0000) * (0.040848 + (pt-1.000000)* -0.001780) + \
     (abs(eta) >= 2.8000 && abs(eta) < 3.0000) * (pt >= 10.0000 && pt < 100.0000) * (0.024824 + (pt-10.000000)* 0.000029) + \
     (abs(eta) >= 2.8000 && abs(eta) < 3.0000) * (pt >= 100.0000) * (0.027445*pt/100.000000) + \
     (abs(eta) >= 3.0000 && abs(eta) < 3.2000) * (pt >= 0.0000 && pt < 1.0000) * (0.04532425) + \
     (abs(eta) >= 3.0000 && abs(eta) < 3.2000) * (pt >= 1.0000 && pt < 10.0000) * (0.045324 + (pt-1.000000)* -0.002497) + \
     (abs(eta) >= 3.0000 && abs(eta) < 3.2000) * (pt >= 10.0000 && pt < 100.0000) * (0.022851 + (pt-10.000000)* 0.000024) + \
     (abs(eta) >= 3.0000 && abs(eta) < 3.2000) * (pt >= 100.0000) * (0.025053*pt/100.000000) + \
     (abs(eta) >= 3.2000 && abs(eta) < 3.4000) * (pt >= 0.0000 && pt < 1.0000) * (0.06418925) + \
     (abs(eta) >= 3.2000 && abs(eta) < 3.4000) * (pt >= 1.0000 && pt < 10.0000) * (0.064189 + (pt-1.000000)* -0.004055) + \
     (abs(eta) >= 3.2000 && abs(eta) < 3.4000) * (pt >= 10.0000 && pt < 100.0000) * (0.027691 + (pt-10.000000)* 0.000034) + \
     (abs(eta) >= 3.2000 && abs(eta) < 3.4000) * (pt >= 100.0000) * (0.030710*pt/100.000000) + \
     (abs(eta) >= 3.4000 && abs(eta) < 3.6000) * (pt >= 0.0000 && pt < 1.0000) * (0.07682500) + \
     (abs(eta) >= 3.4000 && abs(eta) < 3.6000) * (pt >= 1.0000 && pt < 10.0000) * (0.076825 + (pt-1.000000)* -0.004510) + \
     (abs(eta) >= 3.4000 && abs(eta) < 3.6000) * (pt >= 10.0000 && pt < 100.0000) * (0.036234 + (pt-10.000000)* 0.000049) + \
     (abs(eta) >= 3.4000 && abs(eta) < 3.6000) * (pt >= 100.0000) * (0.040629*pt/100.000000) + \
     (abs(eta) >= 3.6000 && abs(eta) < 3.8000) * (pt >= 0.0000 && pt < 1.0000) * (0.09796358) + \
     (abs(eta) >= 3.6000 && abs(eta) < 3.8000) * (pt >= 1.0000 && pt < 10.0000) * (0.097964 + (pt-1.000000)* -0.005758) + \
     (abs(eta) >= 3.6000 && abs(eta) < 3.8000) * (pt >= 10.0000 && pt < 100.0000) * (0.046145 + (pt-10.000000)* 0.000069) + \
     (abs(eta) >= 3.6000 && abs(eta) < 3.8000) * (pt >= 100.0000) * (0.052345*pt/100.000000) + \
     (abs(eta) >= 3.8000 && abs(eta) < 4.0000) * (pt >= 0.0000 && pt < 1.0000) * (0.13415929) + \
     (abs(eta) >= 3.8000 && abs(eta) < 4.0000) * (pt >= 1.0000 && pt < 10.0000) * (0.134159 + (pt-1.000000)* -0.008283) + \
     (abs(eta) >= 3.8000 && abs(eta) < 4.0000) * (pt >= 10.0000 && pt < 100.0000) * (0.059612 + (pt-10.000000)* 0.000111) + \
     (abs(eta) >= 3.8000 && abs(eta) < 4.0000) * (pt >= 100.0000) * (0.069617*pt/100.000000)
  }


}


#################################
# Energy resolution for electrons
#################################

module EnergySmearing ElectronEnergySmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set OutputArray electrons

  # set ResolutionFormula {resolution formula as a function of eta and energy}

  # resolution formula for electrons

  # taking something flat in energy for now, ECAL will take over at high energy anyway.
  # inferred from hep-ex/1306.2016 and 1502.02701
  set ResolutionFormula {

                        (abs(eta) <= 1.5)  * (energy*0.028) +
    (abs(eta) > 1.5  && abs(eta) <= 1.75)  * (energy*0.037) +
    (abs(eta) > 1.75  && abs(eta) <= 2.15) * (energy*0.038) +
    (abs(eta) > 2.15  && abs(eta) <= 3.00) * (energy*0.044) +
    (abs(eta) > 3.00  && abs(eta) <= 4.00) * (energy*0.10)}

}

###############################
# Momentum resolution for muons
###############################

module MomentumSmearing MuonMomentumSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set OutputArray muons
  # resolution formula for muons

  # up to |eta| < 2.8 take measurement from tracking + muon chambers
  # for |eta| > 2.8 and pT < 5.0 take measurement from tracking alone taken from
  # http://mersi.web.cern.ch/mersi/layouts/.private/Baseline_tilted_200_Pixel_1_1_1/index.html
  source muonMomentumResolution.tcl
}


##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronMomentumSmearing/chargedHadrons
  add InputArray ElectronEnergySmearing/electrons
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
          for {set i -139} {$i <= 139} {incr i} {
              add PhiBins [expr {$i * $pi/139}]
          }
  add EtaPhiBins -4.18 $PhiBins
  add EtaPhiBins 4.21 $PhiBins
  set PhiBins {}
          for {set i -142} {$i <= 142} {incr i} {
              add PhiBins [expr {$i * $pi/142}]
          }
  add EtaPhiBins -4.16 $PhiBins
  add EtaPhiBins 4.18 $PhiBins
  set PhiBins {}
          for {set i -145} {$i <= 145} {incr i} {
              add PhiBins [expr {$i * $pi/145}]
          }
  add EtaPhiBins -4.14 $PhiBins
  add EtaPhiBins 4.16 $PhiBins
  set PhiBins {}
          for {set i -148} {$i <= 148} {incr i} {
              add PhiBins [expr {$i * $pi/148}]
          }
  add EtaPhiBins -4.12 $PhiBins
  add EtaPhiBins 4.14 $PhiBins
  set PhiBins {}
          for {set i -151} {$i <= 151} {incr i} {
              add PhiBins [expr {$i * $pi/151}]
          }
  add EtaPhiBins -4.10 $PhiBins
  add EtaPhiBins 4.12 $PhiBins
  set PhiBins {}
          for {set i -155} {$i <= 155} {incr i} {
              add PhiBins [expr {$i * $pi/155}]
          }
  add EtaPhiBins -4.08 $PhiBins
  add EtaPhiBins 4.10 $PhiBins
  set PhiBins {}
          for {set i -158} {$i <= 158} {incr i} {
              add PhiBins [expr {$i * $pi/158}]
          }
  add EtaPhiBins -4.06 $PhiBins
  add EtaPhiBins 4.08 $PhiBins
  set PhiBins {}
          for {set i -161} {$i <= 161} {incr i} {
              add PhiBins [expr {$i * $pi/161}]
          }
  add EtaPhiBins -4.04 $PhiBins
  add EtaPhiBins 4.06 $PhiBins
  set PhiBins {}
          for {set i -164} {$i <= 164} {incr i} {
              add PhiBins [expr {$i * $pi/164}]
          }
  add EtaPhiBins -4.02 $PhiBins
  add EtaPhiBins 4.04 $PhiBins
  set PhiBins {}
          for {set i -167} {$i <= 167} {incr i} {
              add PhiBins [expr {$i * $pi/167}]
          }
  add EtaPhiBins -4.00 $PhiBins
  add EtaPhiBins 4.02 $PhiBins
  set PhiBins {}
          for {set i -170} {$i <= 170} {incr i} {
              add PhiBins [expr {$i * $pi/170}]
          }
  add EtaPhiBins -3.98 $PhiBins
  add EtaPhiBins 4.00 $PhiBins
  set PhiBins {}
          for {set i -173} {$i <= 173} {incr i} {
              add PhiBins [expr {$i * $pi/173}]
          }
  add EtaPhiBins -3.96 $PhiBins
  add EtaPhiBins 3.98 $PhiBins
  set PhiBins {}
          for {set i -176} {$i <= 176} {incr i} {
              add PhiBins [expr {$i * $pi/176}]
          }
  add EtaPhiBins -3.95 $PhiBins
  add EtaPhiBins 3.96 $PhiBins
  set PhiBins {}
          for {set i -180} {$i <= 180} {incr i} {
              add PhiBins [expr {$i * $pi/180}]
          }
  add EtaPhiBins -3.93 $PhiBins
  add EtaPhiBins 3.95 $PhiBins
  set PhiBins {}
          for {set i -183} {$i <= 183} {incr i} {
              add PhiBins [expr {$i * $pi/183}]
          }
  add EtaPhiBins -3.91 $PhiBins
  add EtaPhiBins 3.93 $PhiBins
  set PhiBins {}
          for {set i -186} {$i <= 186} {incr i} {
              add PhiBins [expr {$i * $pi/186}]
          }
  add EtaPhiBins -3.90 $PhiBins
  add EtaPhiBins 3.91 $PhiBins
  set PhiBins {}
          for {set i -189} {$i <= 189} {incr i} {
              add PhiBins [expr {$i * $pi/189}]
          }
  add EtaPhiBins -3.88 $PhiBins
  add EtaPhiBins 3.90 $PhiBins
  set PhiBins {}
          for {set i -192} {$i <= 192} {incr i} {
              add PhiBins [expr {$i * $pi/192}]
          }
  add EtaPhiBins -3.86 $PhiBins
  add EtaPhiBins 3.88 $PhiBins
  set PhiBins {}
          for {set i -195} {$i <= 195} {incr i} {
              add PhiBins [expr {$i * $pi/195}]
          }
  add EtaPhiBins -3.85 $PhiBins
  add EtaPhiBins 3.86 $PhiBins
  set PhiBins {}
          for {set i -198} {$i <= 198} {incr i} {
              add PhiBins [expr {$i * $pi/198}]
          }
  add EtaPhiBins -3.83 $PhiBins
  add EtaPhiBins 3.85 $PhiBins
  set PhiBins {}
          for {set i -202} {$i <= 202} {incr i} {
              add PhiBins [expr {$i * $pi/202}]
          }
  add EtaPhiBins -3.82 $PhiBins
  add EtaPhiBins 3.83 $PhiBins
  set PhiBins {}
          for {set i -205} {$i <= 205} {incr i} {
              add PhiBins [expr {$i * $pi/205}]
          }
  add EtaPhiBins -3.80 $PhiBins
  add EtaPhiBins 3.82 $PhiBins
  set PhiBins {}
          for {set i -208} {$i <= 208} {incr i} {
              add PhiBins [expr {$i * $pi/208}]
          }
  add EtaPhiBins -3.78 $PhiBins
  add EtaPhiBins 3.80 $PhiBins
  set PhiBins {}
          for {set i -211} {$i <= 211} {incr i} {
              add PhiBins [expr {$i * $pi/211}]
          }
  add EtaPhiBins -3.77 $PhiBins
  add EtaPhiBins 3.78 $PhiBins
  set PhiBins {}
          for {set i -214} {$i <= 214} {incr i} {
              add PhiBins [expr {$i * $pi/214}]
          }
  add EtaPhiBins -3.76 $PhiBins
  add EtaPhiBins 3.77 $PhiBins
  set PhiBins {}
          for {set i -217} {$i <= 217} {incr i} {
              add PhiBins [expr {$i * $pi/217}]
          }
  add EtaPhiBins -3.74 $PhiBins
  add EtaPhiBins 3.76 $PhiBins
  set PhiBins {}
          for {set i -220} {$i <= 220} {incr i} {
              add PhiBins [expr {$i * $pi/220}]
          }
  add EtaPhiBins -3.73 $PhiBins
  add EtaPhiBins 3.74 $PhiBins
  set PhiBins {}
          for {set i -224} {$i <= 224} {incr i} {
              add PhiBins [expr {$i * $pi/224}]
          }
  add EtaPhiBins -3.71 $PhiBins
  add EtaPhiBins 3.73 $PhiBins
  set PhiBins {}
          for {set i -227} {$i <= 227} {incr i} {
              add PhiBins [expr {$i * $pi/227}]
          }
  add EtaPhiBins -3.70 $PhiBins
  add EtaPhiBins 3.71 $PhiBins
  set PhiBins {}
          for {set i -230} {$i <= 230} {incr i} {
              add PhiBins [expr {$i * $pi/230}]
          }
  add EtaPhiBins -3.69 $PhiBins
  add EtaPhiBins 3.70 $PhiBins
  set PhiBins {}
          for {set i -233} {$i <= 233} {incr i} {
              add PhiBins [expr {$i * $pi/233}]
          }
  add EtaPhiBins -3.67 $PhiBins
  add EtaPhiBins 3.69 $PhiBins
  set PhiBins {}
          for {set i -236} {$i <= 236} {incr i} {
              add PhiBins [expr {$i * $pi/236}]
          }
  add EtaPhiBins -3.66 $PhiBins
  add EtaPhiBins 3.67 $PhiBins
  set PhiBins {}
          for {set i -239} {$i <= 239} {incr i} {
              add PhiBins [expr {$i * $pi/239}]
          }
  add EtaPhiBins -3.65 $PhiBins
  add EtaPhiBins 3.66 $PhiBins
  set PhiBins {}
          for {set i -242} {$i <= 242} {incr i} {
              add PhiBins [expr {$i * $pi/242}]
          }
  add EtaPhiBins -3.63 $PhiBins
  add EtaPhiBins 3.65 $PhiBins
  set PhiBins {}
          for {set i -246} {$i <= 246} {incr i} {
              add PhiBins [expr {$i * $pi/246}]
          }
  add EtaPhiBins -3.62 $PhiBins
  add EtaPhiBins 3.63 $PhiBins
  set PhiBins {}
          for {set i -249} {$i <= 249} {incr i} {
              add PhiBins [expr {$i * $pi/249}]
          }
  add EtaPhiBins -3.61 $PhiBins
  add EtaPhiBins 3.62 $PhiBins
  set PhiBins {}
          for {set i -252} {$i <= 252} {incr i} {
              add PhiBins [expr {$i * $pi/252}]
          }
  add EtaPhiBins -3.59 $PhiBins
  add EtaPhiBins 3.61 $PhiBins
  set PhiBins {}
          for {set i -255} {$i <= 255} {incr i} {
              add PhiBins [expr {$i * $pi/255}]
          }
  add EtaPhiBins -3.58 $PhiBins
  add EtaPhiBins 3.59 $PhiBins
  set PhiBins {}
          for {set i -258} {$i <= 258} {incr i} {
              add PhiBins [expr {$i * $pi/258}]
          }
  add EtaPhiBins -3.57 $PhiBins
  add EtaPhiBins 3.58 $PhiBins
  set PhiBins {}
          for {set i -261} {$i <= 261} {incr i} {
              add PhiBins [expr {$i * $pi/261}]
          }
  add EtaPhiBins -3.56 $PhiBins
  add EtaPhiBins 3.57 $PhiBins
  set PhiBins {}
          for {set i -264} {$i <= 264} {incr i} {
              add PhiBins [expr {$i * $pi/264}]
          }
  add EtaPhiBins -3.55 $PhiBins
  add EtaPhiBins 3.56 $PhiBins
  set PhiBins {}
          for {set i -268} {$i <= 268} {incr i} {
              add PhiBins [expr {$i * $pi/268}]
          }
  add EtaPhiBins -3.54 $PhiBins
  add EtaPhiBins 3.55 $PhiBins
  set PhiBins {}
          for {set i -271} {$i <= 271} {incr i} {
              add PhiBins [expr {$i * $pi/271}]
          }
  add EtaPhiBins -3.52 $PhiBins
  add EtaPhiBins 3.54 $PhiBins
  set PhiBins {}
          for {set i -274} {$i <= 274} {incr i} {
              add PhiBins [expr {$i * $pi/274}]
          }
  add EtaPhiBins -3.51 $PhiBins
  add EtaPhiBins 3.52 $PhiBins
  set PhiBins {}
          for {set i -277} {$i <= 277} {incr i} {
              add PhiBins [expr {$i * $pi/277}]
          }
  add EtaPhiBins -3.50 $PhiBins
  add EtaPhiBins 3.51 $PhiBins
  set PhiBins {}
          for {set i -280} {$i <= 280} {incr i} {
              add PhiBins [expr {$i * $pi/280}]
          }
  add EtaPhiBins -3.49 $PhiBins
  add EtaPhiBins 3.50 $PhiBins
  set PhiBins {}
          for {set i -283} {$i <= 283} {incr i} {
              add PhiBins [expr {$i * $pi/283}]
          }
  add EtaPhiBins -3.48 $PhiBins
  add EtaPhiBins 3.49 $PhiBins
  set PhiBins {}
          for {set i -286} {$i <= 286} {incr i} {
              add PhiBins [expr {$i * $pi/286}]
          }
  add EtaPhiBins -3.47 $PhiBins
  add EtaPhiBins 3.48 $PhiBins
  set PhiBins {}
          for {set i -290} {$i <= 290} {incr i} {
              add PhiBins [expr {$i * $pi/290}]
          }
  add EtaPhiBins -3.46 $PhiBins
  add EtaPhiBins 3.47 $PhiBins
  set PhiBins {}
          for {set i -293} {$i <= 293} {incr i} {
              add PhiBins [expr {$i * $pi/293}]
          }
  add EtaPhiBins -3.45 $PhiBins
  add EtaPhiBins 3.46 $PhiBins
  set PhiBins {}
          for {set i -296} {$i <= 296} {incr i} {
              add PhiBins [expr {$i * $pi/296}]
          }
  add EtaPhiBins -3.44 $PhiBins
  add EtaPhiBins 3.45 $PhiBins
  set PhiBins {}
          for {set i -299} {$i <= 299} {incr i} {
              add PhiBins [expr {$i * $pi/299}]
          }
  add EtaPhiBins -3.43 $PhiBins
  add EtaPhiBins 3.44 $PhiBins
  set PhiBins {}
          for {set i -302} {$i <= 302} {incr i} {
              add PhiBins [expr {$i * $pi/302}]
          }
  add EtaPhiBins -3.41 $PhiBins
  add EtaPhiBins 3.43 $PhiBins
  set PhiBins {}
          for {set i -305} {$i <= 305} {incr i} {
              add PhiBins [expr {$i * $pi/305}]
          }
  add EtaPhiBins -3.40 $PhiBins
  add EtaPhiBins 3.41 $PhiBins
  set PhiBins {}
          for {set i -308} {$i <= 308} {incr i} {
              add PhiBins [expr {$i * $pi/308}]
          }
  add EtaPhiBins -3.39 $PhiBins
  add EtaPhiBins 3.40 $PhiBins
  set PhiBins {}
          for {set i -312} {$i <= 312} {incr i} {
              add PhiBins [expr {$i * $pi/312}]
          }
  add EtaPhiBins -3.38 $PhiBins
  add EtaPhiBins 3.39 $PhiBins
  set PhiBins {}
          for {set i -315} {$i <= 315} {incr i} {
              add PhiBins [expr {$i * $pi/315}]
          }
  add EtaPhiBins -3.37 $PhiBins
  add EtaPhiBins 3.38 $PhiBins
  set PhiBins {}
          for {set i -318} {$i <= 318} {incr i} {
              add PhiBins [expr {$i * $pi/318}]
          }
  add EtaPhiBins -3.36 $PhiBins
  add EtaPhiBins 3.37 $PhiBins
  set PhiBins {}
          for {set i -321} {$i <= 321} {incr i} {
              add PhiBins [expr {$i * $pi/321}]
          }
  add EtaPhiBins -3.35 $PhiBins
  add EtaPhiBins 3.36 $PhiBins
  set PhiBins {}
          for {set i -324} {$i <= 324} {incr i} {
              add PhiBins [expr {$i * $pi/324}]
          }
  add EtaPhiBins -3.35 $PhiBins
  add EtaPhiBins 3.35 $PhiBins
  set PhiBins {}
          for {set i -327} {$i <= 327} {incr i} {
              add PhiBins [expr {$i * $pi/327}]
          }
  add EtaPhiBins -3.34 $PhiBins
  add EtaPhiBins 3.35 $PhiBins
  set PhiBins {}
          for {set i -330} {$i <= 330} {incr i} {
              add PhiBins [expr {$i * $pi/330}]
          }
  add EtaPhiBins -3.33 $PhiBins
  add EtaPhiBins 3.34 $PhiBins
  set PhiBins {}
          for {set i -334} {$i <= 334} {incr i} {
              add PhiBins [expr {$i * $pi/334}]
          }
  add EtaPhiBins -3.32 $PhiBins
  add EtaPhiBins 3.33 $PhiBins
  set PhiBins {}
          for {set i -337} {$i <= 337} {incr i} {
              add PhiBins [expr {$i * $pi/337}]
          }
  add EtaPhiBins -3.31 $PhiBins
  add EtaPhiBins 3.32 $PhiBins
  set PhiBins {}
          for {set i -340} {$i <= 340} {incr i} {
              add PhiBins [expr {$i * $pi/340}]
          }
  add EtaPhiBins -3.30 $PhiBins
  add EtaPhiBins 3.31 $PhiBins
  set PhiBins {}
          for {set i -343} {$i <= 343} {incr i} {
              add PhiBins [expr {$i * $pi/343}]
          }
  add EtaPhiBins -3.29 $PhiBins
  add EtaPhiBins 3.30 $PhiBins
  set PhiBins {}
          for {set i -346} {$i <= 346} {incr i} {
              add PhiBins [expr {$i * $pi/346}]
          }
  add EtaPhiBins -3.28 $PhiBins
  add EtaPhiBins 3.29 $PhiBins
  set PhiBins {}
          for {set i -349} {$i <= 349} {incr i} {
              add PhiBins [expr {$i * $pi/349}]
          }
  add EtaPhiBins -3.27 $PhiBins
  add EtaPhiBins 3.28 $PhiBins
  set PhiBins {}
          for {set i -352} {$i <= 352} {incr i} {
              add PhiBins [expr {$i * $pi/352}]
          }
  add EtaPhiBins -3.26 $PhiBins
  add EtaPhiBins 3.27 $PhiBins
  set PhiBins {}
          for {set i -356} {$i <= 356} {incr i} {
              add PhiBins [expr {$i * $pi/356}]
          }
  add EtaPhiBins -3.25 $PhiBins
  add EtaPhiBins 3.26 $PhiBins
  set PhiBins {}
          for {set i -359} {$i <= 359} {incr i} {
              add PhiBins [expr {$i * $pi/359}]
          }
  add EtaPhiBins -3.24 $PhiBins
  add EtaPhiBins 3.25 $PhiBins
  set PhiBins {}
          for {set i -362} {$i <= 362} {incr i} {
              add PhiBins [expr {$i * $pi/362}]
          }
  add EtaPhiBins -3.24 $PhiBins
  add EtaPhiBins 3.24 $PhiBins
  set PhiBins {}
          for {set i -365} {$i <= 365} {incr i} {
              add PhiBins [expr {$i * $pi/365}]
          }
  add EtaPhiBins -3.23 $PhiBins
  add EtaPhiBins 3.24 $PhiBins
  set PhiBins {}
          for {set i -368} {$i <= 368} {incr i} {
              add PhiBins [expr {$i * $pi/368}]
          }
  add EtaPhiBins -3.22 $PhiBins
  add EtaPhiBins 3.23 $PhiBins
  set PhiBins {}
          for {set i -371} {$i <= 371} {incr i} {
              add PhiBins [expr {$i * $pi/371}]
          }
  add EtaPhiBins -3.21 $PhiBins
  add EtaPhiBins 3.22 $PhiBins
  set PhiBins {}
          for {set i -374} {$i <= 374} {incr i} {
              add PhiBins [expr {$i * $pi/374}]
          }
  add EtaPhiBins -3.20 $PhiBins
  add EtaPhiBins 3.21 $PhiBins
  set PhiBins {}
          for {set i -378} {$i <= 378} {incr i} {
              add PhiBins [expr {$i * $pi/378}]
          }
  add EtaPhiBins -3.19 $PhiBins
  add EtaPhiBins 3.20 $PhiBins
  set PhiBins {}
          for {set i -381} {$i <= 381} {incr i} {
              add PhiBins [expr {$i * $pi/381}]
          }
  add EtaPhiBins -3.19 $PhiBins
  add EtaPhiBins 3.19 $PhiBins
  set PhiBins {}
          for {set i -384} {$i <= 384} {incr i} {
              add PhiBins [expr {$i * $pi/384}]
          }
  add EtaPhiBins -3.18 $PhiBins
  add EtaPhiBins 3.19 $PhiBins
  set PhiBins {}
          for {set i -387} {$i <= 387} {incr i} {
              add PhiBins [expr {$i * $pi/387}]
          }
  add EtaPhiBins -3.17 $PhiBins
  add EtaPhiBins 3.18 $PhiBins
  set PhiBins {}
          for {set i -390} {$i <= 390} {incr i} {
              add PhiBins [expr {$i * $pi/390}]
          }
  add EtaPhiBins -3.16 $PhiBins
  add EtaPhiBins 3.17 $PhiBins
  set PhiBins {}
          for {set i -393} {$i <= 393} {incr i} {
              add PhiBins [expr {$i * $pi/393}]
          }
  add EtaPhiBins -3.15 $PhiBins
  add EtaPhiBins 3.16 $PhiBins
  set PhiBins {}
          for {set i -396} {$i <= 396} {incr i} {
              add PhiBins [expr {$i * $pi/396}]
          }
  add EtaPhiBins -3.15 $PhiBins
  add EtaPhiBins 3.15 $PhiBins
  set PhiBins {}
          for {set i -400} {$i <= 400} {incr i} {
              add PhiBins [expr {$i * $pi/400}]
          }
  add EtaPhiBins -3.14 $PhiBins
  add EtaPhiBins 3.15 $PhiBins
  set PhiBins {}
          for {set i -403} {$i <= 403} {incr i} {
              add PhiBins [expr {$i * $pi/403}]
          }
  add EtaPhiBins -3.13 $PhiBins
  add EtaPhiBins 3.14 $PhiBins
  set PhiBins {}
          for {set i -406} {$i <= 406} {incr i} {
              add PhiBins [expr {$i * $pi/406}]
          }
  add EtaPhiBins -3.12 $PhiBins
  add EtaPhiBins 3.13 $PhiBins
  set PhiBins {}
          for {set i -409} {$i <= 409} {incr i} {
              add PhiBins [expr {$i * $pi/409}]
          }
  add EtaPhiBins -3.11 $PhiBins
  add EtaPhiBins 3.12 $PhiBins
  set PhiBins {}
          for {set i -412} {$i <= 412} {incr i} {
              add PhiBins [expr {$i * $pi/412}]
          }
  add EtaPhiBins -3.11 $PhiBins
  add EtaPhiBins 3.11 $PhiBins
  set PhiBins {}
          for {set i -415} {$i <= 415} {incr i} {
              add PhiBins [expr {$i * $pi/415}]
          }
  add EtaPhiBins -3.10 $PhiBins
  add EtaPhiBins 3.11 $PhiBins
  set PhiBins {}
          for {set i -418} {$i <= 418} {incr i} {
              add PhiBins [expr {$i * $pi/418}]
          }
  add EtaPhiBins -3.09 $PhiBins
  add EtaPhiBins 3.10 $PhiBins
  set PhiBins {}
          for {set i -422} {$i <= 422} {incr i} {
              add PhiBins [expr {$i * $pi/422}]
          }
  add EtaPhiBins -3.08 $PhiBins
  add EtaPhiBins 3.09 $PhiBins
  set PhiBins {}
          for {set i -425} {$i <= 425} {incr i} {
              add PhiBins [expr {$i * $pi/425}]
          }
  add EtaPhiBins -3.08 $PhiBins
  add EtaPhiBins 3.08 $PhiBins
  set PhiBins {}
          for {set i -428} {$i <= 428} {incr i} {
              add PhiBins [expr {$i * $pi/428}]
          }
  add EtaPhiBins -3.07 $PhiBins
  add EtaPhiBins 3.08 $PhiBins
  set PhiBins {}
          for {set i -431} {$i <= 431} {incr i} {
              add PhiBins [expr {$i * $pi/431}]
          }
  add EtaPhiBins -3.06 $PhiBins
  add EtaPhiBins 3.07 $PhiBins
  set PhiBins {}
          for {set i -434} {$i <= 434} {incr i} {
              add PhiBins [expr {$i * $pi/434}]
          }
  add EtaPhiBins -3.06 $PhiBins
  add EtaPhiBins 3.06 $PhiBins
  set PhiBins {}
          for {set i -437} {$i <= 437} {incr i} {
              add PhiBins [expr {$i * $pi/437}]
          }
  add EtaPhiBins -3.05 $PhiBins
  add EtaPhiBins 3.06 $PhiBins
  set PhiBins {}
          for {set i -440} {$i <= 440} {incr i} {
              add PhiBins [expr {$i * $pi/440}]
          }
  add EtaPhiBins -3.04 $PhiBins
  add EtaPhiBins 3.05 $PhiBins
  set PhiBins {}
          for {set i -444} {$i <= 444} {incr i} {
              add PhiBins [expr {$i * $pi/444}]
          }
  add EtaPhiBins -3.03 $PhiBins
  add EtaPhiBins 3.04 $PhiBins
  set PhiBins {}
          for {set i -447} {$i <= 447} {incr i} {
              add PhiBins [expr {$i * $pi/447}]
          }
  add EtaPhiBins -3.03 $PhiBins
  add EtaPhiBins 3.03 $PhiBins
  set PhiBins {}
          for {set i -450} {$i <= 450} {incr i} {
              add PhiBins [expr {$i * $pi/450}]
          }
  add EtaPhiBins -3.02 $PhiBins
  add EtaPhiBins 3.03 $PhiBins
  set PhiBins {}
          for {set i -453} {$i <= 453} {incr i} {
              add PhiBins [expr {$i * $pi/453}]
          }
  add EtaPhiBins -3.01 $PhiBins
  add EtaPhiBins 3.02 $PhiBins
  set PhiBins {}
          for {set i -456} {$i <= 456} {incr i} {
              add PhiBins [expr {$i * $pi/456}]
          }
  add EtaPhiBins -3.01 $PhiBins
  add EtaPhiBins 3.01 $PhiBins
  set PhiBins {}
          for {set i -459} {$i <= 459} {incr i} {
              add PhiBins [expr {$i * $pi/459}]
          }
  add EtaPhiBins -3.00 $PhiBins
  add EtaPhiBins 3.01 $PhiBins
  set PhiBins {}
          for {set i -462} {$i <= 462} {incr i} {
              add PhiBins [expr {$i * $pi/462}]
          }
  add EtaPhiBins -2.99 $PhiBins
  add EtaPhiBins 3.00 $PhiBins
  set PhiBins {}
          for {set i -466} {$i <= 466} {incr i} {
              add PhiBins [expr {$i * $pi/466}]
          }
  add EtaPhiBins -2.99 $PhiBins
  add EtaPhiBins 2.99 $PhiBins
  set PhiBins {}
          for {set i -161} {$i <= 161} {incr i} {
              add PhiBins [expr {$i * $pi/161}]
          }
  add EtaPhiBins -2.97 $PhiBins
  add EtaPhiBins 2.99 $PhiBins
  set PhiBins {}
          for {set i -164} {$i <= 164} {incr i} {
              add PhiBins [expr {$i * $pi/164}]
          }
  add EtaPhiBins -2.95 $PhiBins
  add EtaPhiBins 2.97 $PhiBins
  set PhiBins {}
          for {set i -167} {$i <= 167} {incr i} {
              add PhiBins [expr {$i * $pi/167}]
          }
  add EtaPhiBins -2.93 $PhiBins
  add EtaPhiBins 2.95 $PhiBins
  set PhiBins {}
          for {set i -170} {$i <= 170} {incr i} {
              add PhiBins [expr {$i * $pi/170}]
          }
  add EtaPhiBins -2.91 $PhiBins
  add EtaPhiBins 2.93 $PhiBins
  set PhiBins {}
          for {set i -174} {$i <= 174} {incr i} {
              add PhiBins [expr {$i * $pi/174}]
          }
  add EtaPhiBins -2.90 $PhiBins
  add EtaPhiBins 2.91 $PhiBins
  set PhiBins {}
          for {set i -177} {$i <= 177} {incr i} {
              add PhiBins [expr {$i * $pi/177}]
          }
  add EtaPhiBins -2.88 $PhiBins
  add EtaPhiBins 2.90 $PhiBins
  set PhiBins {}
          for {set i -180} {$i <= 180} {incr i} {
              add PhiBins [expr {$i * $pi/180}]
          }
  add EtaPhiBins -2.86 $PhiBins
  add EtaPhiBins 2.88 $PhiBins
  set PhiBins {}
          for {set i -183} {$i <= 183} {incr i} {
              add PhiBins [expr {$i * $pi/183}]
          }
  add EtaPhiBins -2.84 $PhiBins
  add EtaPhiBins 2.86 $PhiBins
  set PhiBins {}
          for {set i -186} {$i <= 186} {incr i} {
              add PhiBins [expr {$i * $pi/186}]
          }
  add EtaPhiBins -2.83 $PhiBins
  add EtaPhiBins 2.84 $PhiBins
  set PhiBins {}
          for {set i -189} {$i <= 189} {incr i} {
              add PhiBins [expr {$i * $pi/189}]
          }
  add EtaPhiBins -2.81 $PhiBins
  add EtaPhiBins 2.83 $PhiBins
  set PhiBins {}
          for {set i -192} {$i <= 192} {incr i} {
              add PhiBins [expr {$i * $pi/192}]
          }
  add EtaPhiBins -2.80 $PhiBins
  add EtaPhiBins 2.81 $PhiBins
  set PhiBins {}
          for {set i -196} {$i <= 196} {incr i} {
              add PhiBins [expr {$i * $pi/196}]
          }
  add EtaPhiBins -2.78 $PhiBins
  add EtaPhiBins 2.80 $PhiBins
  set PhiBins {}
          for {set i -199} {$i <= 199} {incr i} {
              add PhiBins [expr {$i * $pi/199}]
          }
  add EtaPhiBins -2.76 $PhiBins
  add EtaPhiBins 2.78 $PhiBins
  set PhiBins {}
          for {set i -202} {$i <= 202} {incr i} {
              add PhiBins [expr {$i * $pi/202}]
          }
  add EtaPhiBins -2.75 $PhiBins
  add EtaPhiBins 2.76 $PhiBins
  set PhiBins {}
          for {set i -205} {$i <= 205} {incr i} {
              add PhiBins [expr {$i * $pi/205}]
          }
  add EtaPhiBins -2.73 $PhiBins
  add EtaPhiBins 2.75 $PhiBins
  set PhiBins {}
          for {set i -208} {$i <= 208} {incr i} {
              add PhiBins [expr {$i * $pi/208}]
          }
  add EtaPhiBins -2.72 $PhiBins
  add EtaPhiBins 2.73 $PhiBins
  set PhiBins {}
          for {set i -211} {$i <= 211} {incr i} {
              add PhiBins [expr {$i * $pi/211}]
          }
  add EtaPhiBins -2.70 $PhiBins
  add EtaPhiBins 2.72 $PhiBins
  set PhiBins {}
          for {set i -214} {$i <= 214} {incr i} {
              add PhiBins [expr {$i * $pi/214}]
          }
  add EtaPhiBins -2.69 $PhiBins
  add EtaPhiBins 2.70 $PhiBins
  set PhiBins {}
          for {set i -218} {$i <= 218} {incr i} {
              add PhiBins [expr {$i * $pi/218}]
          }
  add EtaPhiBins -2.67 $PhiBins
  add EtaPhiBins 2.69 $PhiBins
  set PhiBins {}
          for {set i -221} {$i <= 221} {incr i} {
              add PhiBins [expr {$i * $pi/221}]
          }
  add EtaPhiBins -2.66 $PhiBins
  add EtaPhiBins 2.67 $PhiBins
  set PhiBins {}
          for {set i -224} {$i <= 224} {incr i} {
              add PhiBins [expr {$i * $pi/224}]
          }
  add EtaPhiBins -2.65 $PhiBins
  add EtaPhiBins 2.66 $PhiBins
  set PhiBins {}
          for {set i -227} {$i <= 227} {incr i} {
              add PhiBins [expr {$i * $pi/227}]
          }
  add EtaPhiBins -2.63 $PhiBins
  add EtaPhiBins 2.65 $PhiBins
  set PhiBins {}
          for {set i -230} {$i <= 230} {incr i} {
              add PhiBins [expr {$i * $pi/230}]
          }
  add EtaPhiBins -2.62 $PhiBins
  add EtaPhiBins 2.63 $PhiBins
  set PhiBins {}
          for {set i -233} {$i <= 233} {incr i} {
              add PhiBins [expr {$i * $pi/233}]
          }
  add EtaPhiBins -2.61 $PhiBins
  add EtaPhiBins 2.62 $PhiBins
  set PhiBins {}
          for {set i -236} {$i <= 236} {incr i} {
              add PhiBins [expr {$i * $pi/236}]
          }
  add EtaPhiBins -2.59 $PhiBins
  add EtaPhiBins 2.61 $PhiBins
  set PhiBins {}
          for {set i -240} {$i <= 240} {incr i} {
              add PhiBins [expr {$i * $pi/240}]
          }
  add EtaPhiBins -2.58 $PhiBins
  add EtaPhiBins 2.59 $PhiBins
  set PhiBins {}
          for {set i -243} {$i <= 243} {incr i} {
              add PhiBins [expr {$i * $pi/243}]
          }
  add EtaPhiBins -2.57 $PhiBins
  add EtaPhiBins 2.58 $PhiBins
  set PhiBins {}
          for {set i -246} {$i <= 246} {incr i} {
              add PhiBins [expr {$i * $pi/246}]
          }
  add EtaPhiBins -2.56 $PhiBins
  add EtaPhiBins 2.57 $PhiBins
  set PhiBins {}
          for {set i -249} {$i <= 249} {incr i} {
              add PhiBins [expr {$i * $pi/249}]
          }
  add EtaPhiBins -2.54 $PhiBins
  add EtaPhiBins 2.56 $PhiBins
  set PhiBins {}
          for {set i -252} {$i <= 252} {incr i} {
              add PhiBins [expr {$i * $pi/252}]
          }
  add EtaPhiBins -2.53 $PhiBins
  add EtaPhiBins 2.54 $PhiBins
  set PhiBins {}
          for {set i -255} {$i <= 255} {incr i} {
              add PhiBins [expr {$i * $pi/255}]
          }
  add EtaPhiBins -2.52 $PhiBins
  add EtaPhiBins 2.53 $PhiBins
  set PhiBins {}
          for {set i -258} {$i <= 258} {incr i} {
              add PhiBins [expr {$i * $pi/258}]
          }
  add EtaPhiBins -2.51 $PhiBins
  add EtaPhiBins 2.52 $PhiBins
  set PhiBins {}
          for {set i -262} {$i <= 262} {incr i} {
              add PhiBins [expr {$i * $pi/262}]
          }
  add EtaPhiBins -2.49 $PhiBins
  add EtaPhiBins 2.51 $PhiBins
  set PhiBins {}
          for {set i -265} {$i <= 265} {incr i} {
              add PhiBins [expr {$i * $pi/265}]
          }
  add EtaPhiBins -2.48 $PhiBins
  add EtaPhiBins 2.49 $PhiBins
  set PhiBins {}
          for {set i -268} {$i <= 268} {incr i} {
              add PhiBins [expr {$i * $pi/268}]
          }
  add EtaPhiBins -2.47 $PhiBins
  add EtaPhiBins 2.48 $PhiBins
  set PhiBins {}
          for {set i -271} {$i <= 271} {incr i} {
              add PhiBins [expr {$i * $pi/271}]
          }
  add EtaPhiBins -2.46 $PhiBins
  add EtaPhiBins 2.47 $PhiBins
  set PhiBins {}
          for {set i -274} {$i <= 274} {incr i} {
              add PhiBins [expr {$i * $pi/274}]
          }
  add EtaPhiBins -2.45 $PhiBins
  add EtaPhiBins 2.46 $PhiBins
  set PhiBins {}
          for {set i -277} {$i <= 277} {incr i} {
              add PhiBins [expr {$i * $pi/277}]
          }
  add EtaPhiBins -2.44 $PhiBins
  add EtaPhiBins 2.45 $PhiBins
  set PhiBins {}
          for {set i -280} {$i <= 280} {incr i} {
              add PhiBins [expr {$i * $pi/280}]
          }
  add EtaPhiBins -2.43 $PhiBins
  add EtaPhiBins 2.44 $PhiBins
  set PhiBins {}
          for {set i -284} {$i <= 284} {incr i} {
              add PhiBins [expr {$i * $pi/284}]
          }
  add EtaPhiBins -2.42 $PhiBins
  add EtaPhiBins 2.43 $PhiBins
  set PhiBins {}
          for {set i -287} {$i <= 287} {incr i} {
              add PhiBins [expr {$i * $pi/287}]
          }
  add EtaPhiBins -2.40 $PhiBins
  add EtaPhiBins 2.42 $PhiBins
  set PhiBins {}
          for {set i -290} {$i <= 290} {incr i} {
              add PhiBins [expr {$i * $pi/290}]
          }
  add EtaPhiBins -2.39 $PhiBins
  add EtaPhiBins 2.40 $PhiBins
  set PhiBins {}
          for {set i -293} {$i <= 293} {incr i} {
              add PhiBins [expr {$i * $pi/293}]
          }
  add EtaPhiBins -2.38 $PhiBins
  add EtaPhiBins 2.39 $PhiBins
  set PhiBins {}
          for {set i -296} {$i <= 296} {incr i} {
              add PhiBins [expr {$i * $pi/296}]
          }
  add EtaPhiBins -2.37 $PhiBins
  add EtaPhiBins 2.38 $PhiBins
  set PhiBins {}
          for {set i -299} {$i <= 299} {incr i} {
              add PhiBins [expr {$i * $pi/299}]
          }
  add EtaPhiBins -2.36 $PhiBins
  add EtaPhiBins 2.37 $PhiBins
  set PhiBins {}
          for {set i -302} {$i <= 302} {incr i} {
              add PhiBins [expr {$i * $pi/302}]
          }
  add EtaPhiBins -2.35 $PhiBins
  add EtaPhiBins 2.36 $PhiBins
  set PhiBins {}
          for {set i -306} {$i <= 306} {incr i} {
              add PhiBins [expr {$i * $pi/306}]
          }
  add EtaPhiBins -2.34 $PhiBins
  add EtaPhiBins 2.35 $PhiBins
  set PhiBins {}
          for {set i -309} {$i <= 309} {incr i} {
              add PhiBins [expr {$i * $pi/309}]
          }
  add EtaPhiBins -2.33 $PhiBins
  add EtaPhiBins 2.34 $PhiBins
  set PhiBins {}
          for {set i -312} {$i <= 312} {incr i} {
              add PhiBins [expr {$i * $pi/312}]
          }
  add EtaPhiBins -2.32 $PhiBins
  add EtaPhiBins 2.33 $PhiBins
  set PhiBins {}
          for {set i -315} {$i <= 315} {incr i} {
              add PhiBins [expr {$i * $pi/315}]
          }
  add EtaPhiBins -2.31 $PhiBins
  add EtaPhiBins 2.32 $PhiBins
  set PhiBins {}
          for {set i -318} {$i <= 318} {incr i} {
              add PhiBins [expr {$i * $pi/318}]
          }
  add EtaPhiBins -2.30 $PhiBins
  add EtaPhiBins 2.31 $PhiBins
  set PhiBins {}
          for {set i -321} {$i <= 321} {incr i} {
              add PhiBins [expr {$i * $pi/321}]
          }
  add EtaPhiBins -2.29 $PhiBins
  add EtaPhiBins 2.30 $PhiBins
  set PhiBins {}
          for {set i -324} {$i <= 324} {incr i} {
              add PhiBins [expr {$i * $pi/324}]
          }
  add EtaPhiBins -2.28 $PhiBins
  add EtaPhiBins 2.29 $PhiBins
  set PhiBins {}
          for {set i -328} {$i <= 328} {incr i} {
              add PhiBins [expr {$i * $pi/328}]
          }
  add EtaPhiBins -2.27 $PhiBins
  add EtaPhiBins 2.28 $PhiBins
  set PhiBins {}
          for {set i -331} {$i <= 331} {incr i} {
              add PhiBins [expr {$i * $pi/331}]
          }
  add EtaPhiBins -2.27 $PhiBins
  add EtaPhiBins 2.27 $PhiBins
  set PhiBins {}
          for {set i -334} {$i <= 334} {incr i} {
              add PhiBins [expr {$i * $pi/334}]
          }
  add EtaPhiBins -2.26 $PhiBins
  add EtaPhiBins 2.27 $PhiBins
  set PhiBins {}
          for {set i -337} {$i <= 337} {incr i} {
              add PhiBins [expr {$i * $pi/337}]
          }
  add EtaPhiBins -2.25 $PhiBins
  add EtaPhiBins 2.26 $PhiBins
  set PhiBins {}
          for {set i -340} {$i <= 340} {incr i} {
              add PhiBins [expr {$i * $pi/340}]
          }
  add EtaPhiBins -2.24 $PhiBins
  add EtaPhiBins 2.25 $PhiBins
  set PhiBins {}
          for {set i -343} {$i <= 343} {incr i} {
              add PhiBins [expr {$i * $pi/343}]
          }
  add EtaPhiBins -2.23 $PhiBins
  add EtaPhiBins 2.24 $PhiBins
  set PhiBins {}
          for {set i -346} {$i <= 346} {incr i} {
              add PhiBins [expr {$i * $pi/346}]
          }
  add EtaPhiBins -2.22 $PhiBins
  add EtaPhiBins 2.23 $PhiBins
  set PhiBins {}
          for {set i -350} {$i <= 350} {incr i} {
              add PhiBins [expr {$i * $pi/350}]
          }
  add EtaPhiBins -2.21 $PhiBins
  add EtaPhiBins 2.22 $PhiBins
  set PhiBins {}
          for {set i -353} {$i <= 353} {incr i} {
              add PhiBins [expr {$i * $pi/353}]
          }
  add EtaPhiBins -2.20 $PhiBins
  add EtaPhiBins 2.21 $PhiBins
  set PhiBins {}
          for {set i -356} {$i <= 356} {incr i} {
              add PhiBins [expr {$i * $pi/356}]
          }
  add EtaPhiBins -2.20 $PhiBins
  add EtaPhiBins 2.20 $PhiBins
  set PhiBins {}
          for {set i -252} {$i <= 252} {incr i} {
              add PhiBins [expr {$i * $pi/252}]
          }
  add EtaPhiBins -2.19 $PhiBins
  add EtaPhiBins 2.20 $PhiBins
  set PhiBins {}
          for {set i -256} {$i <= 256} {incr i} {
              add PhiBins [expr {$i * $pi/256}]
          }
  add EtaPhiBins -2.18 $PhiBins
  add EtaPhiBins 2.19 $PhiBins
  set PhiBins {}
          for {set i -259} {$i <= 259} {incr i} {
              add PhiBins [expr {$i * $pi/259}]
          }
  add EtaPhiBins -2.17 $PhiBins
  add EtaPhiBins 2.18 $PhiBins
  set PhiBins {}
          for {set i -262} {$i <= 262} {incr i} {
              add PhiBins [expr {$i * $pi/262}]
          }
  add EtaPhiBins -2.15 $PhiBins
  add EtaPhiBins 2.17 $PhiBins
  set PhiBins {}
          for {set i -265} {$i <= 265} {incr i} {
              add PhiBins [expr {$i * $pi/265}]
          }
  add EtaPhiBins -2.14 $PhiBins
  add EtaPhiBins 2.15 $PhiBins
  set PhiBins {}
          for {set i -268} {$i <= 268} {incr i} {
              add PhiBins [expr {$i * $pi/268}]
          }
  add EtaPhiBins -2.13 $PhiBins
  add EtaPhiBins 2.14 $PhiBins
  set PhiBins {}
          for {set i -271} {$i <= 271} {incr i} {
              add PhiBins [expr {$i * $pi/271}]
          }
  add EtaPhiBins -2.12 $PhiBins
  add EtaPhiBins 2.13 $PhiBins
  set PhiBins {}
          for {set i -274} {$i <= 274} {incr i} {
              add PhiBins [expr {$i * $pi/274}]
          }
  add EtaPhiBins -2.11 $PhiBins
  add EtaPhiBins 2.12 $PhiBins
  set PhiBins {}
          for {set i -278} {$i <= 278} {incr i} {
              add PhiBins [expr {$i * $pi/278}]
          }
  add EtaPhiBins -2.10 $PhiBins
  add EtaPhiBins 2.11 $PhiBins
  set PhiBins {}
          for {set i -281} {$i <= 281} {incr i} {
              add PhiBins [expr {$i * $pi/281}]
          }
  add EtaPhiBins -2.09 $PhiBins
  add EtaPhiBins 2.10 $PhiBins
  set PhiBins {}
          for {set i -284} {$i <= 284} {incr i} {
              add PhiBins [expr {$i * $pi/284}]
          }
  add EtaPhiBins -2.08 $PhiBins
  add EtaPhiBins 2.09 $PhiBins
  set PhiBins {}
          for {set i -287} {$i <= 287} {incr i} {
              add PhiBins [expr {$i * $pi/287}]
          }
  add EtaPhiBins -2.07 $PhiBins
  add EtaPhiBins 2.08 $PhiBins
  set PhiBins {}
          for {set i -290} {$i <= 290} {incr i} {
              add PhiBins [expr {$i * $pi/290}]
          }
  add EtaPhiBins -2.05 $PhiBins
  add EtaPhiBins 2.07 $PhiBins
  set PhiBins {}
          for {set i -293} {$i <= 293} {incr i} {
              add PhiBins [expr {$i * $pi/293}]
          }
  add EtaPhiBins -2.04 $PhiBins
  add EtaPhiBins 2.05 $PhiBins
  set PhiBins {}
          for {set i -296} {$i <= 296} {incr i} {
              add PhiBins [expr {$i * $pi/296}]
          }
  add EtaPhiBins -2.03 $PhiBins
  add EtaPhiBins 2.04 $PhiBins
  set PhiBins {}
          for {set i -300} {$i <= 300} {incr i} {
              add PhiBins [expr {$i * $pi/300}]
          }
  add EtaPhiBins -2.02 $PhiBins
  add EtaPhiBins 2.03 $PhiBins
  set PhiBins {}
          for {set i -303} {$i <= 303} {incr i} {
              add PhiBins [expr {$i * $pi/303}]
          }
  add EtaPhiBins -2.01 $PhiBins
  add EtaPhiBins 2.02 $PhiBins
  set PhiBins {}
          for {set i -306} {$i <= 306} {incr i} {
              add PhiBins [expr {$i * $pi/306}]
          }
  add EtaPhiBins -2.00 $PhiBins
  add EtaPhiBins 2.01 $PhiBins
  set PhiBins {}
          for {set i -309} {$i <= 309} {incr i} {
              add PhiBins [expr {$i * $pi/309}]
          }
  add EtaPhiBins -1.99 $PhiBins
  add EtaPhiBins 2.00 $PhiBins
  set PhiBins {}
          for {set i -312} {$i <= 312} {incr i} {
              add PhiBins [expr {$i * $pi/312}]
          }
  add EtaPhiBins -1.98 $PhiBins
  add EtaPhiBins 1.99 $PhiBins
  set PhiBins {}
          for {set i -315} {$i <= 315} {incr i} {
              add PhiBins [expr {$i * $pi/315}]
          }
  add EtaPhiBins -1.98 $PhiBins
  add EtaPhiBins 1.98 $PhiBins
  set PhiBins {}
          for {set i -318} {$i <= 318} {incr i} {
              add PhiBins [expr {$i * $pi/318}]
          }
  add EtaPhiBins -1.97 $PhiBins
  add EtaPhiBins 1.98 $PhiBins
  set PhiBins {}
          for {set i -322} {$i <= 322} {incr i} {
              add PhiBins [expr {$i * $pi/322}]
          }
  add EtaPhiBins -1.96 $PhiBins
  add EtaPhiBins 1.97 $PhiBins
  set PhiBins {}
          for {set i -325} {$i <= 325} {incr i} {
              add PhiBins [expr {$i * $pi/325}]
          }
  add EtaPhiBins -1.95 $PhiBins
  add EtaPhiBins 1.96 $PhiBins
  set PhiBins {}
          for {set i -328} {$i <= 328} {incr i} {
              add PhiBins [expr {$i * $pi/328}]
          }
  add EtaPhiBins -1.94 $PhiBins
  add EtaPhiBins 1.95 $PhiBins
  set PhiBins {}
          for {set i -331} {$i <= 331} {incr i} {
              add PhiBins [expr {$i * $pi/331}]
          }
  add EtaPhiBins -1.93 $PhiBins
  add EtaPhiBins 1.94 $PhiBins
  set PhiBins {}
          for {set i -334} {$i <= 334} {incr i} {
              add PhiBins [expr {$i * $pi/334}]
          }
  add EtaPhiBins -1.92 $PhiBins
  add EtaPhiBins 1.93 $PhiBins
  set PhiBins {}
          for {set i -337} {$i <= 337} {incr i} {
              add PhiBins [expr {$i * $pi/337}]
          }
  add EtaPhiBins -1.91 $PhiBins
  add EtaPhiBins 1.92 $PhiBins
  set PhiBins {}
          for {set i -340} {$i <= 340} {incr i} {
              add PhiBins [expr {$i * $pi/340}]
          }
  add EtaPhiBins -1.90 $PhiBins
  add EtaPhiBins 1.91 $PhiBins
  set PhiBins {}
          for {set i -344} {$i <= 344} {incr i} {
              add PhiBins [expr {$i * $pi/344}]
          }
  add EtaPhiBins -1.89 $PhiBins
  add EtaPhiBins 1.90 $PhiBins
  set PhiBins {}
          for {set i -347} {$i <= 347} {incr i} {
              add PhiBins [expr {$i * $pi/347}]
          }
  add EtaPhiBins -1.88 $PhiBins
  add EtaPhiBins 1.89 $PhiBins
  set PhiBins {}
          for {set i -350} {$i <= 350} {incr i} {
              add PhiBins [expr {$i * $pi/350}]
          }
  add EtaPhiBins -1.88 $PhiBins
  add EtaPhiBins 1.88 $PhiBins
  set PhiBins {}
          for {set i -353} {$i <= 353} {incr i} {
              add PhiBins [expr {$i * $pi/353}]
          }
  add EtaPhiBins -1.87 $PhiBins
  add EtaPhiBins 1.88 $PhiBins
  set PhiBins {}
          for {set i -356} {$i <= 356} {incr i} {
              add PhiBins [expr {$i * $pi/356}]
          }
  add EtaPhiBins -1.86 $PhiBins
  add EtaPhiBins 1.87 $PhiBins
  set PhiBins {}
          for {set i -359} {$i <= 359} {incr i} {
              add PhiBins [expr {$i * $pi/359}]
          }
  add EtaPhiBins -1.85 $PhiBins
  add EtaPhiBins 1.86 $PhiBins
  set PhiBins {}
          for {set i -362} {$i <= 362} {incr i} {
              add PhiBins [expr {$i * $pi/362}]
          }
  add EtaPhiBins -1.84 $PhiBins
  add EtaPhiBins 1.85 $PhiBins
  set PhiBins {}
          for {set i -365} {$i <= 365} {incr i} {
              add PhiBins [expr {$i * $pi/365}]
          }
  add EtaPhiBins -1.83 $PhiBins
  add EtaPhiBins 1.84 $PhiBins
  set PhiBins {}
          for {set i -369} {$i <= 369} {incr i} {
              add PhiBins [expr {$i * $pi/369}]
          }
  add EtaPhiBins -1.83 $PhiBins
  add EtaPhiBins 1.83 $PhiBins
  set PhiBins {}
          for {set i -372} {$i <= 372} {incr i} {
              add PhiBins [expr {$i * $pi/372}]
          }
  add EtaPhiBins -1.82 $PhiBins
  add EtaPhiBins 1.83 $PhiBins
  set PhiBins {}
          for {set i -375} {$i <= 375} {incr i} {
              add PhiBins [expr {$i * $pi/375}]
          }
  add EtaPhiBins -1.81 $PhiBins
  add EtaPhiBins 1.82 $PhiBins
  set PhiBins {}
          for {set i -378} {$i <= 378} {incr i} {
              add PhiBins [expr {$i * $pi/378}]
          }
  add EtaPhiBins -1.80 $PhiBins
  add EtaPhiBins 1.81 $PhiBins
  set PhiBins {}
          for {set i -381} {$i <= 381} {incr i} {
              add PhiBins [expr {$i * $pi/381}]
          }
  add EtaPhiBins -1.79 $PhiBins
  add EtaPhiBins 1.80 $PhiBins
  set PhiBins {}
          for {set i -384} {$i <= 384} {incr i} {
              add PhiBins [expr {$i * $pi/384}]
          }
  add EtaPhiBins -1.79 $PhiBins
  add EtaPhiBins 1.79 $PhiBins
  set PhiBins {}
          for {set i -387} {$i <= 387} {incr i} {
              add PhiBins [expr {$i * $pi/387}]
          }
  add EtaPhiBins -1.78 $PhiBins
  add EtaPhiBins 1.79 $PhiBins
  set PhiBins {}
          for {set i -391} {$i <= 391} {incr i} {
              add PhiBins [expr {$i * $pi/391}]
          }
  add EtaPhiBins -1.77 $PhiBins
  add EtaPhiBins 1.78 $PhiBins
  set PhiBins {}
          for {set i -394} {$i <= 394} {incr i} {
              add PhiBins [expr {$i * $pi/394}]
          }
  add EtaPhiBins -1.76 $PhiBins
  add EtaPhiBins 1.77 $PhiBins
  set PhiBins {}
          for {set i -397} {$i <= 397} {incr i} {
              add PhiBins [expr {$i * $pi/397}]
          }
  add EtaPhiBins -1.76 $PhiBins
  add EtaPhiBins 1.76 $PhiBins
  set PhiBins {}
          for {set i -400} {$i <= 400} {incr i} {
              add PhiBins [expr {$i * $pi/400}]
          }
  add EtaPhiBins -1.75 $PhiBins
  add EtaPhiBins 1.76 $PhiBins
  set PhiBins {}
          for {set i -403} {$i <= 403} {incr i} {
              add PhiBins [expr {$i * $pi/403}]
          }
  add EtaPhiBins -1.74 $PhiBins
  add EtaPhiBins 1.75 $PhiBins
  set PhiBins {}
          for {set i -406} {$i <= 406} {incr i} {
              add PhiBins [expr {$i * $pi/406}]
          }
  add EtaPhiBins -1.73 $PhiBins
  add EtaPhiBins 1.74 $PhiBins
  set PhiBins {}
          for {set i -409} {$i <= 409} {incr i} {
              add PhiBins [expr {$i * $pi/409}]
          }
  add EtaPhiBins -1.73 $PhiBins
  add EtaPhiBins 1.73 $PhiBins
  set PhiBins {}
          for {set i -413} {$i <= 413} {incr i} {
              add PhiBins [expr {$i * $pi/413}]
          }
  add EtaPhiBins -1.72 $PhiBins
  add EtaPhiBins 1.73 $PhiBins
  set PhiBins {}
          for {set i -416} {$i <= 416} {incr i} {
              add PhiBins [expr {$i * $pi/416}]
          }
  add EtaPhiBins -1.71 $PhiBins
  add EtaPhiBins 1.72 $PhiBins
  set PhiBins {}
          for {set i -419} {$i <= 419} {incr i} {
              add PhiBins [expr {$i * $pi/419}]
          }
  add EtaPhiBins -1.71 $PhiBins
  add EtaPhiBins 1.71 $PhiBins
  set PhiBins {}
          for {set i -422} {$i <= 422} {incr i} {
              add PhiBins [expr {$i * $pi/422}]
          }
  add EtaPhiBins -1.70 $PhiBins
  add EtaPhiBins 1.71 $PhiBins
  set PhiBins {}
          for {set i -425} {$i <= 425} {incr i} {
              add PhiBins [expr {$i * $pi/425}]
          }
  add EtaPhiBins -1.69 $PhiBins
  add EtaPhiBins 1.70 $PhiBins
  set PhiBins {}
          for {set i -428} {$i <= 428} {incr i} {
              add PhiBins [expr {$i * $pi/428}]
          }
  add EtaPhiBins -1.69 $PhiBins
  add EtaPhiBins 1.69 $PhiBins
  set PhiBins {}
          for {set i -431} {$i <= 431} {incr i} {
              add PhiBins [expr {$i * $pi/431}]
          }
  add EtaPhiBins -1.68 $PhiBins
  add EtaPhiBins 1.69 $PhiBins
  set PhiBins {}
          for {set i -435} {$i <= 435} {incr i} {
              add PhiBins [expr {$i * $pi/435}]
          }
  add EtaPhiBins -1.67 $PhiBins
  add EtaPhiBins 1.68 $PhiBins
  set PhiBins {}
          for {set i -438} {$i <= 438} {incr i} {
              add PhiBins [expr {$i * $pi/438}]
          }
  add EtaPhiBins -1.67 $PhiBins
  add EtaPhiBins 1.67 $PhiBins
  set PhiBins {}
          for {set i -441} {$i <= 441} {incr i} {
              add PhiBins [expr {$i * $pi/441}]
          }
  add EtaPhiBins -1.66 $PhiBins
  add EtaPhiBins 1.67 $PhiBins
  set PhiBins {}
          for {set i -444} {$i <= 444} {incr i} {
              add PhiBins [expr {$i * $pi/444}]
          }
  add EtaPhiBins -1.65 $PhiBins
  add EtaPhiBins 1.66 $PhiBins
  set PhiBins {}
          for {set i -447} {$i <= 447} {incr i} {
              add PhiBins [expr {$i * $pi/447}]
          }
  add EtaPhiBins -1.65 $PhiBins
  add EtaPhiBins 1.65 $PhiBins
  set PhiBins {}
          for {set i -450} {$i <= 450} {incr i} {
              add PhiBins [expr {$i * $pi/450}]
          }
  add EtaPhiBins -1.64 $PhiBins
  add EtaPhiBins 1.65 $PhiBins
  set PhiBins {}
          for {set i -453} {$i <= 453} {incr i} {
              add PhiBins [expr {$i * $pi/453}]
          }
  add EtaPhiBins -1.63 $PhiBins
  add EtaPhiBins 1.64 $PhiBins
  set PhiBins {}
          for {set i -457} {$i <= 457} {incr i} {
              add PhiBins [expr {$i * $pi/457}]
          }
  add EtaPhiBins -1.63 $PhiBins
  add EtaPhiBins 1.63 $PhiBins
  set PhiBins {}
          for {set i -460} {$i <= 460} {incr i} {
              add PhiBins [expr {$i * $pi/460}]
          }
  add EtaPhiBins -1.62 $PhiBins
  add EtaPhiBins 1.63 $PhiBins
  set PhiBins {}
          for {set i -463} {$i <= 463} {incr i} {
              add PhiBins [expr {$i * $pi/463}]
          }
  add EtaPhiBins -1.61 $PhiBins
  add EtaPhiBins 1.62 $PhiBins
  set PhiBins {}
          for {set i -466} {$i <= 466} {incr i} {
              add PhiBins [expr {$i * $pi/466}]
          }
  add EtaPhiBins -1.61 $PhiBins
  add EtaPhiBins 1.61 $PhiBins
  set PhiBins {}
          for {set i -469} {$i <= 469} {incr i} {
              add PhiBins [expr {$i * $pi/469}]
          }
  add EtaPhiBins -1.60 $PhiBins
  add EtaPhiBins 1.61 $PhiBins
  set PhiBins {}
          for {set i -472} {$i <= 472} {incr i} {
              add PhiBins [expr {$i * $pi/472}]
          }
  add EtaPhiBins -1.60 $PhiBins
  add EtaPhiBins 1.60 $PhiBins
  set PhiBins {}
          for {set i -475} {$i <= 475} {incr i} {
              add PhiBins [expr {$i * $pi/475}]
          }
  add EtaPhiBins -1.59 $PhiBins
  add EtaPhiBins 1.60 $PhiBins
  set PhiBins {}
          for {set i -479} {$i <= 479} {incr i} {
              add PhiBins [expr {$i * $pi/479}]
          }
  add EtaPhiBins -1.58 $PhiBins
  add EtaPhiBins 1.59 $PhiBins
  set PhiBins {}
          for {set i -482} {$i <= 482} {incr i} {
              add PhiBins [expr {$i * $pi/482}]
          }
  add EtaPhiBins -1.58 $PhiBins
  add EtaPhiBins 1.58 $PhiBins
  set PhiBins {}
          for {set i -485} {$i <= 485} {incr i} {
              add PhiBins [expr {$i * $pi/485}]
          }
  add EtaPhiBins -1.57 $PhiBins
  add EtaPhiBins 1.58 $PhiBins
  set PhiBins {}
          for {set i -488} {$i <= 488} {incr i} {
              add PhiBins [expr {$i * $pi/488}]
          }
  add EtaPhiBins -1.57 $PhiBins
  add EtaPhiBins 1.57 $PhiBins
  set PhiBins {}
          for {set i -491} {$i <= 491} {incr i} {
              add PhiBins [expr {$i * $pi/491}]
          }
  add EtaPhiBins -1.56 $PhiBins
  add EtaPhiBins 1.57 $PhiBins
  set PhiBins {}
          for {set i -494} {$i <= 494} {incr i} {
              add PhiBins [expr {$i * $pi/494}]
          }
  add EtaPhiBins -1.55 $PhiBins
  add EtaPhiBins 1.56 $PhiBins
  set PhiBins {}
          for {set i -497} {$i <= 497} {incr i} {
              add PhiBins [expr {$i * $pi/497}]
          }
  add EtaPhiBins -1.55 $PhiBins
  add EtaPhiBins 1.55 $PhiBins
  set PhiBins {}
          for {set i -501} {$i <= 501} {incr i} {
              add PhiBins [expr {$i * $pi/501}]
          }
  add EtaPhiBins -1.54 $PhiBins
  add EtaPhiBins 1.55 $PhiBins
  set PhiBins {}
          for {set i -504} {$i <= 504} {incr i} {
              add PhiBins [expr {$i * $pi/504}]
          }
  add EtaPhiBins -1.54 $PhiBins
  add EtaPhiBins 1.54 $PhiBins
  set PhiBins {}
          for {set i -507} {$i <= 507} {incr i} {
              add PhiBins [expr {$i * $pi/507}]
          }
  add EtaPhiBins -1.53 $PhiBins
  add EtaPhiBins 1.54 $PhiBins
  set PhiBins {}
          for {set i -510} {$i <= 510} {incr i} {
              add PhiBins [expr {$i * $pi/510}]
          }
  add EtaPhiBins -1.53 $PhiBins
  add EtaPhiBins 1.53 $PhiBins
  set PhiBins {}
          for {set i -513} {$i <= 513} {incr i} {
              add PhiBins [expr {$i * $pi/513}]
          }
  add EtaPhiBins -1.52 $PhiBins
  add EtaPhiBins 1.53 $PhiBins
  set PhiBins {}
          for {set i -516} {$i <= 516} {incr i} {
              add PhiBins [expr {$i * $pi/516}]
          }
  add EtaPhiBins -1.51 $PhiBins
  add EtaPhiBins 1.52 $PhiBins
  set PhiBins {}
          for {set i -519} {$i <= 519} {incr i} {
              add PhiBins [expr {$i * $pi/519}]
          }
  add EtaPhiBins -1.51 $PhiBins
  add EtaPhiBins 1.51 $PhiBins
  set PhiBins {}
          for {set i -523} {$i <= 523} {incr i} {
              add PhiBins [expr {$i * $pi/523}]
          }
  add EtaPhiBins -1.50 $PhiBins
  add EtaPhiBins 1.51 $PhiBins
  set PhiBins {}
          for {set i -526} {$i <= 526} {incr i} {
              add PhiBins [expr {$i * $pi/526}]
          }
  add EtaPhiBins -1.50 $PhiBins
  add EtaPhiBins 1.50 $PhiBins
  set PhiBins {}
          for {set i -529} {$i <= 529} {incr i} {
              add PhiBins [expr {$i * $pi/529}]
          }
  add EtaPhiBins -1.50 $PhiBins
  add EtaPhiBins 1.50 $PhiBins

  set EtaPhiRes 0.02
  set EtaMax 1.5

  set nbins_phi [expr {$pi/$EtaPhiRes} ]
  set nbins_phi [expr {int($nbins_phi)} ]

  set PhiBins {}
  for {set i -$nbins_phi} {$i <= $nbins_phi} {incr i} {
    add PhiBins [expr {$i * $pi/$nbins_phi}]
  }

  set nbins_eta [expr {$EtaMax/$EtaPhiRes} ]
  set nbins_eta [expr {int($nbins_eta)} ]

  set nbins_eta_m [expr {int($nbins_eta - 1)} ]

  for {set i -$nbins_eta_m} {$i <= $nbins_eta} {incr i} {
    set eta [expr {$i * $EtaPhiRes}]
    add EtaPhiBins $eta $PhiBins
  }

  # take present CMS granularity for HF

  # 0.175 x (0.175 - 0.35) resolution in eta,phi in the HF 3.0 < |eta| < 5.0
  set PhiBins {}
  for {set i -18} {$i <= 18} {incr i} {
    add PhiBins [expr {$i * $pi/18.0}]
  }

  foreach eta {-5 -4.7 -4.525 -4.35 -4.21 4.35 4.525 4.7 5} {
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

  # for the ECAL barrel (|eta| < 1.5), see hep-ex/1306.2016 and 1502.02701
  # for the endcaps (1.5 < |eta| < 3.0), we take HGCAL  see LHCC-P-008, Fig. 3.39, p.117

  set ResolutionFormula {  (abs(eta) <= 1.50)                    * sqrt(energy^2*0.009^2 + energy*0.12^2 + 0.45^2) +
                           (abs(eta) > 1.50 && abs(eta) <= 1.75) * sqrt(energy^2*0.006^2 + energy*0.20^2) + \
                           (abs(eta) > 1.75 && abs(eta) <= 2.15) * sqrt(energy^2*0.007^2 + energy*0.21^2) + \
                           (abs(eta) > 2.15 && abs(eta) <= 3.00) * sqrt(energy^2*0.008^2 + energy*0.24^2) + \
                           (abs(eta) >= 3.0 && abs(eta) <= 4.176)  * sqrt(energy^2*0.01^2 + energy*0.67^2) + \
                           (abs(eta) >= 4.176 && abs(eta) <= 5.0)  * sqrt(energy^2*0.10^2 + energy*1.82^2)}

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

  set PhiBins {}
          for {set i -139} {$i <= 139} {incr i} {
              add PhiBins [expr {$i * $pi/139}]
          }
  add EtaPhiBins -4.18 $PhiBins
  add EtaPhiBins 4.21 $PhiBins
  set PhiBins {}
          for {set i -142} {$i <= 142} {incr i} {
              add PhiBins [expr {$i * $pi/142}]
          }
  add EtaPhiBins -4.16 $PhiBins
  add EtaPhiBins 4.18 $PhiBins
  set PhiBins {}
          for {set i -145} {$i <= 145} {incr i} {
              add PhiBins [expr {$i * $pi/145}]
          }
  add EtaPhiBins -4.14 $PhiBins
  add EtaPhiBins 4.16 $PhiBins
  set PhiBins {}
          for {set i -148} {$i <= 148} {incr i} {
              add PhiBins [expr {$i * $pi/148}]
          }
  add EtaPhiBins -4.12 $PhiBins
  add EtaPhiBins 4.14 $PhiBins
  set PhiBins {}
          for {set i -151} {$i <= 151} {incr i} {
              add PhiBins [expr {$i * $pi/151}]
          }
  add EtaPhiBins -4.10 $PhiBins
  add EtaPhiBins 4.12 $PhiBins
  set PhiBins {}
          for {set i -155} {$i <= 155} {incr i} {
              add PhiBins [expr {$i * $pi/155}]
          }
  add EtaPhiBins -4.08 $PhiBins
  add EtaPhiBins 4.10 $PhiBins
  set PhiBins {}
          for {set i -158} {$i <= 158} {incr i} {
              add PhiBins [expr {$i * $pi/158}]
          }
  add EtaPhiBins -4.06 $PhiBins
  add EtaPhiBins 4.08 $PhiBins
  set PhiBins {}
          for {set i -161} {$i <= 161} {incr i} {
              add PhiBins [expr {$i * $pi/161}]
          }
  add EtaPhiBins -4.04 $PhiBins
  add EtaPhiBins 4.06 $PhiBins
  set PhiBins {}
          for {set i -164} {$i <= 164} {incr i} {
              add PhiBins [expr {$i * $pi/164}]
          }
  add EtaPhiBins -4.02 $PhiBins
  add EtaPhiBins 4.04 $PhiBins
  set PhiBins {}
          for {set i -167} {$i <= 167} {incr i} {
              add PhiBins [expr {$i * $pi/167}]
          }
  add EtaPhiBins -4.00 $PhiBins
  add EtaPhiBins 4.02 $PhiBins
  set PhiBins {}
          for {set i -170} {$i <= 170} {incr i} {
              add PhiBins [expr {$i * $pi/170}]
          }
  add EtaPhiBins -3.98 $PhiBins
  add EtaPhiBins 4.00 $PhiBins
  set PhiBins {}
          for {set i -173} {$i <= 173} {incr i} {
              add PhiBins [expr {$i * $pi/173}]
          }
  add EtaPhiBins -3.96 $PhiBins
  add EtaPhiBins 3.98 $PhiBins
  set PhiBins {}
          for {set i -176} {$i <= 176} {incr i} {
              add PhiBins [expr {$i * $pi/176}]
          }
  add EtaPhiBins -3.95 $PhiBins
  add EtaPhiBins 3.96 $PhiBins
  set PhiBins {}
          for {set i -180} {$i <= 180} {incr i} {
              add PhiBins [expr {$i * $pi/180}]
          }
  add EtaPhiBins -3.93 $PhiBins
  add EtaPhiBins 3.95 $PhiBins
  set PhiBins {}
          for {set i -183} {$i <= 183} {incr i} {
              add PhiBins [expr {$i * $pi/183}]
          }
  add EtaPhiBins -3.91 $PhiBins
  add EtaPhiBins 3.93 $PhiBins
  set PhiBins {}
          for {set i -186} {$i <= 186} {incr i} {
              add PhiBins [expr {$i * $pi/186}]
          }
  add EtaPhiBins -3.90 $PhiBins
  add EtaPhiBins 3.91 $PhiBins
  set PhiBins {}
          for {set i -189} {$i <= 189} {incr i} {
              add PhiBins [expr {$i * $pi/189}]
          }
  add EtaPhiBins -3.88 $PhiBins
  add EtaPhiBins 3.90 $PhiBins
  set PhiBins {}
          for {set i -192} {$i <= 192} {incr i} {
              add PhiBins [expr {$i * $pi/192}]
          }
  add EtaPhiBins -3.86 $PhiBins
  add EtaPhiBins 3.88 $PhiBins
  set PhiBins {}
          for {set i -195} {$i <= 195} {incr i} {
              add PhiBins [expr {$i * $pi/195}]
          }
  add EtaPhiBins -3.85 $PhiBins
  add EtaPhiBins 3.86 $PhiBins
  set PhiBins {}
          for {set i -198} {$i <= 198} {incr i} {
              add PhiBins [expr {$i * $pi/198}]
          }
  add EtaPhiBins -3.83 $PhiBins
  add EtaPhiBins 3.85 $PhiBins
  set PhiBins {}
          for {set i -202} {$i <= 202} {incr i} {
              add PhiBins [expr {$i * $pi/202}]
          }
  add EtaPhiBins -3.82 $PhiBins
  add EtaPhiBins 3.83 $PhiBins
  set PhiBins {}
          for {set i -205} {$i <= 205} {incr i} {
              add PhiBins [expr {$i * $pi/205}]
          }
  add EtaPhiBins -3.80 $PhiBins
  add EtaPhiBins 3.82 $PhiBins
  set PhiBins {}
          for {set i -208} {$i <= 208} {incr i} {
              add PhiBins [expr {$i * $pi/208}]
          }
  add EtaPhiBins -3.78 $PhiBins
  add EtaPhiBins 3.80 $PhiBins
  set PhiBins {}
          for {set i -211} {$i <= 211} {incr i} {
              add PhiBins [expr {$i * $pi/211}]
          }
  add EtaPhiBins -3.77 $PhiBins
  add EtaPhiBins 3.78 $PhiBins
  set PhiBins {}
          for {set i -214} {$i <= 214} {incr i} {
              add PhiBins [expr {$i * $pi/214}]
          }
  add EtaPhiBins -3.76 $PhiBins
  add EtaPhiBins 3.77 $PhiBins
  set PhiBins {}
          for {set i -217} {$i <= 217} {incr i} {
              add PhiBins [expr {$i * $pi/217}]
          }
  add EtaPhiBins -3.74 $PhiBins
  add EtaPhiBins 3.76 $PhiBins
  set PhiBins {}
          for {set i -220} {$i <= 220} {incr i} {
              add PhiBins [expr {$i * $pi/220}]
          }
  add EtaPhiBins -3.73 $PhiBins
  add EtaPhiBins 3.74 $PhiBins
  set PhiBins {}
          for {set i -224} {$i <= 224} {incr i} {
              add PhiBins [expr {$i * $pi/224}]
          }
  add EtaPhiBins -3.71 $PhiBins
  add EtaPhiBins 3.73 $PhiBins
  set PhiBins {}
          for {set i -227} {$i <= 227} {incr i} {
              add PhiBins [expr {$i * $pi/227}]
          }
  add EtaPhiBins -3.70 $PhiBins
  add EtaPhiBins 3.71 $PhiBins
  set PhiBins {}
          for {set i -230} {$i <= 230} {incr i} {
              add PhiBins [expr {$i * $pi/230}]
          }
  add EtaPhiBins -3.69 $PhiBins
  add EtaPhiBins 3.70 $PhiBins
  set PhiBins {}
          for {set i -233} {$i <= 233} {incr i} {
              add PhiBins [expr {$i * $pi/233}]
          }
  add EtaPhiBins -3.67 $PhiBins
  add EtaPhiBins 3.69 $PhiBins
  set PhiBins {}
          for {set i -236} {$i <= 236} {incr i} {
              add PhiBins [expr {$i * $pi/236}]
          }
  add EtaPhiBins -3.66 $PhiBins
  add EtaPhiBins 3.67 $PhiBins
  set PhiBins {}
          for {set i -239} {$i <= 239} {incr i} {
              add PhiBins [expr {$i * $pi/239}]
          }
  add EtaPhiBins -3.65 $PhiBins
  add EtaPhiBins 3.66 $PhiBins
  set PhiBins {}
          for {set i -242} {$i <= 242} {incr i} {
              add PhiBins [expr {$i * $pi/242}]
          }
  add EtaPhiBins -3.63 $PhiBins
  add EtaPhiBins 3.65 $PhiBins
  set PhiBins {}
          for {set i -246} {$i <= 246} {incr i} {
              add PhiBins [expr {$i * $pi/246}]
          }
  add EtaPhiBins -3.62 $PhiBins
  add EtaPhiBins 3.63 $PhiBins
  set PhiBins {}
          for {set i -249} {$i <= 249} {incr i} {
              add PhiBins [expr {$i * $pi/249}]
          }
  add EtaPhiBins -3.61 $PhiBins
  add EtaPhiBins 3.62 $PhiBins
  set PhiBins {}
          for {set i -252} {$i <= 252} {incr i} {
              add PhiBins [expr {$i * $pi/252}]
          }
  add EtaPhiBins -3.59 $PhiBins
  add EtaPhiBins 3.61 $PhiBins
  set PhiBins {}
          for {set i -255} {$i <= 255} {incr i} {
              add PhiBins [expr {$i * $pi/255}]
          }
  add EtaPhiBins -3.58 $PhiBins
  add EtaPhiBins 3.59 $PhiBins
  set PhiBins {}
          for {set i -258} {$i <= 258} {incr i} {
              add PhiBins [expr {$i * $pi/258}]
          }
  add EtaPhiBins -3.57 $PhiBins
  add EtaPhiBins 3.58 $PhiBins
  set PhiBins {}
          for {set i -261} {$i <= 261} {incr i} {
              add PhiBins [expr {$i * $pi/261}]
          }
  add EtaPhiBins -3.56 $PhiBins
  add EtaPhiBins 3.57 $PhiBins
  set PhiBins {}
          for {set i -264} {$i <= 264} {incr i} {
              add PhiBins [expr {$i * $pi/264}]
          }
  add EtaPhiBins -3.55 $PhiBins
  add EtaPhiBins 3.56 $PhiBins
  set PhiBins {}
          for {set i -268} {$i <= 268} {incr i} {
              add PhiBins [expr {$i * $pi/268}]
          }
  add EtaPhiBins -3.54 $PhiBins
  add EtaPhiBins 3.55 $PhiBins
  set PhiBins {}
          for {set i -271} {$i <= 271} {incr i} {
              add PhiBins [expr {$i * $pi/271}]
          }
  add EtaPhiBins -3.52 $PhiBins
  add EtaPhiBins 3.54 $PhiBins
  set PhiBins {}
          for {set i -274} {$i <= 274} {incr i} {
              add PhiBins [expr {$i * $pi/274}]
          }
  add EtaPhiBins -3.51 $PhiBins
  add EtaPhiBins 3.52 $PhiBins
  set PhiBins {}
          for {set i -277} {$i <= 277} {incr i} {
              add PhiBins [expr {$i * $pi/277}]
          }
  add EtaPhiBins -3.50 $PhiBins
  add EtaPhiBins 3.51 $PhiBins
  set PhiBins {}
          for {set i -280} {$i <= 280} {incr i} {
              add PhiBins [expr {$i * $pi/280}]
          }
  add EtaPhiBins -3.49 $PhiBins
  add EtaPhiBins 3.50 $PhiBins
  set PhiBins {}
          for {set i -283} {$i <= 283} {incr i} {
              add PhiBins [expr {$i * $pi/283}]
          }
  add EtaPhiBins -3.48 $PhiBins
  add EtaPhiBins 3.49 $PhiBins
  set PhiBins {}
          for {set i -286} {$i <= 286} {incr i} {
              add PhiBins [expr {$i * $pi/286}]
          }
  add EtaPhiBins -3.47 $PhiBins
  add EtaPhiBins 3.48 $PhiBins
  set PhiBins {}
          for {set i -290} {$i <= 290} {incr i} {
              add PhiBins [expr {$i * $pi/290}]
          }
  add EtaPhiBins -3.46 $PhiBins
  add EtaPhiBins 3.47 $PhiBins
  set PhiBins {}
          for {set i -293} {$i <= 293} {incr i} {
              add PhiBins [expr {$i * $pi/293}]
          }
  add EtaPhiBins -3.45 $PhiBins
  add EtaPhiBins 3.46 $PhiBins
  set PhiBins {}
          for {set i -296} {$i <= 296} {incr i} {
              add PhiBins [expr {$i * $pi/296}]
          }
  add EtaPhiBins -3.44 $PhiBins
  add EtaPhiBins 3.45 $PhiBins
  set PhiBins {}
          for {set i -299} {$i <= 299} {incr i} {
              add PhiBins [expr {$i * $pi/299}]
          }
  add EtaPhiBins -3.43 $PhiBins
  add EtaPhiBins 3.44 $PhiBins
  set PhiBins {}
          for {set i -302} {$i <= 302} {incr i} {
              add PhiBins [expr {$i * $pi/302}]
          }
  add EtaPhiBins -3.41 $PhiBins
  add EtaPhiBins 3.43 $PhiBins
  set PhiBins {}
          for {set i -305} {$i <= 305} {incr i} {
              add PhiBins [expr {$i * $pi/305}]
          }
  add EtaPhiBins -3.40 $PhiBins
  add EtaPhiBins 3.41 $PhiBins
  set PhiBins {}
          for {set i -308} {$i <= 308} {incr i} {
              add PhiBins [expr {$i * $pi/308}]
          }
  add EtaPhiBins -3.39 $PhiBins
  add EtaPhiBins 3.40 $PhiBins
  set PhiBins {}
          for {set i -312} {$i <= 312} {incr i} {
              add PhiBins [expr {$i * $pi/312}]
          }
  add EtaPhiBins -3.38 $PhiBins
  add EtaPhiBins 3.39 $PhiBins
  set PhiBins {}
          for {set i -315} {$i <= 315} {incr i} {
              add PhiBins [expr {$i * $pi/315}]
          }
  add EtaPhiBins -3.37 $PhiBins
  add EtaPhiBins 3.38 $PhiBins
  set PhiBins {}
          for {set i -318} {$i <= 318} {incr i} {
              add PhiBins [expr {$i * $pi/318}]
          }
  add EtaPhiBins -3.36 $PhiBins
  add EtaPhiBins 3.37 $PhiBins
  set PhiBins {}
          for {set i -321} {$i <= 321} {incr i} {
              add PhiBins [expr {$i * $pi/321}]
          }
  add EtaPhiBins -3.35 $PhiBins
  add EtaPhiBins 3.36 $PhiBins
  set PhiBins {}
          for {set i -324} {$i <= 324} {incr i} {
              add PhiBins [expr {$i * $pi/324}]
          }
  add EtaPhiBins -3.35 $PhiBins
  add EtaPhiBins 3.35 $PhiBins
  set PhiBins {}
          for {set i -327} {$i <= 327} {incr i} {
              add PhiBins [expr {$i * $pi/327}]
          }
  add EtaPhiBins -3.34 $PhiBins
  add EtaPhiBins 3.35 $PhiBins
  set PhiBins {}
          for {set i -330} {$i <= 330} {incr i} {
              add PhiBins [expr {$i * $pi/330}]
          }
  add EtaPhiBins -3.33 $PhiBins
  add EtaPhiBins 3.34 $PhiBins
  set PhiBins {}
          for {set i -334} {$i <= 334} {incr i} {
              add PhiBins [expr {$i * $pi/334}]
          }
  add EtaPhiBins -3.32 $PhiBins
  add EtaPhiBins 3.33 $PhiBins
  set PhiBins {}
          for {set i -337} {$i <= 337} {incr i} {
              add PhiBins [expr {$i * $pi/337}]
          }
  add EtaPhiBins -3.31 $PhiBins
  add EtaPhiBins 3.32 $PhiBins
  set PhiBins {}
          for {set i -340} {$i <= 340} {incr i} {
              add PhiBins [expr {$i * $pi/340}]
          }
  add EtaPhiBins -3.30 $PhiBins
  add EtaPhiBins 3.31 $PhiBins
  set PhiBins {}
          for {set i -343} {$i <= 343} {incr i} {
              add PhiBins [expr {$i * $pi/343}]
          }
  add EtaPhiBins -3.29 $PhiBins
  add EtaPhiBins 3.30 $PhiBins
  set PhiBins {}
          for {set i -346} {$i <= 346} {incr i} {
              add PhiBins [expr {$i * $pi/346}]
          }
  add EtaPhiBins -3.28 $PhiBins
  add EtaPhiBins 3.29 $PhiBins
  set PhiBins {}
          for {set i -349} {$i <= 349} {incr i} {
              add PhiBins [expr {$i * $pi/349}]
          }
  add EtaPhiBins -3.27 $PhiBins
  add EtaPhiBins 3.28 $PhiBins
  set PhiBins {}
          for {set i -352} {$i <= 352} {incr i} {
              add PhiBins [expr {$i * $pi/352}]
          }
  add EtaPhiBins -3.26 $PhiBins
  add EtaPhiBins 3.27 $PhiBins
  set PhiBins {}
          for {set i -356} {$i <= 356} {incr i} {
              add PhiBins [expr {$i * $pi/356}]
          }
  add EtaPhiBins -3.25 $PhiBins
  add EtaPhiBins 3.26 $PhiBins
  set PhiBins {}
          for {set i -359} {$i <= 359} {incr i} {
              add PhiBins [expr {$i * $pi/359}]
          }
  add EtaPhiBins -3.24 $PhiBins
  add EtaPhiBins 3.25 $PhiBins
  set PhiBins {}
          for {set i -362} {$i <= 362} {incr i} {
              add PhiBins [expr {$i * $pi/362}]
          }
  add EtaPhiBins -3.24 $PhiBins
  add EtaPhiBins 3.24 $PhiBins
  set PhiBins {}
          for {set i -365} {$i <= 365} {incr i} {
              add PhiBins [expr {$i * $pi/365}]
          }
  add EtaPhiBins -3.23 $PhiBins
  add EtaPhiBins 3.24 $PhiBins
  set PhiBins {}
          for {set i -368} {$i <= 368} {incr i} {
              add PhiBins [expr {$i * $pi/368}]
          }
  add EtaPhiBins -3.22 $PhiBins
  add EtaPhiBins 3.23 $PhiBins
  set PhiBins {}
          for {set i -371} {$i <= 371} {incr i} {
              add PhiBins [expr {$i * $pi/371}]
          }
  add EtaPhiBins -3.21 $PhiBins
  add EtaPhiBins 3.22 $PhiBins
  set PhiBins {}
          for {set i -374} {$i <= 374} {incr i} {
              add PhiBins [expr {$i * $pi/374}]
          }
  add EtaPhiBins -3.20 $PhiBins
  add EtaPhiBins 3.21 $PhiBins
  set PhiBins {}
          for {set i -378} {$i <= 378} {incr i} {
              add PhiBins [expr {$i * $pi/378}]
          }
  add EtaPhiBins -3.19 $PhiBins
  add EtaPhiBins 3.20 $PhiBins
  set PhiBins {}
          for {set i -381} {$i <= 381} {incr i} {
              add PhiBins [expr {$i * $pi/381}]
          }
  add EtaPhiBins -3.19 $PhiBins
  add EtaPhiBins 3.19 $PhiBins
  set PhiBins {}
          for {set i -384} {$i <= 384} {incr i} {
              add PhiBins [expr {$i * $pi/384}]
          }
  add EtaPhiBins -3.18 $PhiBins
  add EtaPhiBins 3.19 $PhiBins
  set PhiBins {}
          for {set i -387} {$i <= 387} {incr i} {
              add PhiBins [expr {$i * $pi/387}]
          }
  add EtaPhiBins -3.17 $PhiBins
  add EtaPhiBins 3.18 $PhiBins
  set PhiBins {}
          for {set i -390} {$i <= 390} {incr i} {
              add PhiBins [expr {$i * $pi/390}]
          }
  add EtaPhiBins -3.16 $PhiBins
  add EtaPhiBins 3.17 $PhiBins
  set PhiBins {}
          for {set i -393} {$i <= 393} {incr i} {
              add PhiBins [expr {$i * $pi/393}]
          }
  add EtaPhiBins -3.15 $PhiBins
  add EtaPhiBins 3.16 $PhiBins
  set PhiBins {}
          for {set i -396} {$i <= 396} {incr i} {
              add PhiBins [expr {$i * $pi/396}]
          }
  add EtaPhiBins -3.15 $PhiBins
  add EtaPhiBins 3.15 $PhiBins
  set PhiBins {}
          for {set i -400} {$i <= 400} {incr i} {
              add PhiBins [expr {$i * $pi/400}]
          }
  add EtaPhiBins -3.14 $PhiBins
  add EtaPhiBins 3.15 $PhiBins
  set PhiBins {}
          for {set i -403} {$i <= 403} {incr i} {
              add PhiBins [expr {$i * $pi/403}]
          }
  add EtaPhiBins -3.13 $PhiBins
  add EtaPhiBins 3.14 $PhiBins
  set PhiBins {}
          for {set i -406} {$i <= 406} {incr i} {
              add PhiBins [expr {$i * $pi/406}]
          }
  add EtaPhiBins -3.12 $PhiBins
  add EtaPhiBins 3.13 $PhiBins
  set PhiBins {}
          for {set i -409} {$i <= 409} {incr i} {
              add PhiBins [expr {$i * $pi/409}]
          }
  add EtaPhiBins -3.11 $PhiBins
  add EtaPhiBins 3.12 $PhiBins
  set PhiBins {}
          for {set i -412} {$i <= 412} {incr i} {
              add PhiBins [expr {$i * $pi/412}]
          }
  add EtaPhiBins -3.11 $PhiBins
  add EtaPhiBins 3.11 $PhiBins
  set PhiBins {}
          for {set i -415} {$i <= 415} {incr i} {
              add PhiBins [expr {$i * $pi/415}]
          }
  add EtaPhiBins -3.10 $PhiBins
  add EtaPhiBins 3.11 $PhiBins
  set PhiBins {}
          for {set i -418} {$i <= 418} {incr i} {
              add PhiBins [expr {$i * $pi/418}]
          }
  add EtaPhiBins -3.09 $PhiBins
  add EtaPhiBins 3.10 $PhiBins
  set PhiBins {}
          for {set i -422} {$i <= 422} {incr i} {
              add PhiBins [expr {$i * $pi/422}]
          }
  add EtaPhiBins -3.08 $PhiBins
  add EtaPhiBins 3.09 $PhiBins
  set PhiBins {}
          for {set i -425} {$i <= 425} {incr i} {
              add PhiBins [expr {$i * $pi/425}]
          }
  add EtaPhiBins -3.08 $PhiBins
  add EtaPhiBins 3.08 $PhiBins
  set PhiBins {}
          for {set i -428} {$i <= 428} {incr i} {
              add PhiBins [expr {$i * $pi/428}]
          }
  add EtaPhiBins -3.07 $PhiBins
  add EtaPhiBins 3.08 $PhiBins
  set PhiBins {}
          for {set i -431} {$i <= 431} {incr i} {
              add PhiBins [expr {$i * $pi/431}]
          }
  add EtaPhiBins -3.06 $PhiBins
  add EtaPhiBins 3.07 $PhiBins
  set PhiBins {}
          for {set i -434} {$i <= 434} {incr i} {
              add PhiBins [expr {$i * $pi/434}]
          }
  add EtaPhiBins -3.06 $PhiBins
  add EtaPhiBins 3.06 $PhiBins
  set PhiBins {}
          for {set i -437} {$i <= 437} {incr i} {
              add PhiBins [expr {$i * $pi/437}]
          }
  add EtaPhiBins -3.05 $PhiBins
  add EtaPhiBins 3.06 $PhiBins
  set PhiBins {}
          for {set i -440} {$i <= 440} {incr i} {
              add PhiBins [expr {$i * $pi/440}]
          }
  add EtaPhiBins -3.04 $PhiBins
  add EtaPhiBins 3.05 $PhiBins
  set PhiBins {}
          for {set i -444} {$i <= 444} {incr i} {
              add PhiBins [expr {$i * $pi/444}]
          }
  add EtaPhiBins -3.03 $PhiBins
  add EtaPhiBins 3.04 $PhiBins
  set PhiBins {}
          for {set i -447} {$i <= 447} {incr i} {
              add PhiBins [expr {$i * $pi/447}]
          }
  add EtaPhiBins -3.03 $PhiBins
  add EtaPhiBins 3.03 $PhiBins
  set PhiBins {}
          for {set i -450} {$i <= 450} {incr i} {
              add PhiBins [expr {$i * $pi/450}]
          }
  add EtaPhiBins -3.02 $PhiBins
  add EtaPhiBins 3.03 $PhiBins
  set PhiBins {}
          for {set i -453} {$i <= 453} {incr i} {
              add PhiBins [expr {$i * $pi/453}]
          }
  add EtaPhiBins -3.01 $PhiBins
  add EtaPhiBins 3.02 $PhiBins
  set PhiBins {}
          for {set i -456} {$i <= 456} {incr i} {
              add PhiBins [expr {$i * $pi/456}]
          }
  add EtaPhiBins -3.01 $PhiBins
  add EtaPhiBins 3.01 $PhiBins
  set PhiBins {}
          for {set i -459} {$i <= 459} {incr i} {
              add PhiBins [expr {$i * $pi/459}]
          }
  add EtaPhiBins -3.00 $PhiBins
  add EtaPhiBins 3.01 $PhiBins
  set PhiBins {}
          for {set i -462} {$i <= 462} {incr i} {
              add PhiBins [expr {$i * $pi/462}]
          }
  add EtaPhiBins -2.99 $PhiBins
  add EtaPhiBins 3.00 $PhiBins
  set PhiBins {}
          for {set i -466} {$i <= 466} {incr i} {
              add PhiBins [expr {$i * $pi/466}]
          }
  add EtaPhiBins -2.99 $PhiBins
  add EtaPhiBins 2.99 $PhiBins
  set PhiBins {}
          for {set i -161} {$i <= 161} {incr i} {
              add PhiBins [expr {$i * $pi/161}]
          }
  add EtaPhiBins -2.97 $PhiBins
  add EtaPhiBins 2.99 $PhiBins
  set PhiBins {}
          for {set i -164} {$i <= 164} {incr i} {
              add PhiBins [expr {$i * $pi/164}]
          }
  add EtaPhiBins -2.95 $PhiBins
  add EtaPhiBins 2.97 $PhiBins
  set PhiBins {}
          for {set i -167} {$i <= 167} {incr i} {
              add PhiBins [expr {$i * $pi/167}]
          }
  add EtaPhiBins -2.93 $PhiBins
  add EtaPhiBins 2.95 $PhiBins
  set PhiBins {}
          for {set i -170} {$i <= 170} {incr i} {
              add PhiBins [expr {$i * $pi/170}]
          }
  add EtaPhiBins -2.91 $PhiBins
  add EtaPhiBins 2.93 $PhiBins
  set PhiBins {}
          for {set i -174} {$i <= 174} {incr i} {
              add PhiBins [expr {$i * $pi/174}]
          }
  add EtaPhiBins -2.90 $PhiBins
  add EtaPhiBins 2.91 $PhiBins
  set PhiBins {}
          for {set i -177} {$i <= 177} {incr i} {
              add PhiBins [expr {$i * $pi/177}]
          }
  add EtaPhiBins -2.88 $PhiBins
  add EtaPhiBins 2.90 $PhiBins
  set PhiBins {}
          for {set i -180} {$i <= 180} {incr i} {
              add PhiBins [expr {$i * $pi/180}]
          }
  add EtaPhiBins -2.86 $PhiBins
  add EtaPhiBins 2.88 $PhiBins
  set PhiBins {}
          for {set i -183} {$i <= 183} {incr i} {
              add PhiBins [expr {$i * $pi/183}]
          }
  add EtaPhiBins -2.84 $PhiBins
  add EtaPhiBins 2.86 $PhiBins
  set PhiBins {}
          for {set i -186} {$i <= 186} {incr i} {
              add PhiBins [expr {$i * $pi/186}]
          }
  add EtaPhiBins -2.83 $PhiBins
  add EtaPhiBins 2.84 $PhiBins
  set PhiBins {}
          for {set i -189} {$i <= 189} {incr i} {
              add PhiBins [expr {$i * $pi/189}]
          }
  add EtaPhiBins -2.81 $PhiBins
  add EtaPhiBins 2.83 $PhiBins
  set PhiBins {}
          for {set i -192} {$i <= 192} {incr i} {
              add PhiBins [expr {$i * $pi/192}]
          }
  add EtaPhiBins -2.80 $PhiBins
  add EtaPhiBins 2.81 $PhiBins
  set PhiBins {}
          for {set i -196} {$i <= 196} {incr i} {
              add PhiBins [expr {$i * $pi/196}]
          }
  add EtaPhiBins -2.78 $PhiBins
  add EtaPhiBins 2.80 $PhiBins
  set PhiBins {}
          for {set i -199} {$i <= 199} {incr i} {
              add PhiBins [expr {$i * $pi/199}]
          }
  add EtaPhiBins -2.76 $PhiBins
  add EtaPhiBins 2.78 $PhiBins
  set PhiBins {}
          for {set i -202} {$i <= 202} {incr i} {
              add PhiBins [expr {$i * $pi/202}]
          }
  add EtaPhiBins -2.75 $PhiBins
  add EtaPhiBins 2.76 $PhiBins
  set PhiBins {}
          for {set i -205} {$i <= 205} {incr i} {
              add PhiBins [expr {$i * $pi/205}]
          }
  add EtaPhiBins -2.73 $PhiBins
  add EtaPhiBins 2.75 $PhiBins
  set PhiBins {}
          for {set i -208} {$i <= 208} {incr i} {
              add PhiBins [expr {$i * $pi/208}]
          }
  add EtaPhiBins -2.72 $PhiBins
  add EtaPhiBins 2.73 $PhiBins
  set PhiBins {}
          for {set i -211} {$i <= 211} {incr i} {
              add PhiBins [expr {$i * $pi/211}]
          }
  add EtaPhiBins -2.70 $PhiBins
  add EtaPhiBins 2.72 $PhiBins
  set PhiBins {}
          for {set i -214} {$i <= 214} {incr i} {
              add PhiBins [expr {$i * $pi/214}]
          }
  add EtaPhiBins -2.69 $PhiBins
  add EtaPhiBins 2.70 $PhiBins
  set PhiBins {}
          for {set i -218} {$i <= 218} {incr i} {
              add PhiBins [expr {$i * $pi/218}]
          }
  add EtaPhiBins -2.67 $PhiBins
  add EtaPhiBins 2.69 $PhiBins
  set PhiBins {}
          for {set i -221} {$i <= 221} {incr i} {
              add PhiBins [expr {$i * $pi/221}]
          }
  add EtaPhiBins -2.66 $PhiBins
  add EtaPhiBins 2.67 $PhiBins
  set PhiBins {}
          for {set i -224} {$i <= 224} {incr i} {
              add PhiBins [expr {$i * $pi/224}]
          }
  add EtaPhiBins -2.65 $PhiBins
  add EtaPhiBins 2.66 $PhiBins
  set PhiBins {}
          for {set i -227} {$i <= 227} {incr i} {
              add PhiBins [expr {$i * $pi/227}]
          }
  add EtaPhiBins -2.63 $PhiBins
  add EtaPhiBins 2.65 $PhiBins
  set PhiBins {}
          for {set i -230} {$i <= 230} {incr i} {
              add PhiBins [expr {$i * $pi/230}]
          }
  add EtaPhiBins -2.62 $PhiBins
  add EtaPhiBins 2.63 $PhiBins
  set PhiBins {}
          for {set i -233} {$i <= 233} {incr i} {
              add PhiBins [expr {$i * $pi/233}]
          }
  add EtaPhiBins -2.61 $PhiBins
  add EtaPhiBins 2.62 $PhiBins
  set PhiBins {}
          for {set i -236} {$i <= 236} {incr i} {
              add PhiBins [expr {$i * $pi/236}]
          }
  add EtaPhiBins -2.59 $PhiBins
  add EtaPhiBins 2.61 $PhiBins
  set PhiBins {}
          for {set i -240} {$i <= 240} {incr i} {
              add PhiBins [expr {$i * $pi/240}]
          }
  add EtaPhiBins -2.58 $PhiBins
  add EtaPhiBins 2.59 $PhiBins
  set PhiBins {}
          for {set i -243} {$i <= 243} {incr i} {
              add PhiBins [expr {$i * $pi/243}]
          }
  add EtaPhiBins -2.57 $PhiBins
  add EtaPhiBins 2.58 $PhiBins
  set PhiBins {}
          for {set i -246} {$i <= 246} {incr i} {
              add PhiBins [expr {$i * $pi/246}]
          }
  add EtaPhiBins -2.56 $PhiBins
  add EtaPhiBins 2.57 $PhiBins
  set PhiBins {}
          for {set i -249} {$i <= 249} {incr i} {
              add PhiBins [expr {$i * $pi/249}]
          }
  add EtaPhiBins -2.54 $PhiBins
  add EtaPhiBins 2.56 $PhiBins
  set PhiBins {}
          for {set i -252} {$i <= 252} {incr i} {
              add PhiBins [expr {$i * $pi/252}]
          }
  add EtaPhiBins -2.53 $PhiBins
  add EtaPhiBins 2.54 $PhiBins
  set PhiBins {}
          for {set i -255} {$i <= 255} {incr i} {
              add PhiBins [expr {$i * $pi/255}]
          }
  add EtaPhiBins -2.52 $PhiBins
  add EtaPhiBins 2.53 $PhiBins
  set PhiBins {}
          for {set i -258} {$i <= 258} {incr i} {
              add PhiBins [expr {$i * $pi/258}]
          }
  add EtaPhiBins -2.51 $PhiBins
  add EtaPhiBins 2.52 $PhiBins
  set PhiBins {}
          for {set i -262} {$i <= 262} {incr i} {
              add PhiBins [expr {$i * $pi/262}]
          }
  add EtaPhiBins -2.49 $PhiBins
  add EtaPhiBins 2.51 $PhiBins
  set PhiBins {}
          for {set i -265} {$i <= 265} {incr i} {
              add PhiBins [expr {$i * $pi/265}]
          }
  add EtaPhiBins -2.48 $PhiBins
  add EtaPhiBins 2.49 $PhiBins
  set PhiBins {}
          for {set i -268} {$i <= 268} {incr i} {
              add PhiBins [expr {$i * $pi/268}]
          }
  add EtaPhiBins -2.47 $PhiBins
  add EtaPhiBins 2.48 $PhiBins
  set PhiBins {}
          for {set i -271} {$i <= 271} {incr i} {
              add PhiBins [expr {$i * $pi/271}]
          }
  add EtaPhiBins -2.46 $PhiBins
  add EtaPhiBins 2.47 $PhiBins
  set PhiBins {}
          for {set i -274} {$i <= 274} {incr i} {
              add PhiBins [expr {$i * $pi/274}]
          }
  add EtaPhiBins -2.45 $PhiBins
  add EtaPhiBins 2.46 $PhiBins
  set PhiBins {}
          for {set i -277} {$i <= 277} {incr i} {
              add PhiBins [expr {$i * $pi/277}]
          }
  add EtaPhiBins -2.44 $PhiBins
  add EtaPhiBins 2.45 $PhiBins
  set PhiBins {}
          for {set i -280} {$i <= 280} {incr i} {
              add PhiBins [expr {$i * $pi/280}]
          }
  add EtaPhiBins -2.43 $PhiBins
  add EtaPhiBins 2.44 $PhiBins
  set PhiBins {}
          for {set i -284} {$i <= 284} {incr i} {
              add PhiBins [expr {$i * $pi/284}]
          }
  add EtaPhiBins -2.42 $PhiBins
  add EtaPhiBins 2.43 $PhiBins
  set PhiBins {}
          for {set i -287} {$i <= 287} {incr i} {
              add PhiBins [expr {$i * $pi/287}]
          }
  add EtaPhiBins -2.40 $PhiBins
  add EtaPhiBins 2.42 $PhiBins
  set PhiBins {}
          for {set i -290} {$i <= 290} {incr i} {
              add PhiBins [expr {$i * $pi/290}]
          }
  add EtaPhiBins -2.39 $PhiBins
  add EtaPhiBins 2.40 $PhiBins
  set PhiBins {}
          for {set i -293} {$i <= 293} {incr i} {
              add PhiBins [expr {$i * $pi/293}]
          }
  add EtaPhiBins -2.38 $PhiBins
  add EtaPhiBins 2.39 $PhiBins
  set PhiBins {}
          for {set i -296} {$i <= 296} {incr i} {
              add PhiBins [expr {$i * $pi/296}]
          }
  add EtaPhiBins -2.37 $PhiBins
  add EtaPhiBins 2.38 $PhiBins
  set PhiBins {}
          for {set i -299} {$i <= 299} {incr i} {
              add PhiBins [expr {$i * $pi/299}]
          }
  add EtaPhiBins -2.36 $PhiBins
  add EtaPhiBins 2.37 $PhiBins
  set PhiBins {}
          for {set i -302} {$i <= 302} {incr i} {
              add PhiBins [expr {$i * $pi/302}]
          }
  add EtaPhiBins -2.35 $PhiBins
  add EtaPhiBins 2.36 $PhiBins
  set PhiBins {}
          for {set i -306} {$i <= 306} {incr i} {
              add PhiBins [expr {$i * $pi/306}]
          }
  add EtaPhiBins -2.34 $PhiBins
  add EtaPhiBins 2.35 $PhiBins
  set PhiBins {}
          for {set i -309} {$i <= 309} {incr i} {
              add PhiBins [expr {$i * $pi/309}]
          }
  add EtaPhiBins -2.33 $PhiBins
  add EtaPhiBins 2.34 $PhiBins
  set PhiBins {}
          for {set i -312} {$i <= 312} {incr i} {
              add PhiBins [expr {$i * $pi/312}]
          }
  add EtaPhiBins -2.32 $PhiBins
  add EtaPhiBins 2.33 $PhiBins
  set PhiBins {}
          for {set i -315} {$i <= 315} {incr i} {
              add PhiBins [expr {$i * $pi/315}]
          }
  add EtaPhiBins -2.31 $PhiBins
  add EtaPhiBins 2.32 $PhiBins
  set PhiBins {}
          for {set i -318} {$i <= 318} {incr i} {
              add PhiBins [expr {$i * $pi/318}]
          }
  add EtaPhiBins -2.30 $PhiBins
  add EtaPhiBins 2.31 $PhiBins
  set PhiBins {}
          for {set i -321} {$i <= 321} {incr i} {
              add PhiBins [expr {$i * $pi/321}]
          }
  add EtaPhiBins -2.29 $PhiBins
  add EtaPhiBins 2.30 $PhiBins
  set PhiBins {}
          for {set i -324} {$i <= 324} {incr i} {
              add PhiBins [expr {$i * $pi/324}]
          }
  add EtaPhiBins -2.28 $PhiBins
  add EtaPhiBins 2.29 $PhiBins
  set PhiBins {}
          for {set i -328} {$i <= 328} {incr i} {
              add PhiBins [expr {$i * $pi/328}]
          }
  add EtaPhiBins -2.27 $PhiBins
  add EtaPhiBins 2.28 $PhiBins
  set PhiBins {}
          for {set i -331} {$i <= 331} {incr i} {
              add PhiBins [expr {$i * $pi/331}]
          }
  add EtaPhiBins -2.27 $PhiBins
  add EtaPhiBins 2.27 $PhiBins
  set PhiBins {}
          for {set i -334} {$i <= 334} {incr i} {
              add PhiBins [expr {$i * $pi/334}]
          }
  add EtaPhiBins -2.26 $PhiBins
  add EtaPhiBins 2.27 $PhiBins
  set PhiBins {}
          for {set i -337} {$i <= 337} {incr i} {
              add PhiBins [expr {$i * $pi/337}]
          }
  add EtaPhiBins -2.25 $PhiBins
  add EtaPhiBins 2.26 $PhiBins
  set PhiBins {}
          for {set i -340} {$i <= 340} {incr i} {
              add PhiBins [expr {$i * $pi/340}]
          }
  add EtaPhiBins -2.24 $PhiBins
  add EtaPhiBins 2.25 $PhiBins
  set PhiBins {}
          for {set i -343} {$i <= 343} {incr i} {
              add PhiBins [expr {$i * $pi/343}]
          }
  add EtaPhiBins -2.23 $PhiBins
  add EtaPhiBins 2.24 $PhiBins
  set PhiBins {}
          for {set i -346} {$i <= 346} {incr i} {
              add PhiBins [expr {$i * $pi/346}]
          }
  add EtaPhiBins -2.22 $PhiBins
  add EtaPhiBins 2.23 $PhiBins
  set PhiBins {}
          for {set i -350} {$i <= 350} {incr i} {
              add PhiBins [expr {$i * $pi/350}]
          }
  add EtaPhiBins -2.21 $PhiBins
  add EtaPhiBins 2.22 $PhiBins
  set PhiBins {}
          for {set i -353} {$i <= 353} {incr i} {
              add PhiBins [expr {$i * $pi/353}]
          }
  add EtaPhiBins -2.20 $PhiBins
  add EtaPhiBins 2.21 $PhiBins
  set PhiBins {}
          for {set i -356} {$i <= 356} {incr i} {
              add PhiBins [expr {$i * $pi/356}]
          }
  add EtaPhiBins -2.20 $PhiBins
  add EtaPhiBins 2.20 $PhiBins
  set PhiBins {}
          for {set i -252} {$i <= 252} {incr i} {
              add PhiBins [expr {$i * $pi/252}]
          }
  add EtaPhiBins -2.19 $PhiBins
  add EtaPhiBins 2.20 $PhiBins
  set PhiBins {}
          for {set i -256} {$i <= 256} {incr i} {
              add PhiBins [expr {$i * $pi/256}]
          }
  add EtaPhiBins -2.18 $PhiBins
  add EtaPhiBins 2.19 $PhiBins
  set PhiBins {}
          for {set i -259} {$i <= 259} {incr i} {
              add PhiBins [expr {$i * $pi/259}]
          }
  add EtaPhiBins -2.17 $PhiBins
  add EtaPhiBins 2.18 $PhiBins
  set PhiBins {}
          for {set i -262} {$i <= 262} {incr i} {
              add PhiBins [expr {$i * $pi/262}]
          }
  add EtaPhiBins -2.15 $PhiBins
  add EtaPhiBins 2.17 $PhiBins
  set PhiBins {}
          for {set i -265} {$i <= 265} {incr i} {
              add PhiBins [expr {$i * $pi/265}]
          }
  add EtaPhiBins -2.14 $PhiBins
  add EtaPhiBins 2.15 $PhiBins
  set PhiBins {}
          for {set i -268} {$i <= 268} {incr i} {
              add PhiBins [expr {$i * $pi/268}]
          }
  add EtaPhiBins -2.13 $PhiBins
  add EtaPhiBins 2.14 $PhiBins
  set PhiBins {}
          for {set i -271} {$i <= 271} {incr i} {
              add PhiBins [expr {$i * $pi/271}]
          }
  add EtaPhiBins -2.12 $PhiBins
  add EtaPhiBins 2.13 $PhiBins
  set PhiBins {}
          for {set i -274} {$i <= 274} {incr i} {
              add PhiBins [expr {$i * $pi/274}]
          }
  add EtaPhiBins -2.11 $PhiBins
  add EtaPhiBins 2.12 $PhiBins
  set PhiBins {}
          for {set i -278} {$i <= 278} {incr i} {
              add PhiBins [expr {$i * $pi/278}]
          }
  add EtaPhiBins -2.10 $PhiBins
  add EtaPhiBins 2.11 $PhiBins
  set PhiBins {}
          for {set i -281} {$i <= 281} {incr i} {
              add PhiBins [expr {$i * $pi/281}]
          }
  add EtaPhiBins -2.09 $PhiBins
  add EtaPhiBins 2.10 $PhiBins
  set PhiBins {}
          for {set i -284} {$i <= 284} {incr i} {
              add PhiBins [expr {$i * $pi/284}]
          }
  add EtaPhiBins -2.08 $PhiBins
  add EtaPhiBins 2.09 $PhiBins
  set PhiBins {}
          for {set i -287} {$i <= 287} {incr i} {
              add PhiBins [expr {$i * $pi/287}]
          }
  add EtaPhiBins -2.07 $PhiBins
  add EtaPhiBins 2.08 $PhiBins
  set PhiBins {}
          for {set i -290} {$i <= 290} {incr i} {
              add PhiBins [expr {$i * $pi/290}]
          }
  add EtaPhiBins -2.05 $PhiBins
  add EtaPhiBins 2.07 $PhiBins
  set PhiBins {}
          for {set i -293} {$i <= 293} {incr i} {
              add PhiBins [expr {$i * $pi/293}]
          }
  add EtaPhiBins -2.04 $PhiBins
  add EtaPhiBins 2.05 $PhiBins
  set PhiBins {}
          for {set i -296} {$i <= 296} {incr i} {
              add PhiBins [expr {$i * $pi/296}]
          }
  add EtaPhiBins -2.03 $PhiBins
  add EtaPhiBins 2.04 $PhiBins
  set PhiBins {}
          for {set i -300} {$i <= 300} {incr i} {
              add PhiBins [expr {$i * $pi/300}]
          }
  add EtaPhiBins -2.02 $PhiBins
  add EtaPhiBins 2.03 $PhiBins
  set PhiBins {}
          for {set i -303} {$i <= 303} {incr i} {
              add PhiBins [expr {$i * $pi/303}]
          }
  add EtaPhiBins -2.01 $PhiBins
  add EtaPhiBins 2.02 $PhiBins
  set PhiBins {}
          for {set i -306} {$i <= 306} {incr i} {
              add PhiBins [expr {$i * $pi/306}]
          }
  add EtaPhiBins -2.00 $PhiBins
  add EtaPhiBins 2.01 $PhiBins
  set PhiBins {}
          for {set i -309} {$i <= 309} {incr i} {
              add PhiBins [expr {$i * $pi/309}]
          }
  add EtaPhiBins -1.99 $PhiBins
  add EtaPhiBins 2.00 $PhiBins
  set PhiBins {}
          for {set i -312} {$i <= 312} {incr i} {
              add PhiBins [expr {$i * $pi/312}]
          }
  add EtaPhiBins -1.98 $PhiBins
  add EtaPhiBins 1.99 $PhiBins
  set PhiBins {}
          for {set i -315} {$i <= 315} {incr i} {
              add PhiBins [expr {$i * $pi/315}]
          }
  add EtaPhiBins -1.98 $PhiBins
  add EtaPhiBins 1.98 $PhiBins
  set PhiBins {}
          for {set i -318} {$i <= 318} {incr i} {
              add PhiBins [expr {$i * $pi/318}]
          }
  add EtaPhiBins -1.97 $PhiBins
  add EtaPhiBins 1.98 $PhiBins
  set PhiBins {}
          for {set i -322} {$i <= 322} {incr i} {
              add PhiBins [expr {$i * $pi/322}]
          }
  add EtaPhiBins -1.96 $PhiBins
  add EtaPhiBins 1.97 $PhiBins
  set PhiBins {}
          for {set i -325} {$i <= 325} {incr i} {
              add PhiBins [expr {$i * $pi/325}]
          }
  add EtaPhiBins -1.95 $PhiBins
  add EtaPhiBins 1.96 $PhiBins
  set PhiBins {}
          for {set i -328} {$i <= 328} {incr i} {
              add PhiBins [expr {$i * $pi/328}]
          }
  add EtaPhiBins -1.94 $PhiBins
  add EtaPhiBins 1.95 $PhiBins
  set PhiBins {}
          for {set i -331} {$i <= 331} {incr i} {
              add PhiBins [expr {$i * $pi/331}]
          }
  add EtaPhiBins -1.93 $PhiBins
  add EtaPhiBins 1.94 $PhiBins
  set PhiBins {}
          for {set i -334} {$i <= 334} {incr i} {
              add PhiBins [expr {$i * $pi/334}]
          }
  add EtaPhiBins -1.92 $PhiBins
  add EtaPhiBins 1.93 $PhiBins
  set PhiBins {}
          for {set i -337} {$i <= 337} {incr i} {
              add PhiBins [expr {$i * $pi/337}]
          }
  add EtaPhiBins -1.91 $PhiBins
  add EtaPhiBins 1.92 $PhiBins
  set PhiBins {}
          for {set i -340} {$i <= 340} {incr i} {
              add PhiBins [expr {$i * $pi/340}]
          }
  add EtaPhiBins -1.90 $PhiBins
  add EtaPhiBins 1.91 $PhiBins
  set PhiBins {}
          for {set i -344} {$i <= 344} {incr i} {
              add PhiBins [expr {$i * $pi/344}]
          }
  add EtaPhiBins -1.89 $PhiBins
  add EtaPhiBins 1.90 $PhiBins
  set PhiBins {}
          for {set i -347} {$i <= 347} {incr i} {
              add PhiBins [expr {$i * $pi/347}]
          }
  add EtaPhiBins -1.88 $PhiBins
  add EtaPhiBins 1.89 $PhiBins
  set PhiBins {}
          for {set i -350} {$i <= 350} {incr i} {
              add PhiBins [expr {$i * $pi/350}]
          }
  add EtaPhiBins -1.88 $PhiBins
  add EtaPhiBins 1.88 $PhiBins
  set PhiBins {}
          for {set i -353} {$i <= 353} {incr i} {
              add PhiBins [expr {$i * $pi/353}]
          }
  add EtaPhiBins -1.87 $PhiBins
  add EtaPhiBins 1.88 $PhiBins
  set PhiBins {}
          for {set i -356} {$i <= 356} {incr i} {
              add PhiBins [expr {$i * $pi/356}]
          }
  add EtaPhiBins -1.86 $PhiBins
  add EtaPhiBins 1.87 $PhiBins
  set PhiBins {}
          for {set i -359} {$i <= 359} {incr i} {
              add PhiBins [expr {$i * $pi/359}]
          }
  add EtaPhiBins -1.85 $PhiBins
  add EtaPhiBins 1.86 $PhiBins
  set PhiBins {}
          for {set i -362} {$i <= 362} {incr i} {
              add PhiBins [expr {$i * $pi/362}]
          }
  add EtaPhiBins -1.84 $PhiBins
  add EtaPhiBins 1.85 $PhiBins
  set PhiBins {}
          for {set i -365} {$i <= 365} {incr i} {
              add PhiBins [expr {$i * $pi/365}]
          }
  add EtaPhiBins -1.83 $PhiBins
  add EtaPhiBins 1.84 $PhiBins
  set PhiBins {}
          for {set i -369} {$i <= 369} {incr i} {
              add PhiBins [expr {$i * $pi/369}]
          }
  add EtaPhiBins -1.83 $PhiBins
  add EtaPhiBins 1.83 $PhiBins
  set PhiBins {}
          for {set i -372} {$i <= 372} {incr i} {
              add PhiBins [expr {$i * $pi/372}]
          }
  add EtaPhiBins -1.82 $PhiBins
  add EtaPhiBins 1.83 $PhiBins
  set PhiBins {}
          for {set i -375} {$i <= 375} {incr i} {
              add PhiBins [expr {$i * $pi/375}]
          }
  add EtaPhiBins -1.81 $PhiBins
  add EtaPhiBins 1.82 $PhiBins
  set PhiBins {}
          for {set i -378} {$i <= 378} {incr i} {
              add PhiBins [expr {$i * $pi/378}]
          }
  add EtaPhiBins -1.80 $PhiBins
  add EtaPhiBins 1.81 $PhiBins
  set PhiBins {}
          for {set i -381} {$i <= 381} {incr i} {
              add PhiBins [expr {$i * $pi/381}]
          }
  add EtaPhiBins -1.79 $PhiBins
  add EtaPhiBins 1.80 $PhiBins
  set PhiBins {}
          for {set i -384} {$i <= 384} {incr i} {
              add PhiBins [expr {$i * $pi/384}]
          }
  add EtaPhiBins -1.79 $PhiBins
  add EtaPhiBins 1.79 $PhiBins
  set PhiBins {}
          for {set i -387} {$i <= 387} {incr i} {
              add PhiBins [expr {$i * $pi/387}]
          }
  add EtaPhiBins -1.78 $PhiBins
  add EtaPhiBins 1.79 $PhiBins
  set PhiBins {}
          for {set i -391} {$i <= 391} {incr i} {
              add PhiBins [expr {$i * $pi/391}]
          }
  add EtaPhiBins -1.77 $PhiBins
  add EtaPhiBins 1.78 $PhiBins
  set PhiBins {}
          for {set i -394} {$i <= 394} {incr i} {
              add PhiBins [expr {$i * $pi/394}]
          }
  add EtaPhiBins -1.76 $PhiBins
  add EtaPhiBins 1.77 $PhiBins
  set PhiBins {}
          for {set i -397} {$i <= 397} {incr i} {
              add PhiBins [expr {$i * $pi/397}]
          }
  add EtaPhiBins -1.76 $PhiBins
  add EtaPhiBins 1.76 $PhiBins
  set PhiBins {}
          for {set i -400} {$i <= 400} {incr i} {
              add PhiBins [expr {$i * $pi/400}]
          }
  add EtaPhiBins -1.75 $PhiBins
  add EtaPhiBins 1.76 $PhiBins
  set PhiBins {}
          for {set i -403} {$i <= 403} {incr i} {
              add PhiBins [expr {$i * $pi/403}]
          }
  add EtaPhiBins -1.74 $PhiBins
  add EtaPhiBins 1.75 $PhiBins
  set PhiBins {}
          for {set i -406} {$i <= 406} {incr i} {
              add PhiBins [expr {$i * $pi/406}]
          }
  add EtaPhiBins -1.73 $PhiBins
  add EtaPhiBins 1.74 $PhiBins
  set PhiBins {}
          for {set i -409} {$i <= 409} {incr i} {
              add PhiBins [expr {$i * $pi/409}]
          }
  add EtaPhiBins -1.73 $PhiBins
  add EtaPhiBins 1.73 $PhiBins
  set PhiBins {}
          for {set i -413} {$i <= 413} {incr i} {
              add PhiBins [expr {$i * $pi/413}]
          }
  add EtaPhiBins -1.72 $PhiBins
  add EtaPhiBins 1.73 $PhiBins
  set PhiBins {}
          for {set i -416} {$i <= 416} {incr i} {
              add PhiBins [expr {$i * $pi/416}]
          }
  add EtaPhiBins -1.71 $PhiBins
  add EtaPhiBins 1.72 $PhiBins
  set PhiBins {}
          for {set i -419} {$i <= 419} {incr i} {
              add PhiBins [expr {$i * $pi/419}]
          }
  add EtaPhiBins -1.71 $PhiBins
  add EtaPhiBins 1.71 $PhiBins
  set PhiBins {}
          for {set i -422} {$i <= 422} {incr i} {
              add PhiBins [expr {$i * $pi/422}]
          }
  add EtaPhiBins -1.70 $PhiBins
  add EtaPhiBins 1.71 $PhiBins
  set PhiBins {}
          for {set i -425} {$i <= 425} {incr i} {
              add PhiBins [expr {$i * $pi/425}]
          }
  add EtaPhiBins -1.69 $PhiBins
  add EtaPhiBins 1.70 $PhiBins
  set PhiBins {}
          for {set i -428} {$i <= 428} {incr i} {
              add PhiBins [expr {$i * $pi/428}]
          }
  add EtaPhiBins -1.69 $PhiBins
  add EtaPhiBins 1.69 $PhiBins
  set PhiBins {}
          for {set i -431} {$i <= 431} {incr i} {
              add PhiBins [expr {$i * $pi/431}]
          }
  add EtaPhiBins -1.68 $PhiBins
  add EtaPhiBins 1.69 $PhiBins
  set PhiBins {}
          for {set i -435} {$i <= 435} {incr i} {
              add PhiBins [expr {$i * $pi/435}]
          }
  add EtaPhiBins -1.67 $PhiBins
  add EtaPhiBins 1.68 $PhiBins
  set PhiBins {}
          for {set i -438} {$i <= 438} {incr i} {
              add PhiBins [expr {$i * $pi/438}]
          }
  add EtaPhiBins -1.67 $PhiBins
  add EtaPhiBins 1.67 $PhiBins
  set PhiBins {}
          for {set i -441} {$i <= 441} {incr i} {
              add PhiBins [expr {$i * $pi/441}]
          }
  add EtaPhiBins -1.66 $PhiBins
  add EtaPhiBins 1.67 $PhiBins
  set PhiBins {}
          for {set i -444} {$i <= 444} {incr i} {
              add PhiBins [expr {$i * $pi/444}]
          }
  add EtaPhiBins -1.65 $PhiBins
  add EtaPhiBins 1.66 $PhiBins
  set PhiBins {}
          for {set i -447} {$i <= 447} {incr i} {
              add PhiBins [expr {$i * $pi/447}]
          }
  add EtaPhiBins -1.65 $PhiBins
  add EtaPhiBins 1.65 $PhiBins
  set PhiBins {}
          for {set i -450} {$i <= 450} {incr i} {
              add PhiBins [expr {$i * $pi/450}]
          }
  add EtaPhiBins -1.64 $PhiBins
  add EtaPhiBins 1.65 $PhiBins
  set PhiBins {}
          for {set i -453} {$i <= 453} {incr i} {
              add PhiBins [expr {$i * $pi/453}]
          }
  add EtaPhiBins -1.63 $PhiBins
  add EtaPhiBins 1.64 $PhiBins
  set PhiBins {}
          for {set i -457} {$i <= 457} {incr i} {
              add PhiBins [expr {$i * $pi/457}]
          }
  add EtaPhiBins -1.63 $PhiBins
  add EtaPhiBins 1.63 $PhiBins
  set PhiBins {}
          for {set i -460} {$i <= 460} {incr i} {
              add PhiBins [expr {$i * $pi/460}]
          }
  add EtaPhiBins -1.62 $PhiBins
  add EtaPhiBins 1.63 $PhiBins
  set PhiBins {}
          for {set i -463} {$i <= 463} {incr i} {
              add PhiBins [expr {$i * $pi/463}]
          }
  add EtaPhiBins -1.61 $PhiBins
  add EtaPhiBins 1.62 $PhiBins
  set PhiBins {}
          for {set i -466} {$i <= 466} {incr i} {
              add PhiBins [expr {$i * $pi/466}]
          }
  add EtaPhiBins -1.61 $PhiBins
  add EtaPhiBins 1.61 $PhiBins
  set PhiBins {}
          for {set i -469} {$i <= 469} {incr i} {
              add PhiBins [expr {$i * $pi/469}]
          }
  add EtaPhiBins -1.60 $PhiBins
  add EtaPhiBins 1.61 $PhiBins
  set PhiBins {}
          for {set i -472} {$i <= 472} {incr i} {
              add PhiBins [expr {$i * $pi/472}]
          }
  add EtaPhiBins -1.60 $PhiBins
  add EtaPhiBins 1.60 $PhiBins
  set PhiBins {}
          for {set i -475} {$i <= 475} {incr i} {
              add PhiBins [expr {$i * $pi/475}]
          }
  add EtaPhiBins -1.59 $PhiBins
  add EtaPhiBins 1.60 $PhiBins
  set PhiBins {}
          for {set i -479} {$i <= 479} {incr i} {
              add PhiBins [expr {$i * $pi/479}]
          }
  add EtaPhiBins -1.58 $PhiBins
  add EtaPhiBins 1.59 $PhiBins
  set PhiBins {}
          for {set i -482} {$i <= 482} {incr i} {
              add PhiBins [expr {$i * $pi/482}]
          }
  add EtaPhiBins -1.58 $PhiBins
  add EtaPhiBins 1.58 $PhiBins
  set PhiBins {}
          for {set i -485} {$i <= 485} {incr i} {
              add PhiBins [expr {$i * $pi/485}]
          }
  add EtaPhiBins -1.57 $PhiBins
  add EtaPhiBins 1.58 $PhiBins
  set PhiBins {}
          for {set i -488} {$i <= 488} {incr i} {
              add PhiBins [expr {$i * $pi/488}]
          }
  add EtaPhiBins -1.57 $PhiBins
  add EtaPhiBins 1.57 $PhiBins
  set PhiBins {}
          for {set i -491} {$i <= 491} {incr i} {
              add PhiBins [expr {$i * $pi/491}]
          }
  add EtaPhiBins -1.56 $PhiBins
  add EtaPhiBins 1.57 $PhiBins
  set PhiBins {}
          for {set i -494} {$i <= 494} {incr i} {
              add PhiBins [expr {$i * $pi/494}]
          }
  add EtaPhiBins -1.55 $PhiBins
  add EtaPhiBins 1.56 $PhiBins
  set PhiBins {}
          for {set i -497} {$i <= 497} {incr i} {
              add PhiBins [expr {$i * $pi/497}]
          }
  add EtaPhiBins -1.55 $PhiBins
  add EtaPhiBins 1.55 $PhiBins
  set PhiBins {}
          for {set i -501} {$i <= 501} {incr i} {
              add PhiBins [expr {$i * $pi/501}]
          }
  add EtaPhiBins -1.54 $PhiBins
  add EtaPhiBins 1.55 $PhiBins
  set PhiBins {}
          for {set i -504} {$i <= 504} {incr i} {
              add PhiBins [expr {$i * $pi/504}]
          }
  add EtaPhiBins -1.54 $PhiBins
  add EtaPhiBins 1.54 $PhiBins
  set PhiBins {}
          for {set i -507} {$i <= 507} {incr i} {
              add PhiBins [expr {$i * $pi/507}]
          }
  add EtaPhiBins -1.53 $PhiBins
  add EtaPhiBins 1.54 $PhiBins
  set PhiBins {}
          for {set i -510} {$i <= 510} {incr i} {
              add PhiBins [expr {$i * $pi/510}]
          }
  add EtaPhiBins -1.53 $PhiBins
  add EtaPhiBins 1.53 $PhiBins
  set PhiBins {}
          for {set i -513} {$i <= 513} {incr i} {
              add PhiBins [expr {$i * $pi/513}]
          }
  add EtaPhiBins -1.52 $PhiBins
  add EtaPhiBins 1.53 $PhiBins
  set PhiBins {}
          for {set i -516} {$i <= 516} {incr i} {
              add PhiBins [expr {$i * $pi/516}]
          }
  add EtaPhiBins -1.51 $PhiBins
  add EtaPhiBins 1.52 $PhiBins
  set PhiBins {}
          for {set i -519} {$i <= 519} {incr i} {
              add PhiBins [expr {$i * $pi/519}]
          }
  add EtaPhiBins -1.51 $PhiBins
  add EtaPhiBins 1.51 $PhiBins
  set PhiBins {}
          for {set i -523} {$i <= 523} {incr i} {
              add PhiBins [expr {$i * $pi/523}]
          }
  add EtaPhiBins -1.50 $PhiBins
  add EtaPhiBins 1.51 $PhiBins
  set PhiBins {}
          for {set i -526} {$i <= 526} {incr i} {
              add PhiBins [expr {$i * $pi/526}]
          }
  add EtaPhiBins -1.50 $PhiBins
  add EtaPhiBins 1.50 $PhiBins
  set PhiBins {}
          for {set i -529} {$i <= 529} {incr i} {
              add PhiBins [expr {$i * $pi/529}]
          }
  add EtaPhiBins -1.50 $PhiBins
  add EtaPhiBins 1.50 $PhiBins

  set EtaPhiRes 0.1
  set EtaMax 1.5

  set nbins_phi [expr {$pi/$EtaPhiRes} ]
  set nbins_phi [expr {int($nbins_phi)} ]

  set PhiBins {}
  for {set i -$nbins_phi} {$i <= $nbins_phi} {incr i} {
    add PhiBins [expr {$i * $pi/$nbins_phi}]
  }

  set nbins_eta [expr {$EtaMax/$EtaPhiRes} ]
  set nbins_eta [expr {int($nbins_eta)} ]

  set nbins_eta_m [expr {int($nbins_eta - 1)} ]

  for {set i -$nbins_eta_m} {$i <= $nbins_eta} {incr i} {
    set eta [expr {$i * $EtaPhiRes}]
    add EtaPhiBins $eta $PhiBins
  }


  # take present CMS granularity for HF

  # 0.175 x (0.175 - 0.35) resolution in eta,phi in the HF 3.0 < |eta| < 5.0
  set PhiBins {}
  for {set i -18} {$i <= 18} {incr i} {
    add PhiBins [expr {$i * $pi/18.0}]
  }

  foreach eta {-5 -4.7 -4.525 -4.35 -4.21 4.35 4.525 4.7 5} {
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
  set ResolutionFormula {                    (abs(eta) <= 1.5) * sqrt(energy^2*0.05^2 + energy*1.00^2) + \
                                                   (abs(eta) > 1.5 && abs(eta) <= 3.0) * sqrt(energy^2*0.05^2 + energy*1.00^2) + \
                                                   (abs(eta) > 3.0 && abs(eta) <= 4.2108) * sqrt(energy^2*0.05^2 + energy*1.96^2) + \
                                                   (abs(eta) > 4.2108 && abs(eta) <= 5.0) * sqrt(energy^2*0.11^2 + energy*2.40^2)}

}

#################################
# Energy resolution for electrons
#################################

module EnergySmearing PhotonEnergySmearing {
  set InputArray ECal/eflowPhotons
  set OutputArray eflowPhotons

  # adding 1% extra photon smearing
  set ResolutionFormula {energy*0.01}

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
  add InputArray PhotonEnergySmearing/eflowPhotons
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
  set InputArray PhotonEnergySmearing/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for photons
  set EfficiencyFormula {                      (pt <= 10.0) * (0.00) + \
                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.9635) + \
         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 10.0)  * (0.9624) + \
         (abs(eta) > 4.0)                                   * (0.00)}
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
    set EfficiencyFormula {
                          (pt <= 4.0) * (0.00) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 4.0 && pt <= 6.0) * (0.018) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 6.0 && pt <= 8.0) * (0.252) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 8.0 && pt <= 10.0) * (0.480) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 10.0 && pt <= 20.0) * (0.681) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 20.0 && pt <= 35.0) * (0.792) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 35.0 && pt <= 50.0) * (0.862) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 50.0 && pt <= 14000.0) * (0.859) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 4.0 && pt <= 6.0) * (0.016) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 6.0 && pt <= 8.0) * (0.198) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 8.0 && pt <= 10.0) * (0.446) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 10.0 && pt <= 20.0) * (0.598) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 20.0 && pt <= 35.0) * (0.759) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 35.0 && pt <= 50.0) * (0.847) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 50.0 && pt <= 14000.0) * (0.872) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.45) * (pt > 4.0 && pt <= 6.0) * (0.005) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.45) * (pt > 6.0 && pt <= 8.0) * (0.029) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.45) * (pt > 8.0 && pt <= 10.0) * (0.108) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.45) * (pt > 10.0 && pt <= 20.0) * (0.289) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.45) * (pt > 20.0 && pt <= 35.0) * (0.570) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.45) * (pt > 35.0 && pt <= 50.0) * (0.743) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.45) * (pt > 50.0 && pt <= 14000.0) * (0.828) +
                          (abs(eta) > 1.45 && abs(eta) <= 1.55) * (pt > 4.0 && pt <= 6.0) * (0.026) +
                          (abs(eta) > 1.45 && abs(eta) <= 1.55) * (pt > 6.0 && pt <= 8.0) * (0.045) +
                          (abs(eta) > 1.45 && abs(eta) <= 1.55) * (pt > 8.0 && pt <= 10.0) * (0.133) +
                          (abs(eta) > 1.45 && abs(eta) <= 1.55) * (pt > 10.0 && pt <= 20.0) * (0.411) +
                          (abs(eta) > 1.45 && abs(eta) <= 1.55) * (pt > 20.0 && pt <= 35.0) * (0.629) +
                          (abs(eta) > 1.45 && abs(eta) <= 1.55) * (pt > 35.0 && pt <= 50.0) * (0.761) +
                          (abs(eta) > 1.45 && abs(eta) <= 1.55) * (pt > 50.0 && pt <= 14000.0) * (0.752) +
                          (abs(eta) > 1.55 && abs(eta) <= 2.0) * (pt > 4.0 && pt <= 6.0) * (0.061) +
                          (abs(eta) > 1.55 && abs(eta) <= 2.0) * (pt > 6.0 && pt <= 8.0) * (0.191) +
                          (abs(eta) > 1.55 && abs(eta) <= 2.0) * (pt > 8.0 && pt <= 10.0) * (0.337) +
                          (abs(eta) > 1.55 && abs(eta) <= 2.0) * (pt > 10.0 && pt <= 20.0) * (0.475) +
                          (abs(eta) > 1.55 && abs(eta) <= 2.0) * (pt > 20.0 && pt <= 35.0) * (0.605) +
                          (abs(eta) > 1.55 && abs(eta) <= 2.0) * (pt > 35.0 && pt <= 50.0) * (0.713) +
                          (abs(eta) > 1.55 && abs(eta) <= 2.0) * (pt > 50.0 && pt <= 14000.0) * (0.794) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.5) * (pt > 4.0 && pt <= 6.0) * (0.100) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.5) * (pt > 6.0 && pt <= 8.0) * (0.223) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.5) * (pt > 8.0 && pt <= 10.0) * (0.427) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.5) * (pt > 10.0 && pt <= 20.0) * (0.590) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.5) * (pt > 20.0 && pt <= 35.0) * (0.720) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.5) * (pt > 35.0 && pt <= 50.0) * (0.800) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.5) * (pt > 50.0 && pt <= 14000.0) * (0.840) +
                          (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 4.0 && pt <= 6.0) * (0.049) +
                          (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 6.0 && pt <= 8.0) * (0.152) +
                          (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 8.0 && pt <= 10.0) * (0.436) +
                          (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 10.0 && pt <= 20.0) * (0.679) +
                          (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 20.0 && pt <= 35.0) * (0.778) +
                          (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 35.0 && pt <= 50.0) * (0.830) +
                          (abs(eta) > 2.5 && abs(eta) <= 3.0) * (pt > 50.0 && pt <= 14000.0) * (0.919) +

                          (abs(eta) > 3.0 && abs(eta) <= 4.0) * (pt > 4.0 && pt <= 6.0) * (0.049) +
                          (abs(eta) > 3.0  && abs(eta) <= 4.0) * (pt > 6.0 && pt <= 8.0) * (0.152) +
                          (abs(eta) > 3.0  && abs(eta) <= 4.0) * (pt > 8.0 && pt <= 10.0) * (0.436) +
                          (abs(eta) > 3.0  && abs(eta) <= 4.0) * (pt > 10.0 && pt <= 20.0) * (0.679) +
                          (abs(eta) > 3.0  && abs(eta) <= 4.0) * (pt > 20.0 && pt <= 35.0) * (0.778) +
                          (abs(eta) > 3.0  && abs(eta) <= 4.0) * (pt > 35.0 && pt <= 50.0) * (0.830) +
                          (abs(eta) > 3.0  && abs(eta) <= 4.0) * (pt > 50.0 && pt <= 14000.0) * (0.919)}
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

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
  set InputArray MuonMomentumSmearing/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency as a function of eta and pt}

  # efficiency formula for muons
    set EfficiencyFormula {
                          (pt <= 2.0) * (0.00) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 2.0 && pt <= 4.0) * (0.04) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 4.0 && pt <= 6.0) * (0.43) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 6.0 && pt <= 8.0) * (0.53) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 8.0 && pt <= 10.0) * (0.67) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 10.0 && pt <= 20.0) * (0.81) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 20.0 && pt <= 35.0) * (0.90) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 35.0 && pt <= 50.0) * (0.92) +
                          (abs(eta) > 0.0 && abs(eta) <= 0.5) * (pt > 50.0 && pt <= 14000.0) * (0.90) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 2.0 && pt <= 4.0) * (0.05) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 4.0 && pt <= 6.0) * (0.46) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 6.0 && pt <= 8.0) * (0.56) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 8.0 && pt <= 10.0) * (0.65) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 10.0 && pt <= 20.0) * (0.79) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 20.0 && pt <= 35.0) * (0.91) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 35.0 && pt <= 50.0) * (0.93) +
                          (abs(eta) > 0.5 && abs(eta) <= 1.0) * (pt > 50.0 && pt <= 14000.0) * (0.92) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.5) * (pt > 2.0 && pt <= 4.0) * (0.15) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.5) * (pt > 4.0 && pt <= 6.0) * (0.47) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.5) * (pt > 6.0 && pt <= 8.0) * (0.55) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.5) * (pt > 8.0 && pt <= 10.0) * (0.64) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.5) * (pt > 10.0 && pt <= 20.0) * (0.78) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.5) * (pt > 20.0 && pt <= 35.0) * (0.89) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.5) * (pt > 35.0 && pt <= 50.0) * (0.96) +
                          (abs(eta) > 1.0 && abs(eta) <= 1.5) * (pt > 50.0 && pt <= 14000.0) * (0.91) +
                          (abs(eta) > 1.5 && abs(eta) <= 2.0) * (pt > 2.0 && pt <= 4.0) * (0.23) +
                          (abs(eta) > 1.5 && abs(eta) <= 2.0) * (pt > 4.0 && pt <= 6.0) * (0.44) +
                          (abs(eta) > 1.5 && abs(eta) <= 2.0) * (pt > 6.0 && pt <= 8.0) * (0.53) +
                          (abs(eta) > 1.5 && abs(eta) <= 2.0) * (pt > 8.0 && pt <= 10.0) * (0.68) +
                          (abs(eta) > 1.5 && abs(eta) <= 2.0) * (pt > 10.0 && pt <= 20.0) * (0.78) +
                          (abs(eta) > 1.5 && abs(eta) <= 2.0) * (pt > 20.0 && pt <= 35.0) * (0.89) +
                          (abs(eta) > 1.5 && abs(eta) <= 2.0) * (pt > 35.0 && pt <= 50.0) * (0.95) +
                          (abs(eta) > 1.5 && abs(eta) <= 2.0) * (pt > 50.0 && pt <= 14000.0) * (0.88) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.8) * (pt > 2.0 && pt <= 4.0) * (0.22) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.8) * (pt > 4.0 && pt <= 6.0) * (0.36) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.8) * (pt > 6.0 && pt <= 8.0) * (0.44) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.8) * (pt > 8.0 && pt <= 10.0) * (0.57) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.8) * (pt > 10.0 && pt <= 20.0) * (0.63) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.8) * (pt > 20.0 && pt <= 35.0) * (0.71) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.8) * (pt > 35.0 && pt <= 50.0) * (0.76) +
                          (abs(eta) > 2.0 && abs(eta) <= 2.8) * (pt > 50.0 && pt <= 14000.0) * (0.82) +
                          (abs(eta) > 2.8) * (0.00)
      }}

################
# Muon isolation
################

module Isolation MuonIsolation {
  set CandidateInputArray MuonEfficiency/muons
  set IsolationInputArray EFlowFilter/eflow

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
  add InputArray UniqueObjectFinder/jets
  add InputArray UniqueObjectFinder/electrons
  add InputArray UniqueObjectFinder/photons
  add InputArray UniqueObjectFinder/muons
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
  set ParameterR 0.5

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
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.5

  set JetPTMin 20.0
}

##################
# Fat Jet finder
##################

module FastJetFinder FatJetFinder {
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 0.8

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeTrimming 1
  set RTrim 0.2
  set PtFracTrim 0.05

  set ComputePruning 1
  set ZcutPrun 0.1
  set RcutPrun 0.5
  set RPrun 0.8

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

  set JetPTMin 200.0
}




##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder/jets
  set OutputArray jets

  # scale formula for jets
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
  set PartonPTMin 1.0
  set PartonEtaMax 2.5

}

###########
# b-tagging
###########

module BTagging BTagging {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 0

  add EfficiencyFormula {0}      {0.001}

  add EfficiencyFormula {5}      {
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 20.00 && pt <= 30.00) * (0.527) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 30.00 && pt <= 40.00) * (0.598) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 40.00 && pt <= 50.00) * (0.632) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 50.00 && pt <= 60.00) * (0.647) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 60.00 && pt <= 70.00) * (0.652) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 70.00 && pt <= 80.00) * (0.653) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 80.00 && pt <= 90.00) * (0.653) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 90.00 && pt <= 100.00) * (0.655) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 100.00 && pt <= 120.00) * (0.635) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 120.00 && pt <= 140.00) * (0.626) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 140.00 && pt <= 160.00) * (0.614) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 160.00 && pt <= 180.00) * (0.591) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 180.00 && pt <= 200.00) * (0.576) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 200.00 && pt <= 250.00) * (0.550) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 250.00 && pt <= 300.00) * (0.489) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 300.00 && pt <= 350.00) * (0.444) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 350.00 && pt <= 400.00) * (0.408) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 400.00 && pt <= 500.00) * (0.369) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 500.00 && pt <= 600.00) * (0.334) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 600.00 && pt <= 700.00) * (0.269) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 700.00 && pt <= 800.00) * (0.253) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 800.00 && pt <= 1000.00) * (0.247) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 1000.00 && pt <= 1400.00) * (0.230) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 1400.00 && pt <= 2000.00) * (0.209) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 2000.00 && pt <= 3000.00) * (0.209) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 3000.00) * (0.209) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 20.00 && pt <= 30.00) * (0.363) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 30.00 && pt <= 40.00) * (0.465) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 40.00 && pt <= 50.00) * (0.504) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 50.00 && pt <= 60.00) * (0.525) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 60.00 && pt <= 70.00) * (0.531) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 70.00 && pt <= 80.00) * (0.526) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 80.00 && pt <= 90.00) * (0.511) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 90.00 && pt <= 100.00) * (0.506) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 100.00 && pt <= 120.00) * (0.492) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 120.00 && pt <= 140.00) * (0.474) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 140.00 && pt <= 160.00) * (0.452) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 160.00 && pt <= 180.00) * (0.443) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 180.00 && pt <= 200.00) * (0.433) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 200.00 && pt <= 250.00) * (0.385) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 250.00 && pt <= 300.00) * (0.349) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 300.00 && pt <= 350.00) * (0.325) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 350.00 && pt <= 400.00) * (0.316) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 400.00 && pt <= 500.00) * (0.280) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 500.00 && pt <= 600.00) * (0.206) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 600.00 && pt <= 700.00) * (0.206) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 700.00 && pt <= 800.00) * (0.206) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 800.00 && pt <= 1000.00) * (0.206) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 1000.00 && pt <= 1400.00) * (0.206) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 1400.00 && pt <= 2000.00) * (0.206) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 2000.00 && pt <= 3000.00) * (0.206) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 3000.00) * (0.206) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 20.00 && pt <= 30.00) * (0.211) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 30.00 && pt <= 40.00) * (0.278) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 40.00 && pt <= 50.00) * (0.310) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 50.00 && pt <= 60.00) * (0.319) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 60.00 && pt <= 70.00) * (0.327) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 70.00 && pt <= 80.00) * (0.320) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 80.00 && pt <= 90.00) * (0.318) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 90.00 && pt <= 100.00) * (0.314) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 100.00 && pt <= 120.00) * (0.304) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 120.00 && pt <= 140.00) * (0.296) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 140.00 && pt <= 160.00) * (0.271) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 160.00 && pt <= 180.00) * (0.278) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 180.00 && pt <= 200.00) * (0.277) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 200.00 && pt <= 250.00) * (0.255) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 250.00 && pt <= 300.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 300.00 && pt <= 350.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 350.00 && pt <= 400.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 400.00 && pt <= 500.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 500.00 && pt <= 600.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 600.00 && pt <= 700.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 700.00 && pt <= 800.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 800.00 && pt <= 1000.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 1000.00 && pt <= 1400.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 1400.00 && pt <= 2000.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 2000.00 && pt <= 3000.00) * (0.218) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 3000.00) * (0.218)
                                  }

  add EfficiencyFormula {4}      {
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 20.00 && pt <= 30.00) * (0.023) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 30.00 && pt <= 40.00) * (0.025) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 40.00 && pt <= 50.00) * (0.029) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 50.00 && pt <= 60.00) * (0.031) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 60.00 && pt <= 70.00) * (0.034) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 70.00 && pt <= 80.00) * (0.035) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 80.00 && pt <= 90.00) * (0.036) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 90.00 && pt <= 100.00) * (0.038) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 100.00 && pt <= 120.00) * (0.039) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 120.00 && pt <= 140.00) * (0.039) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 140.00 && pt <= 160.00) * (0.042) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 160.00 && pt <= 180.00) * (0.042) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 180.00 && pt <= 200.00) * (0.045) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 200.00 && pt <= 250.00) * (0.040) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 250.00 && pt <= 300.00) * (0.038) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 300.00 && pt <= 350.00) * (0.034) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 350.00 && pt <= 400.00) * (0.031) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 400.00 && pt <= 500.00) * (0.026) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 500.00 && pt <= 600.00) * (0.025) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 600.00 && pt <= 700.00) * (0.019) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 700.00 && pt <= 800.00) * (0.017) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 800.00 && pt <= 1000.00) * (0.018) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 1000.00 && pt <= 1400.00) * (0.016) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 1400.00 && pt <= 2000.00) * (0.015) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 2000.00 && pt <= 3000.00) * (0.015) +
                                  (abs(eta) > 0.00 && abs(eta) <= 1.50) * (pt > 3000.00) * (0.015) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 20.00 && pt <= 30.00) * (0.017) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 30.00 && pt <= 40.00) * (0.018) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 40.00 && pt <= 50.00) * (0.023) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 50.00 && pt <= 60.00) * (0.026) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 60.00 && pt <= 70.00) * (0.029) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 70.00 && pt <= 80.00) * (0.028) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 80.00 && pt <= 90.00) * (0.026) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 90.00 && pt <= 100.00) * (0.032) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 100.00 && pt <= 120.00) * (0.032) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 120.00 && pt <= 140.00) * (0.033) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 140.00 && pt <= 160.00) * (0.034) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 160.00 && pt <= 180.00) * (0.036) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 180.00 && pt <= 200.00) * (0.035) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 200.00 && pt <= 250.00) * (0.035) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 250.00 && pt <= 300.00) * (0.032) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 300.00 && pt <= 350.00) * (0.025) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 350.00 && pt <= 400.00) * (0.027) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 400.00 && pt <= 500.00) * (0.020) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 500.00 && pt <= 600.00) * (0.017) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 600.00 && pt <= 700.00) * (0.017) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 700.00 && pt <= 800.00) * (0.017) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 800.00 && pt <= 1000.00) * (0.017) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 1000.00 && pt <= 1400.00) * (0.017) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 1400.00 && pt <= 2000.00) * (0.017) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 2000.00 && pt <= 3000.00) * (0.017) +
                                  (abs(eta) > 1.50 && abs(eta) <= 2.50) * (pt > 3000.00) * (0.017) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 20.00 && pt <= 30.00) * (0.018) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 30.00 && pt <= 40.00) * (0.017) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 40.00 && pt <= 50.00) * (0.018) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 50.00 && pt <= 60.00) * (0.018) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 60.00 && pt <= 70.00) * (0.020) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 70.00 && pt <= 80.00) * (0.020) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 80.00 && pt <= 90.00) * (0.019) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 90.00 && pt <= 100.00) * (0.020) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 100.00 && pt <= 120.00) * (0.019) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 120.00 && pt <= 140.00) * (0.018) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 140.00 && pt <= 160.00) * (0.021) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 160.00 && pt <= 180.00) * (0.016) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 180.00 && pt <= 200.00) * (0.022) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 200.00 && pt <= 250.00) * (0.023) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 250.00 && pt <= 300.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 300.00 && pt <= 350.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 350.00 && pt <= 400.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 400.00 && pt <= 500.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 500.00 && pt <= 600.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 600.00 && pt <= 700.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 700.00 && pt <= 800.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 800.00 && pt <= 1000.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 1000.00 && pt <= 1400.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 1400.00 && pt <= 2000.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 2000.00 && pt <= 3000.00) * (0.014) +
                                  (abs(eta) > 2.50 && abs(eta) <= 3.50) * (pt > 3000.00) * (0.014)
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

###############################################################################################################
# StatusPidFilter: this module removes all generated particles except electrons, muons, taus, and status == 3 #
###############################################################################################################

module StatusPidFilter GenParticleFilter {

    set InputArray Delphes/allParticles
    set OutputArray filteredParticles
    set PTMin 0.0

}


##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
 # add Branch Delphes/allParticles Particle GenParticle
  add Branch GenParticleFilter/filteredParticles Particle GenParticle


#  add Branch TrackMerger/tracks Track Track
#  add Branch Calorimeter/towers Tower Tower

#  add Branch HCal/eflowTracks EFlowTrack Track
#  add Branch PhotonEnergySmearing/eflowPhotons EFlowPhoton Tower
#  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenMissingET/momentum GenMissingET MissingET

  add Branch UniqueObjectFinder/jets Jet Jet
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/photons Photon Photon
  add Branch UniqueObjectFinder/muons Muon Muon

  add Branch FatJetFinder/jets FatJet Jet

  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}
