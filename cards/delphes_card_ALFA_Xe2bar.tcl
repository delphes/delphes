# ALFA Delphes card, ARC gas variant: Xe2bar; the ARC simulated by the ARCTracing module (configuration B: aerogel n = 1.030, 80 mm SiPM arrays)
set RandomSeed 123

#####################################################################
# FCC-ee ALFA detector model  (ALFA_o1_v01-like, with ARC PID)
#
# Derived from the IDEA card. Tracker, ARC and ECAL numbers follow the
# ALFA full simulation in k4geo (FCCee/ALFA/compact/ALFA_o1_v01) and the
# ALFA talks of Oct 2026 (Reco workshop 26 Oct, ARC pre-demonstrator 16 Oct).
#
#   tracker : all-MAPS. 5-layer vertex (r = 1.3-3.3 cm, 140 mrad acceptance,
#             3 um), 4 double-layer barrel stations (r = 0.77-1.77 m, 3 um,
#             0.2% X0 per plane), 4 double disks per side (|z| = 1.39-2.02 m)
#   ARC     : barrel r = 1.90-2.10 m, endcaps |z| = 2.00-2.20 m, Xe 2 bar +
#             1 cm aerogel; 5.8% X0 passive material in the track fit
#   ECAL    : GRAiNITA, 1.7%/sqrt(E) (+) 1.1%, 0.2 mrad angular at high E
#   HCAL    : SiPM-on-tile in the yoke; resolution is a PLACEHOLDER
#   B field : 3 T
#
# ARC response: the ray-tracing ARCTracing module (modules/ARCTracing), with the
# cell table cards/arc/cells_ALFA_o1_v01.txt and the optical tables
# cards/arc/arc_optics_tables.txt of the k4geo full simulation.
#####################################################################

set B 3.0
## inner radius / half length of the ECAL = where track propagation stops [m]
set R 2.20
set HL 2.30

## ARC envelope [m], as in k4geo ALFA_o1_v01
set ARCBRMIN 1.90
set ARCBRMAX 2.10
set ARCBZMAX 2.20
set ARCEZMIN 2.00
set ARCEZMAX 2.20
set ARCERMIN 0.28
set ARCERMAX 1.90


#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  BeamSpotSmearing

  TruthVertexFinder
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  TrackMergerPre
  TrackSmearing
  ARCTracing
  TrackMerger
  ForwardLooperTracks
  Calorimeter

  TimeSmearingNeutrals
  TimeOfFlightNeutralHadron

  EFlowTrackMerger
  EFlowMerger

  PhotonEfficiency
  PhotonIsolation

  MuonFilter
  TowerMerger

  ElectronFilter
  ElectronEfficiency
  ElectronIsolation

  MuonEfficiency
  MuonIsolation

  MissingET

  NeutrinoFilter
  GenJetFinder
  GenMissingET

  FastJetFinder
  JetEnergyScale

  GenJetFinderDurhamN2
  FastJetFinderDurhamN2

  JetFlavorAssociation

  BTagging
  CTagging
  TauTagging

  TreeWriter
}

#################################
# Beam spot smearing
#################################

module BeamSpotSmearing BeamSpotSmearing {

  set InputArray Delphes/allParticles

  # beam spot size in x, y, z [m] and bunch time spread [s]
  # FCC-ee Z pole example: 5.96E-06, 23.8E-09, 0.397E-03, 36.3E-12
  set SigmaX 0.0
  set SigmaY 0.0
  set SigmaZ 0.0
  set SigmaT 0.0
}


#################################
# Truth Vertex Finder
#################################

module TruthVertexFinder TruthVertexFinder {

  ## below this distance two vertices are assumed to be merged
  set Resolution 1E-06

  set InputArray Delphes/stableParticles
  set VertexOutputArray vertices
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

  # inner radius of the solenoid, in m
  set Radius $::R

  # half-length: z of the solenoid, in m
  set HalfLength $::HL

  # magnetic field, in T
  set Bz $::B
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
    set InputArray ParticlePropagator/chargedHadrons
    set OutputArray chargedHadrons
    set UseMomentumVector true

    set EfficiencyFormula {
        (abs(eta) > 2.6)                                  * (0.000) +
        (pt < 0.1) * (abs(eta) <= 2.6)       * (0.000) +
        (pt >= 0.1) * (abs(eta) <= 2.6)      * (1.000)
    }
}



##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
    set InputArray ParticlePropagator/electrons
    set OutputArray electrons
    set UseMomentumVector true

    set EfficiencyFormula {
        (abs(eta) > 2.6)                                  * (0.000) +
        (pt < 0.1) * (abs(eta) <= 2.6)                    * (0.000) +
        (pt >= 0.1) * (abs(eta) <= 2.6)                   * (1.000)
    }
}


##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
    set InputArray ParticlePropagator/muons
    set OutputArray muons
    set UseMomentumVector true

    set EfficiencyFormula {
        (abs(eta) > 2.6)                                  * (0.000) +
        (pt < 0.1) * (abs(eta) <= 2.6)       * (0.000) +
        (pt >= 0.1) * (abs(eta) <= 2.6)      * (1.000)
    }
}

##############
# Track merger
##############

module Merger TrackMergerPre {
# add InputArray InputArray
  add InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  add InputArray ElectronTrackingEfficiency/electrons
  add InputArray MuonTrackingEfficiency/muons
  set OutputArray tracks
}



########################################
# Smearing for charged tracks
########################################

module TrackCovariance TrackSmearing {

    set InputArray TrackMergerPre/tracks
    set OutputArray tracks

    ## minimum number of hits to accept a track
    set NMinHits 6

    ## magnetic field
    set Bz $::B

    ## scale factors
    set ElectronScaleFactor  {1.25}


    set DetectorGeometry {

      # Layer type 1 = R (barrel) or 2 = z (forward/backward)
      # Layer label
      # Minimum dimension z for barrel or R for forward
      # Maximum dimension z for barrel or R for forward
      # R/z location of layer
      # Thickness (meters)
      # Radiation length (meters)
      # Number of measurements in layers (1D or 2D)
      # Stereo angle (rad) - 0(pi/2) = axial(z) layer - Upper side
      # Stereo angle (rad) - 0(pi/2) = axial(z) layer - Lower side
      # Resolution Upper side (meters) - 0 = no measurement
      # Resolution Lower side (meters) - 0 = no measurement
      # measurement flag = T, scattering only = F

      # barrel  name       zmin   zmax   r        w (m)      X0        n_meas  th_up (rad) th_down (rad)    reso_up (m)   reso_down (m)  flag

      # beam pipe (AlBeMet + paraffin + Au, 0.68% X0 at 90 deg as in k4geo)
      1 PIPE -100 100 0.01 0.00241 0.35276 0 0 0 0 0 0
      # vertex detector: 5 curved MAPS layers, 10 um pixels, 140 mrad acceptance, ~0.05% X0 per layer
      1 VTX -0.092 0.092 0.013 5e-05 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 VTX -0.128 0.128 0.018 5e-05 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 VTX -0.163 0.163 0.023 5e-05 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 VTX -0.199 0.199 0.028 5e-05 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 VTX -0.234 0.234 0.033 5e-05 0.0937 2 0 1.5708 3e-06 3e-06 1
      # outer tracker barrel: 4 stations, each two MAPS planes on the faces of a 5 cm foam stave (k4geo OuterTrackerBarrel_o1_v00);
      # plane thickness per station from the k4geo material scan at 90 deg (T04: 0.46 / 0.39 / 0.40 / 0.52 % X0 per station)
      1 OTB -1.345 1.345 0.767 0.0002155 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 OTB -1.345 1.345 0.82 0.0002155 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 OTB -1.345 1.345 0.893 0.0001827 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 OTB -1.345 1.345 0.946 0.0001827 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 OTB -1.345 1.345 1.587 0.0001874 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 OTB -1.345 1.345 1.64 0.0001874 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 OTB -1.345 1.345 1.714 0.0002436 0.0937 2 0 1.5708 3e-06 3e-06 1
      1 OTB -1.345 1.345 1.767 0.0002436 0.0937 2 0 1.5708 3e-06 3e-06 1
      # ARC barrel as passive material: 20 cm, 6.0% X0 at 90 deg (aerogel n = 1.030) (k4geo T04 material scan)
      1 ARCB -2.2 2.2 2 0.2 3.32 0 0 0 0 0 0
      # outer tracker endcaps: 4 double disks per side (ALFA Delphes card, Jul 2026; engineering design in progress)
      2 OTE 0.15 1.78 -2.045 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 -1.995 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 -1.835 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 -1.785 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 -1.625 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 -1.575 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 -1.415 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 -1.365 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 1.365 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 1.415 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 1.575 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 1.625 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 1.785 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 1.835 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 1.995 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      2 OTE 0.15 1.78 2.045 0.000187 0.0937 2 0 1.5708 3e-06 3e-06 1
      # ARC endcaps as passive material
      2 ARCE 0.28 1.9 -2.1 0.2 3.32 0 0 0 0 0 0
      2 ARCE 0.28 1.9 2.1 0.2 3.32 0 0 0 0 0 0

    }

}

#########################################
# ARC: Array of RICH Cells, per-track PID
#########################################

module ARCTracing ARCTracing {
  ## ARC fast simulation by ray tracing (ARCTracing): photons generated along the helix in the aerogel and the gas from
  ## the k4geo optical tables, traced through the cell optics of the k4geo geometry (mirror, aerogel slab, SiPM array,
  ## pixels) and reconstructed as in the full-simulation analysis. Configuration B: aerogel n = 1.030, 80 mm SiPM arrays; gas Xe2bar.
  set InputArray TrackSmearing/tracks
  set OutputArray tracks
  set Bz $::B

  set CellFile cards/arc/cells_ALFA_o1_v01.txt
  set OpticsFile cards/arc/arc_optics_tables.txt

  ## radiators: index tables (n - 1 scaled by *IndexScale). Aerogel attenuation as in the Geant4 material of the full
  ## simulation (checked with a ddsim DEBUG dump): Rayleigh scattering with the BelleII table, no absorption (DD4hep keeps
  ## the first RAYLEIGH entry of Aerogel_PFRICH and drops ABSLENGTH); scattered photons are propagated, not killed.
  set GasIndexTable RINDEX__Xenon_ARC
  set GasIndexScale 1.0
  set GasAbsorptionTable ""
  set AerogelIndexTable RINDEX__Aerogel_PFRICH
  set AerogelIndexScale 1.000000
  set AerogelAbsorptionTable ""
  set AerogelRayleighTable RAYLEIGH__Aerogel_BelleII
  set AerogelThickness 0.0        ;# mm, 0: from the geometry (10 mm)
  set AerogelSlabInPath true      ;# the photons cross the aerogel again on the way to the sensor

  ## photon detection and readout
  set PDETable ARC_SiPM_QuantumEfficiency
  set PDEScale 1.0
  set MirrorReflectivityTable REFLECTIVITY__ARC_Mirror
  set YieldScale 1.0
  set PixelPitch 0.8              ;# mm
  set SensorSizeScale 1.0         ;# design knob: SiPM array size relative to the geometry
  set SensorShift 0.0             ;# mm along the cell axis
  set MirrorRadiusScale 1.0       ;# design knob: curvature radius, apex fixed
  set DarkCountsPerSensor 0.0

  ## reconstruction (as the full-simulation analysis): emission at mid-radiator, window around the hypothesis, 3-sigma clipping
  set RecoEnergy 3.0              ;# eV, the indices used by the reconstruction
  set RecoRadiusGas 2000.0
  set RecoRadiusAerogel 1925.0
  set RecoZGas 2110.0
  set RecoZAerogel 2025.0
  set Window 0.025
  set MinPe 3.0
  set GasSigmaSingleFloor 1.0E-3
  set AerogelSigmaSingleFloor 1.5E-3
  set GasSigmaCorr 0.0
  set AerogelSigmaCorr 0.0
  set SigmaTrack 0.0
  set SharedCells true            ;# every track is reconstructed against the pooled photons of its cell
  set UseTrueTrack false
  set UseTrueMomentumForPID false
  set TrueHypothesisAngle false
}

##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ARCTracing/tracks
  set OutputArray tracks
}


######################
# Looper Selection
######################

module Efficiency ForwardLooperTracks  {
  set InputArray TrackMerger/tracks
  set OutputArray tracks
  set UseMomentumVector false

  ## select looping tracks that end up in position |eta| > 3.142 (lost by calo)
  set EfficiencyFormula {
    (abs(eta) > 3.0 )                                 * (1.000) +
    (abs(eta) <= 3.0 )                                * (0.000)
  }

}


#############
# Calorimeter
#############
module DualReadoutCalorimeter Calorimeter {
    set ParticleInputArray ParticlePropagator/stableParticles
    set TrackInputArray TrackMerger/tracks

    set TowerOutputArray towers
    set PhotonOutputArray photons

    set EFlowTrackOutputArray eflowTracks
    set EFlowPhotonOutputArray eflowPhotons
    set EFlowNeutralHadronOutputArray eflowNeutralHadrons

    set ECalMinSignificance 2.0
    set HCalMinSignificance 2.5

    set SmearLogNormal false

    set SmearTowerCenter true
    #set SmearTowerCenter false

    ## GRAiNITA barrel photon angular resolution, CLUE clustering (ALFA Reco workshop 26 Oct 2026): a = 0.80, b = 0.55, c = 0.145 mrad; endcap assumed equal
    set PhotonAngularSmearing true

    set PhotonThetaResolutionFormula {
      (abs(eta) <= 0.92)                   * sqrt(0.80^2/energy + 0.55^2/energy^2 + 0.145^2)*1E-3 +
      (abs(eta) > 0.92 && abs(eta) <= 3.0) * sqrt(0.80^2/energy + 0.55^2/energy^2 + 0.145^2)*1E-3
    }

    set PhotonPhiResolutionFormula {
      (abs(eta) <= 0.92)                   * sqrt(0.80^2/energy + 0.55^2/energy^2 + 0.145^2)*1E-3 +
      (abs(eta) > 0.92 && abs(eta) <= 3.0) * sqrt(0.80^2/energy + 0.55^2/energy^2 + 0.145^2)*1E-3
    }

    set PhotonPointing true

    set PointingThetaResolutionFormula { (abs(eta) <= 3.0) * sqrt(111.9^2/energy + 21.0^2)*1E-3 }
    set PointingPhiResolutionFormula   { (abs(eta) <= 3.0) * sqrt(139.4^2/energy + 23.8^2)*1E-3 }

    ## alternative optimistic pointing BDT-corrected
    # set PointingThetaResolutionFormula { (abs(eta) <= 3.0) * sqrt(51.76^2/energy + 11.67^2)*1E-3 }
    # set PointingPhiResolutionFormula   { (abs(eta) <= 3.0) * sqrt(75.76^2/energy + 73.66^2/energy^2 + 16.73^2)*1E-3 }
    set pi [expr {acos(-1)}]

    # Lists of the edges of each tower in eta and phi;
    # each list starts with the lower edge of the first tower;
    # the list ends with the higher edged of the last tower.
    # Barrel:  deta=0.02 towers up to |eta| <= 0.88 ( up to 45°)
    # Endcaps: deta=0.02 towers up to |eta| <= 3.0 (8.6° = 100 mrad)
    # Cell size: about 6 cm x 6 cm

    set EtaPhiRes 0.02
    set EtaMax 3.0

    set pi [expr {acos(-1)}]

    set PhiBins [expr {int($pi/$EtaPhiRes) * 2}]

    set nbins_eta [expr {$EtaMax/$EtaPhiRes} ]
    set nbins_eta [expr {int($nbins_eta)} ]

    for {set i -$nbins_eta} {$i <= $nbins_eta} {incr i} {
      set eta [expr {$i * $EtaPhiRes}]
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
    add EnergyFraction {130} {0.3 0.7}
    add EnergyFraction {3122} {0.3 0.7}


    ## GRAiNITA ECAL: sigma_E/E = 1.7%/sqrt(E) (+) 1.1% at reconstruction level with CLUE (ALFA Reco workshop 26 Oct 2026).
    ## Endcap assumed identical (no endcap design yet).
    set ECalResolutionFormula {
    (abs(eta) <= 0.92 )                     * sqrt(energy^2*0.011^2 + energy*0.017^2 + 0.002^2)+
    (abs(eta) > 0.92 && abs(eta) <= 3.0)    * sqrt(energy^2*0.011^2 + energy*0.017^2 + 0.002^2)
    }


    ## HCAL: SiPM-on-tile scintillator in the instrumented yoke (HGCAL CE-H like).
    ## PLACEHOLDER: 45%/sqrt(E) (+) 3% (CALICE AHCAL-like) until the standalone Geant4 study provides numbers.
    set HCalResolutionFormula {
    (abs(eta) <= 0.92 )                     * sqrt(energy^2*0.03^2 + energy*0.45^2 + 0.05^2)+
    (abs(eta) > 0.92 && abs(eta) <= 3.0)    * sqrt(energy^2*0.03^2 + energy*0.45^2 + 0.05^2)
    }
}

########################################
#   Time Smearing Neutrals
########################################

module TimeSmearing TimeSmearingNeutrals {
  set InputArray Calorimeter/eflowNeutralHadrons
  set OutputArray eflowNeutralHadrons

  # assume constant 30 ps resolution for now
  set TimeResolution {
                       (abs(eta) > 0.0 && abs(eta) <= 3.0)* 30E-12
                     }
}

########################################
#   Time Of Flight Measurement
########################################

module TimeOfFlight TimeOfFlightNeutralHadron {
  set InputArray TimeSmearingNeutrals/eflowNeutralHadrons
  set VertexInputArray TruthVertexFinder/vertices

  set OutputArray eflowNeutralHadrons

  # 0: assume vertex time tV from MC Truth (ideal case)
  # 1: assume vertex time tV = 0
  # 2: calculate vertex time as vertex TOF, assuming tPV=0

  ## TBF (add option to take hard vertex time)
  set VertexTimeMode 1
}


############################
# Energy flow track merger
############################

module Merger EFlowTrackMerger {
# add InputArray InputArray
  add InputArray Calorimeter/eflowTracks
  add InputArray ForwardLooperTracks/tracks
  set OutputArray eflowTracks
}



####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray EFlowTrackMerger/eflowTracks
  add InputArray Calorimeter/eflowPhotons
  add InputArray TimeOfFlightNeutralHadron/eflowNeutralHadrons
  set OutputArray eflow
}


####################
# Tower merger
####################

module Merger TowerMerger {
# add InputArray InputArray
  add InputArray Calorimeter/towers
  add InputArray MuonFilter/muons
  set OutputArray towers
}


###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray Calorimeter/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}
  # efficiency formula for photons
  set EfficiencyFormula {
        (energy < 2.0)                                        * (0.000)+
        (energy >= 2.0) * (abs(eta) <= 0.92)                  * (0.99) +
        (energy >= 2.0) * (abs(eta) >0.88 && abs(eta) <= 3.0) * (0.99) +
        (abs(eta) > 3.0)                                      * (0.000)
  }
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

  set PTRatioMax 9999.
}

#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray EFlowTrackMerger/eflowTracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}

#################
# Muon filter
#################

module PdgCodeFilter MuonFilter {
  set InputArray EFlowTrackMerger/eflowTracks
  set OutputArray muons
  set Invert true
  add PdgCode {13}
  add PdgCode {-13}
}


#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray ElectronFilter/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
  set EfficiencyFormula {
        (energy < 2.0)                                         * (0.000)+
        (energy >= 2.0) * (abs(eta) <= 0.92)                   * (0.99) +
        (energy >= 2.0) * (abs(eta) >0.88 && abs(eta) <= 3.0)  * (0.99) +
        (abs(eta) > 3.0)                                       * (0.000)
  }
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

  set PTRatioMax 9999
}

#################
# Muon efficiency
#################

module Efficiency MuonEfficiency {
  set InputArray MuonFilter/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency as a function of eta and pt}

  # efficiency formula for muons
  set EfficiencyFormula {
        (energy < 2.0)                                         * (0.000)+
        (energy >= 2.0) * (abs(eta) <= 0.92)                   * (0.99) +
        (energy >= 2.0) * (abs(eta) >0.88 && abs(eta) <= 3.0)  * (0.99) +
        (abs(eta) > 3.0)                                       * (0.000)
  }
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

  set PTRatioMax 9999.
}

###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
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

###################################
# Gen Jet finder Durham exclusive
###################################

module FastJetFinder GenJetFinderDurhamN2 {

  set InputArray NeutrinoFilter/filteredParticles
  set OutputArray jets

  # algorithm: 11 ee-durham kT algorithm
  # ref: https://indico.cern.ch/event/1173562/contributions/4929025/attachments/2470068/4237859/2022-06-FCC-jets.pdf
  # to run exclusive njet mode set NJets to int
  # to run exclusive dcut mode set DCut to float
  # if DCut > 0 will run in dcut mode

  set JetAlgorithm 11
  set ExclusiveClustering true
  set NJets 2
  # set DCut 10.0
}

################################
# Jet finder Durham exclusive
################################

module FastJetFinder FastJetFinderDurhamN2 {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 11 ee-durham kT algorithm
  # ref: https://indico.cern.ch/event/1173562/contributions/4929025/attachments/2470068/4237859/2022-06-FCC-jets.pdf
  # to run exclusive njet mode set NJets to int
  # to run exclusive dcut mode set DCut to float
  # if DCut > 0 will run in dcut mode

  set JetAlgorithm 11
  set ExclusiveClustering true
  set NJets 2
  # set DCut 10.0

}




#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray NeutrinoFilter/filteredParticles
  set OutputArray jets

  set JetAlgorithm 10
  set ParameterR 1.5
  set ParameterP -1.0
  set JetPTMin 1.0

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
  set JetAlgorithm 10
  set ParameterR 1.5
  set ParameterP -1.0
  set JetPTMin 1.0


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
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociation {

  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 3.0
}

###########
# b-tagging
###########

module BTagging BTagging {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 0

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.005}

  # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormula {4} {0.01}

  # efficiency formula for b-jets
  add EfficiencyFormula {5} {0.85}
}

###########
# c-tagging
###########

module BTagging CTagging {
  set JetInputArray JetEnergyScale/jets

  set BitNumber 1

  # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.01}

  # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormula {5} {0.05}

  # efficiency formula for b-jets
  add EfficiencyFormula {4} {0.80}
}


#############
# tau-tagging
#############

module TauTagging TauTagging {
  set ParticleInputArray Delphes/allParticles
  set PartonInputArray Delphes/partons
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5
  set TauPTMin 1.0
  set TauEtaMax 3.0

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {0.01}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.85}
}



##################
# ROOT tree writer
##################

# Tracks, towers and eflow objects are not stored by default in the output.
# If needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
    # add Branch InputArray BranchName BranchClass

    add Branch Delphes/allParticles Particle GenParticle
    add Branch TruthVertexFinder/vertices GenVertex Vertex

    add Branch EFlowTrackMerger/eflowTracks EFlowTrack Track
    add Branch TrackMerger/tracks Track Track
    add Branch Calorimeter/eflowPhotons EFlowPhoton Tower
    add Branch TimeOfFlightNeutralHadron/eflowNeutralHadrons EFlowNeutralHadron Tower

    add Branch EFlowMerger/eflow ParticleFlowCandidate ParticleFlowCandidate
    add Branch Calorimeter/towers Tower Tower

    add Branch ElectronEfficiency/electrons Electron Electron
    add Branch MuonEfficiency/muons Muon Muon
    add Branch PhotonEfficiency/photons Photon Photon

    add Branch JetEnergyScale/jets Jet Jet
    add Branch MissingET/momentum MissingET MissingET

    add Branch GenJetFinder/jets GenJet Jet
    add Branch GenMissingET/momentum GenMissingET MissingET

    add Branch GenJetFinderDurhamN2/jets GenJetDurhamN2 Jet
    add Branch FastJetFinderDurhamN2/jets JetDurhamN2 Jet

    # add Info InfoName InfoValue
    add Info Bz $::B
}
