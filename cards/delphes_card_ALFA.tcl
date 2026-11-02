# ALFA Delphes card with the ARC simulated by the ARCTracing module (ray tracing of the Cherenkov photons through the ARC cells).
# The ARC configuration is chosen with the two switches below (set ARCGas, set ARCArrays): everything that depends on it
# (gas index table, SiPM array size, ARC material in the track fit) follows. The cell geometry and the optical tables are
# parameters of the ARCTracing module at the end of this card (k4geo ALFA_o1_v01: the ALFA_o1_v00 detector with the CLD ARC
# cells, 10 mm aerogel n = 1.030, 0.8 mm pixels); no external file is read.
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
#   ARC     : barrel r = 1.90-2.10 m, endcaps |z| = 2.00-2.20 m, gas + 1 cm
#             aerogel; passive material in the track fit from the ARC
#             configuration below (4.4-6.2% X0 at 90 deg)
#   ECAL    : GRAiNITA, 1.7%/sqrt(E) (+) 1.1%, 0.2 mrad angular at high E
#   HCAL    : SiPM-on-tile in the yoke; resolution is a PLACEHOLDER
#   B field : 3 T
#
# ARC response: the ray-tracing ARCTracing module (modules/ARCTracing). Its geometry (vessel, cells, mirrors, sensors)
# and the radiator and detection tables are parameters of the module in this card: change them to describe another ARC;
# a new gas is one more "add OpticsTable {name energy n ...}" line and "set GasIndexTable name".
#####################################################################

set B 3.0

#####################################################################
# ARC configuration
#   ARCGas    : Xe2bar (xenon at 2 bar, n = 1.001313 at 3 eV) | C4F10 (1 bar, n = 1.001331)
#   ARCArrays : baseline  (configuration B: 80 x 80 mm SiPM arrays, the CLD ARC cells)
#               casestudy (configuration C, the ARC case study by Apranik Fatehi: the largest arrays that fit the cells,
#                          89.6 mm barrel / 108.8 mm endcap)
#####################################################################
set ARCGas Xe2bar
set ARCArrays baseline

if {$ARCGas == "Xe2bar"} {
  set ARCGasTable RINDEX__Xenon_ARC
} elseif {$ARCGas == "C4F10"} {
  set ARCGasTable RINDEX__C4F10_PFRICH
} else {
  error "ARCGas must be Xe2bar or C4F10"
}
if {$ARCArrays == "baseline"} {
  set ARCSensorBarrel 80.0
  set ARCSensorEndcap 80.0
} elseif {$ARCArrays == "casestudy"} {
  set ARCSensorBarrel 89.6
  set ARCSensorEndcap 108.8
} else {
  error "ARCArrays must be baseline or casestudy"
}
## radiation length [m] of the 20 cm of ARC as one passive layer of the track fit (k4geo material scans at 90 deg:
## Xe 2 bar 6.0 % X0 with the baseline arrays, 6.2 % with the case-study arrays, C4F10 4.4 %; C4F10 with the case-study
## arrays from the same +0.13 % X0 of the larger arrays)
set ARCX0Table [dict create Xe2bar,baseline 3.32 Xe2bar,casestudy 3.25 C4F10,baseline 4.59 C4F10,casestudy 4.46]
set ARCX0 [dict get $ARCX0Table $ARCGas,$ARCArrays]
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


    set DetectorGeometry [string map [list @ARCX0@ $::ARCX0] {

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
      # ARC barrel as passive material: 20 cm, X0 of the ARC configuration (ARCX0 above)
      1 ARCB -2.2 2.2 2 0.2 @ARCX0@ 0 0 0 0 0 0
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
      2 ARCE 0.28 1.9 -2.1 0.2 @ARCX0@ 0 0 0 0 0 0
      2 ARCE 0.28 1.9 2.1 0.2 @ARCX0@ 0 0 0 0 0 0

    }]

}

#########################################
# ARC: Array of RICH Cells, per-track PID
#########################################

module ARCTracing ARCTracing {
  ## ARC fast simulation by ray tracing (ARCTracing): photons generated along the helix in the aerogel and the gas from
  ## the k4geo optical tables, traced through the cell optics of the k4geo geometry (mirror, aerogel slab, SiPM array,
  ## pixels) and reconstructed as in the full-simulation analysis. Configuration: ARCGas and ARCArrays at the top of the card.
  set InputArray TrackSmearing/tracks
  set OutputArray tracks
  set Bz $::B

  ## ARC geometry, mm. Three ways to give the cells, which add up (each is optional):
  ##  - the construction of the k4geo ARC driver (ARC_geo_o1_v01): barrel of hexagonal pyramid cells, endcaps of hexagonal
  ##    prisms, from the numbers below and the BarrelCellTypes / EndcapCells lines at the end of this module (ALFA_o1_v01: 1203
  ##    cells); leave out BarrelCellTypes or EndcapCells for no barrel or no endcaps;
  ##  - single cells of any layout: add Cell {barrel|endcap numbers...}, one line per cell (origin, axes, mirror, sensor,
  ##    aerogel, cell shape; set DumpCellsFile <file> writes every cell of a configuration in this form, a starting point);
  ##  - a file of such lines without "add Cell" and the braces (set CellFile <file>; a name may follow the region).
  ## Barrel cells must have their axis pointing outward from the beam line and endcap cells along the beam axis: this decides
  ## where a track enters the ARC and where the reconstruction assumes the emission (RecoRadius*, RecoZ*).
  set BarrelInnerRadius 1900.0
  set BarrelOuterRadius 2100.0
  set EndcapZ 2100.0
  set EndcapLength 200.0
  set VesselWallThickness 10.0
  set HexagonSide 148.15
  set BarrelPhiStep 13.333        ;# deg, 27 cells around phi
  set BarrelSectors 27
  set BarrelMaxCell 17            ;# barrel cell numbers 1..17 at z > 0, 2..17 at z < 0
  set EndcapSectors 6
  set EndcapPhiStep 60.0          ;# deg between endcap sectors
  set MirrorThickness 2.0
  set MirrorZOrigin 370.0
  set BarrelMirrorZShrink 6.0
  set MirrorCapAngle 0.897143     ;# rad, 3.14 / 3.5: the mirrors are spherical caps
  set SensorThickness 2.0
  set SensorZOffset 5.0
  set CellAerogelThickness 10.0
  set BarrelSensorSize $::ARCSensorBarrel
  set EndcapSensorSize $::ARCSensorEndcap

  ## radiators: index tables (n - 1 scaled by *IndexScale). Aerogel attenuation as in the Geant4 material of the full
  ## simulation (checked with a ddsim DEBUG dump): Rayleigh scattering with the BelleII table, no absorption (DD4hep keeps
  ## the first RAYLEIGH entry of Aerogel_PFRICH and drops ABSLENGTH); scattered photons are propagated, not killed.
  set GasIndexTable $::ARCGasTable
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
  ## assumed emission point, mm: barrel radius and endcap |z| of the middle of each radiator (change them with the geometry)
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
  set MultiCell true              ;# photons reach the mirrors and arrays of neighbouring cells; each track reconstructed in its main cell
  set RecoMidPathFallback false   ;# study option: gas emission at mid-path when the fixed radius is beyond the main cell's mirror
  ## hypothesis test: expected angle and spread of each hypothesis from test photons traced and reconstructed like the data,
  ## a missing ring weighted by P(N < MinPe); TestPhotons = detected test photons per hypothesis and radiator (at most 4x generated)
  set CalibratedHypotheses true
  set TestPhotons 30
  set AerogelRayleigh true        ;# Rayleigh scattering in the aerogel (table AerogelRayleighTable)
  set AerogelAbsorption true      ;# absorption in the aerogel (table AerogelAbsorptionTable; none here, as the Geant4 material)
  ## validation switches: true track, true momentum for the hypotheses, window around the true angle
  set UseTrueTrack false
  set UseTrueMomentumForPID false
  set TrueHypothesisAngle false
  ## for studies: set DumpCellsFile <file> (the cells as "add Cell" lines), DebugFile <file> (per-track, per-photon and per-hypothesis text
  ## dump), DebugMonoEnergy <eV> (all photons at one energy: no chromatic error), DebugMidEmission true (emission at mid-path)

  ## ---- ARC cell optics (k4geo ARC_geo_o1_v01 construction; mm, rad) ----
  ## barrel cell types: column row mirror-radius mirror-x mirror-z sensor-offset sensor-tilt
  add BarrelCellTypes {0 1 369.008 1.99739 0.484417 2.04862 -0.0204611}
  add BarrelCellTypes {1 1 367.441 -44.1887 2.48451 -20.3282 -0.0931162}
  add BarrelCellTypes {1 2 368.913 -25.092 1.09283 -12.7519 -0.0445222}
  add BarrelCellTypes {2 1 358.554 -81.6915 11.3855 -35.1757 -0.132262}
  add BarrelCellTypes {2 2 363.157 -63.806 6.37531 -28.171 -0.124024}
  add BarrelCellTypes {3 1 346.543 -114.685 23.4365 -47.8341 -0.134323}
  add BarrelCellTypes {3 2 350.086 -95.9577 19.3522 -39.7032 -0.150924}
  add BarrelCellTypes {4 1 330.16 -127.363 42.25 -45.9939 -0.20144}
  add BarrelCellTypes {4 2 336.429 -116.496 33.904 -42.3608 -0.18349}
  add BarrelCellTypes {5 1 325.436 -128.286 50.3257 -32.4965 -0.154886}
  add BarrelCellTypes {5 2 320.575 -123.091 51.37 -36.6375 -0.21376}
  add BarrelCellTypes {6 1 326.104 -139.797 55.998 -29.5479 -0.134243}
  add BarrelCellTypes {6 2 314.895 -129.993 62.4498 -31.5672 -0.229024}
  add BarrelCellTypes {7 1 310.634 -134.094 76.318 -19.8518 -0.259362}
  add BarrelCellTypes {7 2 319.311 -137.129 66.4953 -26.6767 -0.14106}
  add BarrelCellTypes {8 1 350.381 -215.767 78.5187 -82.7883 -0.403075}
  add BarrelCellTypes {8 2 321.726 -141.744 69.8931 -19.6984 -0.0928119}
  add BarrelCellTypes {9 2 254.983 -221.951 135.631 -120.379 -1.38132}
  ## endcap cells of one sector: row column x/side y/apothem mirrored-copy mirror-radius mirror-x mirror-z sensor-offset sensor-tilt
  add EndcapCells {1 2 0 4 0 358.139 -84.6481 11.5558 -40.7355 -0.103914}
  add EndcapCells {1 3 0 6 0 344.967 -107.123 25.0228 -44.7855 -0.137308}
  add EndcapCells {1 4 0 8 0 333.608 -124.894 36.4653 -45.6515 -0.112783}
  add EndcapCells {1 5 0 10 0 323.165 -134.206 53.0051 -44.3281 -0.20109}
  add EndcapCells {1 6 0 12 0 324.725 -138.884 57.7007 -34.8118 -0.159627}
  add EndcapCells {2 2 -1.5 3 0 359.63 -69.3705 10.0598 -30.4485 -0.118494}
  add EndcapCells {2 3 -1.5 5 1 349.695 -99.9593 19.8683 -43.2295 -0.118265}
  add EndcapCells {2 4 -1.5 7 1 339.376 -115.136 30.6836 -42.2966 -0.118536}
  add EndcapCells {2 5 -1.5 9 1 324.686 -124.64 46.9601 -38.3 -0.186657}
  add EndcapCells {2 6 -1.5 11 1 312.448 -125.109 62.1203 -29.028 -0.251983}
  add EndcapCells {2 7 -1.5 13 1 315.405 -130.506 64.6394 -24.065 -0.124305}
  add EndcapCells {3 3 -3 6 0 340.541 -113.788 28.9459 -41.9997 -0.119164}
  add EndcapCells {3 4 -3 8 1 327.512 -126.093 42.7482 -43.3654 -0.143649}
  add EndcapCells {3 5 -3 10 1 313.826 -118.079 57.0635 -25.4509 -0.186348}
  add EndcapCells {3 6 -3 12 1 305.166 -117.67 68.4136 -16.9257 -0.211619}
  add EndcapCells {4 5 -4.5 9 0 317.135 -121.418 53.2098 -28.8431 -0.129972}
  add EndcapCells {4 6 -4.5 11 1 306.595 -116.272 64.8111 -15.9027 -0.19651}
  ## ---- optical tables (Geant4 materials of the geometry): name, then energy [eV] and value pairs ----
  add OpticsTable {ARC_SiPM_QuantumEfficiency 1.1000 0 1.3000 0 1.5500 0.13 1.7710 0.26 2.1010 0.4 2.4650 0.51 2.6550 0.53 2.7430 0.56 2.8370 0.57 2.9380 0.59 3.0020 0.6 3.0840 0.6 3.1630 0.59 3.3780 0.52 3.6470 0.47 3.8500 0.46 4.1470 0.41 4.3200 0.2 4.5000 0.1 4.7000 0 4.9000 0}
  add OpticsTable {RINDEX__C4F10_PFRICH 1.8785 1.0013 1.9239 1.0013 2.1054 1.00131 2.3247 1.00131 2.5950 1.00132 2.9365 1.00133 3.3814 1.00134 3.9852 1.00136 4.8516 1.0014 6.1992 1.00149}
  add OpticsTable {RINDEX__Xenon_ARC 1.8000 1.001257 2.0000 1.001264 2.2000 1.001272 2.4000 1.00128 2.6000 1.00129 2.8000 1.001301 3.0000 1.001313 3.2000 1.001325 3.4000 1.001339 3.6000 1.001354 3.8000 1.001371 4.0000 1.001389 4.2000 1.001408 4.4000 1.00143 4.6000 1.001453 4.8000 1.001478 5.0000 1.001505 5.2000 1.001535 5.4000 1.001568 5.6000 1.001605 5.8000 1.001645 6.0000 1.00169 6.2000 1.00174 6.4000 1.001797 6.6000 1.001862 6.8000 1.001938 7.0000 1.002028}
  add OpticsTable {ABSLENGTH__Xenon_ARC 1.0000 10000 7.5000 10000}
  add OpticsTable {REFLECTIVITY__ARC_Mirror 1.0000 1 7.5000 1}
  add OpticsTable {ABSLENGTH__C4F10_PFRICH 1.7712 6000 1.9239 6000 2.1054 6000 2.3247 6000 2.5950 6000 2.9365 6000 3.3814 6000 3.9852 6000 4.8516 6000 6.1992 6000}
  add OpticsTable {RINDEX__Aerogel_PFRICH 1.8785 1.02876 1.9667 1.02882 2.0549 1.0289 2.1431 1.02898 2.2313 1.02905 2.3194 1.02913 2.4076 1.02921 2.4958 1.0293 2.5840 1.0294 2.6721 1.02949 2.7603 1.02958 2.8485 1.02969 2.9367 1.0298 3.0248 1.02991 3.1130 1.03002 3.2012 1.03014 3.2894 1.03025 3.3775 1.03037 3.4657 1.03051 3.5539 1.03064 3.6421 1.03078 3.7303 1.03092 3.8184 1.03107 3.9066 1.03121 3.9948 1.03137 4.0830 1.03152 4.1711 1.03169 4.2593 1.03186 4.3475 1.03204 4.4357 1.03221 4.5239 1.03239 4.6120 1.03258 4.7002 1.03278 4.7884 1.03297 4.8766 1.03317 4.9647 1.03339 5.0529 1.0336 5.1411 1.03382 5.2293 1.03405 5.3174 1.03429 5.4056 1.03452 5.4938 1.03477 5.5820 1.03502 5.6701 1.03527 5.7583 1.03553 5.8465 1.03581 5.9347 1.03609 6.0229 1.03637 6.1110 1.03666 6.1992 1.03697}
  add OpticsTable {ABSLENGTH__Aerogel_PFRICH 1.8785 140 1.9667 141.973 2.0549 143.776 2.1431 145.431 2.2313 146.955 2.3194 148.364 2.4076 149.669 2.4958 150.882 2.5840 152.012 2.6721 153.067 2.7603 154.055 2.8485 154.982 2.9367 155.854 3.0248 156.674 3.1130 157.448 3.2012 158.18 3.2894 158.872 3.3775 159.528 3.4657 160.15 3.5539 160.742 3.6421 147.916 3.7303 128.139 3.8184 111.378 3.9066 97.121 3.9948 84.948 4.0830 74.518 4.1711 65.552 4.2593 57.819 4.3475 51.13 4.4357 45.327 4.5239 40.278 4.6120 35.873 4.7002 32.019 4.7884 28.641 4.8766 25.67 4.9647 23.054 5.0529 20.742 5.1411 18.698 5.2293 16.884 5.3174 15.272 5.4056 13.837 5.4938 12.557 5.5820 11.413 5.6701 10.389 5.7583 9.47 5.8465 8.645 5.9347 7.902 6.0229 7.233 6.1110 6.629 6.1992 6.082}
  add OpticsTable {RAYLEIGH__Aerogel_BelleII 1.8800 379.7 1.9700 316 2.0500 265.2 2.1400 224.2 2.2300 190.8 2.3200 163.4 2.4100 140.7 2.5000 121.9 2.5800 106.1 2.6700 92.7 2.7600 81.4 2.8500 71.8 2.9400 63.6 3.0200 56.5 3.1100 50.3 3.2000 45 3.2900 40.4 3.3800 36.3 3.4700 32.8 3.5500 29.6 3.6400 26.9 3.7300 24.4 3.8200 22.2 3.9100 20.3 3.9900 18.6 4.0800 17 4.1700 15.6 4.2600 14.4 4.3500 13.2 4.4400 12.2 4.5200 11.3 4.6100 10.5 4.7000 9.7 4.7900 9 4.8800 8.4 4.9600 7.8 5.0500 7.3 5.1400 6.8 5.2300 6.3 5.3200 5.9 5.4100 5.5 5.4900 5.2 5.5800 4.9 5.6700 4.6 5.7600 4.3 5.8500 4 5.9300 3.8 6.0200 3.6 6.1100 3.4 6.2000 3.2}
  add OpticsTable {RAYLEIGH__Aerogel_PFRICH 1.8785 281.107 1.9667 233.984 2.0549 196.334 2.1431 165.962 2.2313 141.242 2.3194 120.958 2.4076 104.188 2.4958 90.226 2.5840 78.527 2.6721 68.663 2.7603 60.301 2.8485 53.174 2.9367 47.07 3.0248 41.816 3.1130 37.277 3.2012 33.336 3.2894 29.903 3.3775 26.9 3.4657 24.265 3.5539 21.946 3.6421 19.896 3.7303 18.08 3.8184 16.468 3.9066 15.03 3.9948 13.746 4.0830 12.596 4.1711 11.564 4.2593 10.637 4.3475 9.799 4.4357 9.043 4.5239 8.358 4.6120 7.738 4.7002 7.172 4.7884 6.659 4.8766 6.191 4.9647 5.762 5.0529 5.37 5.1411 5.011 5.2293 4.681 5.3174 4.379 5.4056 4.1 5.4938 3.844 5.5820 3.606 5.6701 3.386 5.7583 3.184 5.8465 2.996 5.9347 2.822 6.0229 2.66 6.1110 2.51 6.1992 2.37}
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
