set RandomSeed 123

####################################################################                                l
# FCC-ee IDEA detector model
#
# Authors: Elisa Fontanesi, Lorenzo Pezzotti, Massimiliano Antonello, Michele Selvaggi
# email: efontane@bo.infn.it,
#        lorenzo.pezzotti01@universitadipavia.it,
#        m.antonello@uninsubria.it,
#        michele.selvaggi@cern.ch
#####################################################################

## MOD2: set vtx mode timing to MC truth

set B 2.0

## Drift chamber coordinates
set DCHZMIN -2.125
set DCHZMAX 2.125
set DCHRMIN 0.345
set DCHRMAX 2.02


#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  TruthVertexFinder
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  TrackMergerPre
  TrackSmearing
  ClusterCounting
  TimeSmearing
  TimeOfFlight

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

  JetFlavorAssociation

  BTagging
  TauTagging

  TreeWriter
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
  set Radius 2.25

  # half-length: z of the solenoid, in m
  set HalfLength 2.5

  # magnetic field, in T
  set Bz $B
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
    set InputArray ParticlePropagator/chargedHadrons
    set OutputArray chargedHadrons
    set UseMomentumVector true

    set EfficiencyFormula {
        (abs(eta) > 2.56)                                  * (0.000) +
        (pt < 0.1) * (abs(eta) <= 2.56)                    * (0.000) +
        (pt >= 0.1) * (abs(eta) <= 2.56)                   * (1.000)
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
        (abs(eta) > 2.56)                                  * (0.000) +
        (pt < 0.1) * (abs(eta) <= 2.56)                    * (0.000) +
        (pt >= 0.1) * (abs(eta) <= 2.56)                   * (1.000)
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
        (abs(eta) > 2.56)                                  * (0.000) +
        (pt < 0.1) * (abs(eta) <= 2.56)                    * (0.000) +
        (pt >= 0.1) * (abs(eta) <= 2.56)                   * (1.000)
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
    set Bz $B

    ## uses https://raw.githubusercontent.com/selvaggi/FastTrackCovariance/master/GeoIDEA_BASE.txt
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

      # barrel  name       zmin   zmax   r        w (m)      X0        n_meas  th_up (rad) th_down (rad)    reso_up (m)   reso_down (m)  flag

      1        PIPE       -100    100    0.015    0.001655  0.2805     0        0          0                0             0              0
      1        VTXLOW     -0.12   0.12   0.017    0.00028   0.0937     2        0          1.5708           3e-006        3e-006         1
      1        VTXLOW     -0.16   0.16   0.023    0.00028   0.0937     2        0          1.5708           3e-006        3e-006         1
      1        VTXLOW     -0.16   0.16   0.031    0.00028   0.0937     2        0          1.5708           3e-006        3e-006         1
      1        VTXHIGH    -1      1      0.32     0.00047   0.0937     2        0          1.5708           7e-006        7e-006         1
      1        VTXHIGH    -1.05   1.05   0.34     0.00047   0.0937     2        0          1.5708           7e-006        7e-006         1

      # endcap  name       rmin   rmax   z        w (m)      X0        n_meas   th_up (rad)  th_down (rad)   reso_up (m)   reso_down (m) flag

      2        VTXDSK      0.141  0.3   -0.92     0.00028   0.0937     2        0          1.5708           7e-006        7e-006         1
      2        VTXDSK      0.138  0.3   -0.9      0.00028   0.0937     2        0          1.5708           7e-006        7e-006         1
      2        VTXDSK      0.065  0.3   -0.42     0.00028   0.0937     2        0          1.5708           7e-006        7e-006         1
      2        VTXDSK      0.062  0.3   -0.4      0.00028   0.0937     2        0          1.5708           7e-006        7e-006         1
      2        VTXDSK      0.062  0.3    0.4      0.00028   0.0937     2        0          1.5708           7e-006        7e-006         1
      2        VTXDSK      0.065  0.3    0.42     0.00028   0.0937     2        0          1.5708           7e-006        7e-006         1
      2        VTXDSK      0.138  0.3    0.9      0.00028   0.0937     2        0          1.5708           7e-006        7e-006         1
      2        VTXDSK      0.141  0.3    0.92     0.00028   0.0937     2        0          1.5708           7e-006        7e-006         1

      1 DCHCANI $DCHZMIN $DCHZMAX $DCHRMIN 0.0002 0.237223 0 0 0 0 0 0
      1 DCH -2 2 0.36 0.0147748 1400 1 0.0203738 0 0.0001 0 1
      1 DCH -2 2 0.374775 0.0147748 1400 1 -0.0212097 0 0.0001 0 1
      1 DCH -2 2 0.38955 0.0147748 1400 1 0.0220456 0 0.0001 0 1
      1 DCH -2 2 0.404324 0.0147748 1400 1 -0.0228814 0 0.0001 0 1
      1 DCH -2 2 0.419099 0.0147748 1400 1 0.0237172 0 0.0001 0 1
      1 DCH -2 2 0.433874 0.0147748 1400 1 -0.024553 0 0.0001 0 1
      1 DCH -2 2 0.448649 0.0147748 1400 1 0.0253888 0 0.0001 0 1
      1 DCH -2 2 0.463423 0.0147748 1400 1 -0.0262245 0 0.0001 0 1
      1 DCH -2 2 0.478198 0.0147748 1400 1 0.0270602 0 0.0001 0 1
      1 DCH -2 2 0.492973 0.0147748 1400 1 -0.0278958 0 0.0001 0 1
      1 DCH -2 2 0.507748 0.0147748 1400 1 0.0287314 0 0.0001 0 1
      1 DCH -2 2 0.522523 0.0147748 1400 1 -0.029567 0 0.0001 0 1
      1 DCH -2 2 0.537297 0.0147748 1400 1 0.0304025 0 0.0001 0 1
      1 DCH -2 2 0.552072 0.0147748 1400 1 -0.031238 0 0.0001 0 1
      1 DCH -2 2 0.566847 0.0147748 1400 1 0.0320734 0 0.0001 0 1
      1 DCH -2 2 0.581622 0.0147748 1400 1 -0.0329088 0 0.0001 0 1
      1 DCH -2 2 0.596396 0.0147748 1400 1 0.0337442 0 0.0001 0 1
      1 DCH -2 2 0.611171 0.0147748 1400 1 -0.0345795 0 0.0001 0 1
      1 DCH -2 2 0.625946 0.0147748 1400 1 0.0354147 0 0.0001 0 1
      1 DCH -2 2 0.640721 0.0147748 1400 1 -0.0362499 0 0.0001 0 1
      1 DCH -2 2 0.655495 0.0147748 1400 1 0.0370851 0 0.0001 0 1
      1 DCH -2 2 0.67027 0.0147748 1400 1 -0.0379202 0 0.0001 0 1
      1 DCH -2 2 0.685045 0.0147748 1400 1 0.0387552 0 0.0001 0 1
      1 DCH -2 2 0.69982 0.0147748 1400 1 -0.0395902 0 0.0001 0 1
      1 DCH -2 2 0.714595 0.0147748 1400 1 0.0404252 0 0.0001 0 1
      1 DCH -2 2 0.729369 0.0147748 1400 1 -0.04126 0 0.0001 0 1
      1 DCH -2 2 0.744144 0.0147748 1400 1 0.0420949 0 0.0001 0 1
      1 DCH -2 2 0.758919 0.0147748 1400 1 -0.0429296 0 0.0001 0 1
      1 DCH -2 2 0.773694 0.0147748 1400 1 0.0437643 0 0.0001 0 1
      1 DCH -2 2 0.788468 0.0147748 1400 1 -0.044599 0 0.0001 0 1
      1 DCH -2 2 0.803243 0.0147748 1400 1 0.0454336 0 0.0001 0 1
      1 DCH -2 2 0.818018 0.0147748 1400 1 -0.0462681 0 0.0001 0 1
      1 DCH -2 2 0.832793 0.0147748 1400 1 0.0471025 0 0.0001 0 1
      1 DCH -2 2 0.847568 0.0147748 1400 1 -0.0479369 0 0.0001 0 1
      1 DCH -2 2 0.862342 0.0147748 1400 1 0.0487713 0 0.0001 0 1
      1 DCH -2 2 0.877117 0.0147748 1400 1 -0.0496055 0 0.0001 0 1
      1 DCH -2 2 0.891892 0.0147748 1400 1 0.0504397 0 0.0001 0 1
      1 DCH -2 2 0.906667 0.0147748 1400 1 -0.0512738 0 0.0001 0 1
      1 DCH -2 2 0.921441 0.0147748 1400 1 0.0521079 0 0.0001 0 1
      1 DCH -2 2 0.936216 0.0147748 1400 1 -0.0529418 0 0.0001 0 1
      1 DCH -2 2 0.950991 0.0147748 1400 1 0.0537757 0 0.0001 0 1
      1 DCH -2 2 0.965766 0.0147748 1400 1 -0.0546095 0 0.0001 0 1
      1 DCH -2 2 0.980541 0.0147748 1400 1 0.0554433 0 0.0001 0 1
      1 DCH -2 2 0.995315 0.0147748 1400 1 -0.056277 0 0.0001 0 1
      1 DCH -2 2 1.01009 0.0147748 1400 1 0.0571106 0 0.0001 0 1
      1 DCH -2 2 1.02486 0.0147748 1400 1 -0.0579441 0 0.0001 0 1
      1 DCH -2 2 1.03964 0.0147748 1400 1 0.0587775 0 0.0001 0 1
      1 DCH -2 2 1.05441 0.0147748 1400 1 -0.0596108 0 0.0001 0 1
      1 DCH -2 2 1.06919 0.0147748 1400 1 0.0604441 0 0.0001 0 1
      1 DCH -2 2 1.08396 0.0147748 1400 1 -0.0612773 0 0.0001 0 1
      1 DCH -2 2 1.09874 0.0147748 1400 1 0.0621104 0 0.0001 0 1
      1 DCH -2 2 1.11351 0.0147748 1400 1 -0.0629434 0 0.0001 0 1
      1 DCH -2 2 1.12829 0.0147748 1400 1 0.0637763 0 0.0001 0 1
      1 DCH -2 2 1.14306 0.0147748 1400 1 -0.0646092 0 0.0001 0 1
      1 DCH -2 2 1.15784 0.0147748 1400 1 0.0654419 0 0.0001 0 1
      1 DCH -2 2 1.17261 0.0147748 1400 1 -0.0662746 0 0.0001 0 1
      1 DCH -2 2 1.18739 0.0147748 1400 1 0.0671071 0 0.0001 0 1
      1 DCH -2 2 1.20216 0.0147748 1400 1 -0.0679396 0 0.0001 0 1
      1 DCH -2 2 1.21694 0.0147748 1400 1 0.068772 0 0.0001 0 1
      1 DCH -2 2 1.23171 0.0147748 1400 1 -0.0696042 0 0.0001 0 1
      1 DCH -2 2 1.24649 0.0147748 1400 1 0.0704364 0 0.0001 0 1
      1 DCH -2 2 1.26126 0.0147748 1400 1 -0.0712685 0 0.0001 0 1
      1 DCH -2 2 1.27604 0.0147748 1400 1 0.0721005 0 0.0001 0 1
      1 DCH -2 2 1.29081 0.0147748 1400 1 -0.0729324 0 0.0001 0 1
      1 DCH -2 2 1.30559 0.0147748 1400 1 0.0737642 0 0.0001 0 1
      1 DCH -2 2 1.32036 0.0147748 1400 1 -0.0745958 0 0.0001 0 1
      1 DCH -2 2 1.33514 0.0147748 1400 1 0.0754274 0 0.0001 0 1
      1 DCH -2 2 1.34991 0.0147748 1400 1 -0.0762589 0 0.0001 0 1
      1 DCH -2 2 1.36468 0.0147748 1400 1 0.0770903 0 0.0001 0 1
      1 DCH -2 2 1.37946 0.0147748 1400 1 -0.0779215 0 0.0001 0 1
      1 DCH -2 2 1.39423 0.0147748 1400 1 0.0787527 0 0.0001 0 1
      1 DCH -2 2 1.40901 0.0147748 1400 1 -0.0795837 0 0.0001 0 1
      1 DCH -2 2 1.42378 0.0147748 1400 1 0.0804147 0 0.0001 0 1
      1 DCH -2 2 1.43856 0.0147748 1400 1 -0.0812455 0 0.0001 0 1
      1 DCH -2 2 1.45333 0.0147748 1400 1 0.0820762 0 0.0001 0 1
      1 DCH -2 2 1.46811 0.0147748 1400 1 -0.0829068 0 0.0001 0 1
      1 DCH -2 2 1.48288 0.0147748 1400 1 0.0837373 0 0.0001 0 1
      1 DCH -2 2 1.49766 0.0147748 1400 1 -0.0845677 0 0.0001 0 1
      1 DCH -2 2 1.51243 0.0147748 1400 1 0.0853979 0 0.0001 0 1
      1 DCH -2 2 1.52721 0.0147748 1400 1 -0.086228 0 0.0001 0 1
      1 DCH -2 2 1.54198 0.0147748 1400 1 0.087058 0 0.0001 0 1
      1 DCH -2 2 1.55676 0.0147748 1400 1 -0.0878879 0 0.0001 0 1
      1 DCH -2 2 1.57153 0.0147748 1400 1 0.0887177 0 0.0001 0 1
      1 DCH -2 2 1.58631 0.0147748 1400 1 -0.0895474 0 0.0001 0 1
      1 DCH -2 2 1.60108 0.0147748 1400 1 0.0903769 0 0.0001 0 1
      1 DCH -2 2 1.61586 0.0147748 1400 1 -0.0912063 0 0.0001 0 1
      1 DCH -2 2 1.63063 0.0147748 1400 1 0.0920356 0 0.0001 0 1
      1 DCH -2 2 1.64541 0.0147748 1400 1 -0.0928647 0 0.0001 0 1
      1 DCH -2 2 1.66018 0.0147748 1400 1 0.0936937 0 0.0001 0 1
      1 DCH -2 2 1.67495 0.0147748 1400 1 -0.0945226 0 0.0001 0 1
      1 DCH -2 2 1.68973 0.0147748 1400 1 0.0953514 0 0.0001 0 1
      1 DCH -2 2 1.7045 0.0147748 1400 1 -0.09618 0 0.0001 0 1
      1 DCH -2 2 1.71928 0.0147748 1400 1 0.0970085 0 0.0001 0 1
      1 DCH -2 2 1.73405 0.0147748 1400 1 -0.0978369 0 0.0001 0 1
      1 DCH -2 2 1.74883 0.0147748 1400 1 0.0986651 0 0.0001 0 1
      1 DCH -2 2 1.7636 0.0147748 1400 1 -0.0994932 0 0.0001 0 1
      1 DCH -2 2 1.77838 0.0147748 1400 1 0.100321 0 0.0001 0 1
      1 DCH -2 2 1.79315 0.0147748 1400 1 -0.101149 0 0.0001 0 1
      1 DCH -2 2 1.80793 0.0147748 1400 1 0.101977 0 0.0001 0 1
      1 DCH -2 2 1.8227 0.0147748 1400 1 -0.102804 0 0.0001 0 1
      1 DCH -2 2 1.83748 0.0147748 1400 1 0.103632 0 0.0001 0 1
      1 DCH -2 2 1.85225 0.0147748 1400 1 -0.104459 0 0.0001 0 1
      1 DCH -2 2 1.86703 0.0147748 1400 1 0.105286 0 0.0001 0 1
      1 DCH -2 2 1.8818 0.0147748 1400 1 -0.106113 0 0.0001 0 1
      1 DCH -2 2 1.89658 0.0147748 1400 1 0.10694 0 0.0001 0 1
      1 DCH -2 2 1.91135 0.0147748 1400 1 -0.107766 0 0.0001 0 1
      1 DCH -2 2 1.92613 0.0147748 1400 1 0.108593 0 0.0001 0 1
      1 DCH -2 2 1.9409 0.0147748 1400 1 -0.109419 0 0.0001 0 1
      1 DCH -2 2 1.95568 0.0147748 1400 1 0.110246 0 0.0001 0 1
      1 DCH -2 2 1.97045 0.0147748 1400 1 -0.111072 0 0.0001 0 1
      1 DCH -2 2 1.98523 0.0147748 1400 1 0.111898 0 0.0001 0 1
      1 DCH -2 2 2 0.0147748 1400 1 -0.112723 0 0.0001 0 1
      1 DCHCANO $DCHZMIN $DCHZMAX $DCHRMAX $DCHRMAX 0.02 1.667 0 0 0 0 0 0
      1 BSILWRP -2.35 2.35 2.04 0.00047 0.0937 2 0 1.5708 7e-006 9e-005 1
      1 BSILWRP -2.35 2.35 2.06 0.00047 0.0937 2 0 1.5708 7e-006 9e-005 1
      1 MAG -2.5 2.5 2.25 0.05 0.0658 0 0 0 0 0 0
      1 BPRESH -2.55 2.55 2.45 0.02 1 2 0 1.5708 7e-005 0.01 1
      2 DCHWALL $DCHRMIN $DCHRMAX $DCHZMAX 0.25 5.55 0 0 0 0 0 0
      2 DCHWALL $DCHRMIN $DCHRMAX $DCHZMIN 0.25 5.55 0 0 0 0 0 0
      2 FSILWRP 0.354 2.02 -2.32 0.00047 0.0937 2 0 1.5708 7e-006 9e-005 1
      2 FSILWRP 0.35 2.02 -2.3 0.00047 0.0937 2 0 1.5708 7e-006 9e-005 1
      2 FSILWRP 0.35 2.02 2.3 0.00047 0.0937 2 0 1.5708 7e-006 9e-005 1
      2 FSILWRP 0.354 2.02 2.32 0.00047 0.0937 2 0 1.5708 7e-006 9e-005 1
      2 FRAD 0.38 2.09 2.49 0.0043 0.005612 0 0 0 0 0 0
      2 FRAD 0.38 2.09 -2.49 0.0043 0.005612 0 0 0 0 0 0
      2 FPRESH 0.39 2.43 -2.55 0.02 1 2 0 1.5708 7e-005 0.01 1
      2 FPRESH 0.39 2.43 2.55 0.02 1 2 0 1.5708 7e-005 0.01 1
    }

}

###################
# Cluster Counting
###################

module ClusterCounting ClusterCounting {

  add InputArray TrackSmearing/tracks
  set OutputArray tracks

  set Bz $B

  ## check that these are consistent with DCHCANI/DCHNANO parameters in TrackCovariance module
  set Rmin $DCHRMIN
  set Rmax $DCHRMAX
  set Zmin $DCHZMIN
  set Zmax $DCHZMAX

  # gas mix option:
  # 0:  Helium 90% - Isobutane 10%
  # 1:  Helium 100%
  # 2:  Argon 50% - Ethane 50%
  # 3:  Argon 100%

  set GasOption 0

}


########################################
#   Time Smearing MIP
########################################

module TimeSmearing TimeSmearing {
  set InputArray ClusterCounting/tracks
  set OutputArray tracks

  # assume constant 30 ps resolution for now
  set TimeResolution {
                       (abs(eta) > 0.0 && abs(eta) <= 3.0)* 30E-12
                     }
}


########################################
#   Time Of Flight Measurement
########################################

module TimeOfFlight TimeOfFlight {
  set InputArray TimeSmearing/tracks
  set VertexInputArray TruthVertexFinder/vertices

  set OutputArray tracks

  # 0: assume vertex time tV from MC Truth (ideal case)
  # 1: assume vertex time tV = 0
  # 2: calculate vertex time as vertex TOF, assuming tPV=0

  set VertexTimeMode 0
}

##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray TimeOfFlight/tracks
  set OutputArray tracks
}


######################
# Looper Selection
######################

module Efficiency ForwardLooperTracks  {
  set InputArray TrackMerger/tracks
  set OutputArray tracks
  set UseMomentumVector False

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

  set ECalEnergyMin 0.5
  set HCalEnergyMin 0.5
  set ECalEnergySignificanceMin 3.0
  set HCalEnergySignificanceMin 3.0

  set EnergyMin 0.5
  set EnergySignificanceMin 3.0

  #set SmearTowerCenter true
  set SmearTowerCenter false
    set pi [expr {acos(-1)}]

    # Lists of the edges of each tower in eta and phi;
    # each list starts with the lower edge of the first tower;
    # the list ends with the higher edged of the last tower.
    # Barrel:  deta=0.02 towers up to |eta| <= 0.88 ( up to 45°)
    # Endcaps: deta=0.02 towers up to |eta| <= 3.0 (8.6° = 100 mrad)
    # Cell size: about 6 cm x 6 cm

    #barrel:
    set PhiBins {}
    for {set i -120} {$i <= 120} {incr i} {
        add PhiBins [expr {$i * $pi/120}]
    }
    #deta=0.02 units for |eta| <= 0.88
    for {set i -44} {$i < 45} {incr i} {
        set eta [expr {$i * 0.02}]
        add EtaPhiBins $eta $PhiBins
    }

    #endcaps:
    set PhiBins {}
    for {set i -120} {$i <= 120} {incr i} {
        add PhiBins [expr {$i* $pi/120}]
    }
    #deta=0.02 units for 0.88 < |eta| <= 3.0
    #first, from -3.0 to -0.88
    for {set i 0} {$i <=106} {incr i} {
        set eta [expr {-3.00 + $i * 0.02}]
        add EtaPhiBins $eta $PhiBins
    }
    #same for 0.88 to 3.0
    for  {set i 1} {$i <=106} {incr i} {
        set eta [expr {0.88 + $i * 0.02}]
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


    # set ECalResolutionFormula {resolution formula as a function of eta and energy}
    set ECalResolutionFormula {
    (abs(eta) <= 0.88 )                     * sqrt(energy^2*0.01^2 + energy*0.11^2)+
    (abs(eta) > 0.88 && abs(eta) <= 3.0)    * sqrt(energy^2*0.01^2 + energy*0.11^2)
    }

    # set HCalResolutionFormula {resolution formula as a function of eta and energy}
    set HCalResolutionFormula {
    (abs(eta) <= 0.88 )                     * sqrt(energy^2*0.01^2 + energy*0.30^2)+
    (abs(eta) > 0.88 && abs(eta) <= 3.0)    * sqrt(energy^2*0.01^2 + energy*0.30^2)
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
        (energy >= 2.0) * (abs(eta) <= 0.88)                  * (0.99) +
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
        (energy >= 2.0) * (abs(eta) <= 0.88)                   * (0.99) +
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
        (energy >= 2.0) * (abs(eta) <= 0.88)                   * (0.99) +
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
  add EfficiencyFormula {0} {0.01}

  # efficiency formula for c-jets (misidentification rate)
  add EfficiencyFormula {4} {0.10}

  # efficiency formula for b-jets
  add EfficiencyFormula {5} {0.80}
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
  add EfficiencyFormula {0} {0.001}
  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {0.6}
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
    add Branch TrackSmearing/tracks Track Track
    add Branch Calorimeter/eflowPhotons EFlowPhoton Tower
    add Branch TimeOfFlightNeutralHadron/eflowNeutralHadrons EFlowNeutralHadron Tower

    add Branch EFlowMerger/eflow ParticleFlowCandidate ParticleFlowCandidate

    add Branch ElectronEfficiency/electrons Electron Electron
    add Branch MuonEfficiency/muons Muon Muon
    add Branch PhotonEfficiency/photons Photon Photon

    add Branch JetEnergyScale/jets Jet Jet
    add Branch MissingET/momentum MissingET MissingET

    add Branch GenJetFinder/jets GenJet Jet
    add Branch GenMissingET/momentum GenMissingET MissingET

    # add Info InfoName InfoValue
    add Info Bz $B
}
