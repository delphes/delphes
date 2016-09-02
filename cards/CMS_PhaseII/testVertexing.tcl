set MaxEvents 100
#set RandomSeed 123


#
#  Phase II - Pile-Up
#
#  Main authors: Michele Selvaggi (UCL)
#
#  Released on:
#
#  Version: v02 beta - test TrackSmearing and vertexing
#
#
#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {

  BeamSpotFilter

  PileUpMerger
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  TrackMerger
  TrackSmearing

  TimeSmearing

  VertexFinder
  VertexFinderDA4D

  TreeWriter
}


#######################
# GenBeamSpotFilter
# Saves a particle intended to represent the beamspot
#######################

module BeamSpotFilter BeamSpotFilter {
    set InputArray Delphes/stableParticles
    set OutputArray beamSpotParticle

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
  set MeanPileUp 50

  # 0-poisson, 1-uniform, 2-delta
  set PileUpDistribution 2

  # maximum spread in the beam direction in m
  set ZVertexSpread 0.25

  # maximum spread in time in s
  set TVertexSpread 800E-12

  # vertex smearing formula f(z,t) (z,t need to be respectively given in m,s)

  #set VertexDistributionFormula {exp(-(t^2/(2*(0.063/2.99792458E8*exp(-(z^2/(2*(0.063)^2))))^2)))}
  set VertexDistributionFormula {exp(-(t^2/160e-12^2/2))*exp(-(z^2/0.053^2/2))}

  # taking 5.3 cm x 160 ps

  #set VertexDistributionFormula { (abs(t) <= 160e-12) * (abs(z) <= 0.053) * (1.00) +
  #                                (abs(t) >  160e-12) * (abs(z) <= 0.053) * (0.00) +
  #                   	          (abs(t) <= 160e-12) * (abs(z) > 0.053)  * (0.00) +
  # 				                  (abs(t) >  160e-12) * (abs(z) > 0.053)  * (0.00)}

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
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons

  # add EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for charged hadrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.70) +
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.95) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.60) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.85) +
                         (abs(eta) > 2.5)                                                  * (0.00)}
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for electrons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.73) +
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e2) * (0.95) +
                                           (abs(eta) <= 1.5) * (pt > 1.0e2)                * (0.99) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.50) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e2) * (0.83) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e2)                * (0.90) +
                         (abs(eta) > 2.5)                                                  * (0.00)}
}

##########################
# Muon tracking efficiency
##########################

module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # tracking efficiency formula for muons
  set EfficiencyFormula {                                                    (pt <= 0.1)   * (0.00) +
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.75) +
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e3) * (0.99) +
                                           (abs(eta) <= 1.5) * (pt > 1.0e3 )               * (0.99 * exp(0.5 - pt*5.0e-4)) +

                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.70) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e3) * (0.98) +
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e3)                * (0.98 * exp(0.5 - pt*5.0e-4)) +
                         (abs(eta) > 2.5)                                                  * (0.00)}
}



##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  add InputArray ElectronTrackingEfficiency/electrons
  add InputArray MuonTrackingEfficiency/muons
  set OutputArray tracks
}

########################################
#   Smear tracks
########################################

module TrackSmearing TrackSmearing {
  set InputArray TrackMerger/tracks
  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
  set OutputArray tracks
  set ApplyToPileUp true

  # from http://mersi.web.cern.ch/mersi/layouts/.private/Baseline_tilted_200_Pixel_1_1_1/index.html
  source trackResolution.tcl
}

########################################
#   Time Smearing
########################################

module TimeSmearing TimeSmearing {
  set InputArray TrackSmearing/tracks
  set OutputArray tracks

  # assume 20 ps resolution for now
  set TimeResolution 20E-12
}

##################################
# Primary vertex reconstruction
##################################


module VertexFinderDA4D VertexFinderDA4D {
  set InputArray TimeSmearing/tracks
#  set InputArray TrackSmearing/tracks
#  set InputArray TrackMerger/tracks

  set OutputArray tracks
  set VertexOutputArray vertices

  set Verbose 0
  set MinPT 1.0

  # in mm
  set VertexSpaceSize 0.5

  # in s
  set VertexTimeSize 10E-12

  set UseTc 1
  set BetaMax 0.1
  set BetaStop 1.0
  set CoolingFactor 0.8
  set MaxIterations 100

  # in mm
  set DzCutOff 40
  set D0CutOff 30

}

##################################
# Primary vertex reconstruction
##################################


module VertexFinder VertexFinder {

  set InputArray TrackSmearing/tracks
#  set InputArray TimeSmearing/tracks
#  set InputArray TrackMerger/tracks

  set OutputArray tracks
  set VertexOutputArray vertices

  set MinPT 1.0
  set MinNDF 4
  set SeedMinPT 1.0

  set Sigma     3.0
  set MaxEta    10.0
  set GrowSeeds 1

}

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass

  add Branch PileUpMerger/stableParticles Particle GenParticle
  add Branch TimeSmearing/tracks Track Track

  add Branch VertexFinderDA4D/vertices Vertex4D Vertex
  add Branch VertexFinder/vertices Vertex Vertex

  add Branch PileUpMerger/vertices GenVertex Vertex


}

