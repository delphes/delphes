#######################################
# Order of execution of various modules
#######################################

from delphes_card_CMS import *
import Config.Core as delphes
import math


ExecutionPath = [
    "PileUpMerger",
    "ParticlePropagator",

    "ChargedHadronTrackingEfficiency",
    "ElectronTrackingEfficiency",
    "MuonTrackingEfficiency",

    "ChargedHadronMomentumSmearing",
    "ElectronMomentumSmearing",
    "MuonMomentumSmearing",

    "TrackMerger",

    "ECal",
    "HCal",

    "ElectronFilter",
    "TrackPileUpSubtractor",
    "NeutralTowerMerger",
    "EFlowMergerAllTracks",
    "EFlowMerger",
    "EFlowFilter",

    "NeutrinoFilter",
    "GenJetFinder",
    "GenMissingET",

    "Rho",
    "FastJetFinder",
    "PileUpJetID",
    "JetPileUpSubtractor",

    "JetEnergyScale",

    "PhotonEfficiency",
    "PhotonIsolation",

    "ElectronEfficiency",
    "ElectronIsolation",

    "MuonEfficiency",
    "MuonIsolation",

    "MissingET",

    "JetFlavorAssociation",

    "BTagging",
    "TauTagging",

    "UniqueObjectFinder",

    "ScalarHT",

    "TreeWriter",
]

###############
# PileUp Merger
###############

PileUpMerger = delphes.Module("PileUpMerger",
    InputArray = "Delphes/stableParticles",

    ParticleOutputArray = "stableParticles",
    VertexOutputArray = "vertices",

    PileUpFile = "MinBias.pileup",  # pre-generated minbias input file

    MeanPileUp = 50.,          # average expected pile up
    ZVertexSpread = 0.25,      # maximum spread in the beam direction in m
    TVertexSpread = 800.e-12,  # maximum spread in time in s

    # vertex smearing formula f(z,t) (z,t need to be respectively given in m,s)
    VertexDistributionFormula = "exp(-(t^2/160e-12^2/2))*exp(-(z^2/0.053^2/2))",
)

#################################
# Propagate particles in cylinder
#################################

ParticlePropagator = ParticlePropagator.clone(
    InputArray = "PileUpMerger/stableParticles",
)

#############
#   HCAL
#############

HCal = HCal.clone(
    EnergySignificanceMin = 2.0,
)

####################
# Energy flow merger
####################

EFlowMerger = EFlowMerger.clone(
    InputArray = [
        "TrackPileUpSubtractor/eflowTracks",
        "HCal/eflowTracks",
        "ECal/eflowPhotons",
        "HCal/eflowNeutralHadrons",
    ],
)

######################
# EFlowFilter
######################

EFlowFilter = EFlowFilter.clone(
    InputArray = "EFlowMergerAllTracks/eflow",
)

##########################
# Track pile-up subtractor
##########################

TrackPileUpSubtractor = delphes.Module("TrackPileUpSubtractor",
    InputArray = [
        "HCal/eflowTracks", "eflowTracks",
        "ElectronFilter/electrons", "electrons"
        "MuonMomentumSmearing/muons", "muons",
    ],
    VertexInputArray = "PileUpMerger/vertices",
    # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution
    ZVertexResolution = 0.0001,  # Z vertex resolution in m
)

####################
# Neutral Tower merger
####################

NeutralTowerMerger = delphes.Module("Merger",
    InputArray = [
        "ECal/eflowPhotons",
        "HCal/eflowNeutralHadrons",
    ],
    OutputArray = "towers",
)

####################
# Energy flow merger
####################

EFlowMergerAllTracks = delphes.Module("Merger",
    InputArray = [
        "HCal/eflowTracks",
        "ECal/eflowPhotons",
        "HCal/eflowNeutralHadrons",
    ],
    OutputArray = "eflow",
)

#############
# Rho pile-up
#############

Rho = delphes.Module("FastJetGridMedianEstimator",
    InputArray = "EFlowMerger/eflow",
    RhoOutputArray = "rho",

    # add GridRange rapmin rapmax drap dphi
    # rapmin - the minimum rapidity extent of the grid
    # rapmax - the maximum rapidity extent of the grid
    # drap - the grid spacing in rapidity
    # dphi - the grid spacing in azimuth
    GridRange = [
        -5.0, -2.5, 1.0, 1.0,
        -2.5, 2.5, 1.0, 1.0,
        2.5, 5.0, 1.0, 1.0,
    ],
)

###########################
# Jet Pile-Up ID
###########################

PileUpJetID = delphes.Module("PileUpJetID",
    JetInputArray = "FastJetFinder/jets",
    TrackInputArray = "HCal/eflowTracks",
    NeutralInputArray = "NeutralTowerMerger/towers",

    VertexInputArray = "PileUpMerger/vertices",
    # assume perfect pile-up subtraction for tracks with |z| > fZVertexResolution
    ZVertexResolution = 0.0001,  # Z vertex resolution in m

    OutputArray = "jets",

    UseConstituents = False,
    ParameterR = 0.5,

    JetPTMin = 20.0,
)

###########################
# Jet Pile-Up Subtraction
###########################

JetPileUpSubtractor = delphes.Module("JetPileUpSubtractor",
    JetInputArray = "PileUpJetID/jets",
    RhoInputArray = "Rho/rho",

    OutputArray = "jets",

    JetPTMin = 20.0,
)

##################
# Photon isolation
##################

PhotonIsolation = PhotonIsolation.clone(
    RhoInputArray = "Rho/rho",
)

#####################
# Electron efficiency
#####################

ElectronEfficiency = ElectronEfficiency.clone(
    InputArray = "TrackPileUpSubtractor/electrons",
)

####################
# Electron isolation
####################

ElectronIsolation = ElectronIsolation.clone(
    RhoInputArray = "Rho/rho",
)

#################
# Muon efficiency
#################

MuonEfficiency = MuonEfficiency.clone(
    InputArray = "TrackPileUpSubtractor/muons",
)

################
# Muon isolation
################

MuonIsolation = MuonIsolation.clone(
    RhoInputArray = "Rho/rho",
)

###################
# Missing ET merger
###################

MissingET = MissingET.clone(
    InputArray = ["EFlowMergerAllTracks/eflow"],
)

##################
# Jet Energy Scale
##################

JetEnergyScale = JetEnergyScale.clone(
    InputArray = "JetPileUpSubtractor/jets",
    ScaleFormula = "1.0",  # scale formula for jets
)

##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

TreeWriter = delphes.Module("TreeWriter",
    # InputArray BranchName BranchClass
    Branch = [
        "Delphes/allParticles", "Particle", "GenParticle",

        "TrackMerger/tracks", "Track", "Track",
        "Calorimeter/towers", "Tower", "Tower",

        "HCal/eflowTracks", "EFlowTrack",", ""Track",
        "ECal/eflowPhotons", "EFlowPhoton", "Tower",
        "HCal/eflowNeutralHadrons", "EFlowNeutralHadron", "Tower",

        "GenJetFinder/jets", "GenJet", "Jet",
        "GenMissingET/momentum", "GenMissingET", "MissingET",

        "UniqueObjectFinder/jets", "Jet", "Jet",
        "UniqueObjectFinder/electrons", "Electron", "Electron",
        "UniqueObjectFinder/photons", "Photon", "Photon",
        "UniqueObjectFinder/muons", "Muon", "Muon",

        "MissingET/momentum", "MissingET", "MissingET",
        "ScalarHT/energy", "ScalarHT", "ScalarHT",
        "Rho/rho", "Rho", "Rho",
        "PileUpMerger/vertices", "Vertex", "Vertex",
    ],
)
