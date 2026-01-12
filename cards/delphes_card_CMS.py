#######################################
# Order of execution of various modules
#######################################

import Config.Core as delphes
import math


ExecutionPath = [
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

    "Calorimeter",
    "EFlowMerger",
    "EFlowFilter",

    "PhotonEfficiency",
    "PhotonIsolation",

    "ElectronFilter",
    "ElectronEfficiency",
    "ElectronIsolation",

    "ChargedHadronFilter",

    "MuonEfficiency",
    "MuonIsolation",

    "MissingET",

    "NeutrinoFilter",
    "GenJetFinder",
    "GenMissingET",

    "FastJetFinder",
    "FatJetFinder",

    "JetEnergyScale",

    "JetFlavorAssociation",

    "BTagging",
    "TauTagging",

    "UniqueObjectFinder",

    "ScalarHT",

    "TreeWriter",
]

#################################
# Propagate particles in cylinder
#################################

ParticlePropagator = delphes.Module("ParticlePropagator",
    InputArray = "Delphes/stableParticles",

    OutputArray = "stableParticles",
    ChargedHadronOutputArray = "chargedHadrons",
    ElectronOutputArray = "electrons",
    MuonOutputArray = "muons",

    Radius = 1.29,      # radius of the magnetic field coverage, in m
    HalfLength = 3.00,  # half-length of the magnetic field coverage, in m

    Bz = 3.8,  # magnetic field
)

####################################
# Charged hadron tracking efficiency
####################################

ChargedHadronTrackingEfficiency = delphes.Module("Efficiency",
    InputArray = "ParticlePropagator/chargedHadrons",
    OutputArray = "chargedHadrons",

    # add EfficiencyFormula {efficiency formula as a function of eta and pt}

    # tracking efficiency formula for charged hadrons
    EfficiencyFormula = "                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.70) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0)                  * (0.95) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.60) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0)                  * (0.85) + \
                         (abs(eta) > 2.5)                                                  * (0.00)",
)

##############################
# Electron tracking efficiency
##############################

ElectronTrackingEfficiency = delphes.Module("Efficiency",
    InputArray = "ParticlePropagator/electrons",
    OutputArray = "electrons",

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    # tracking efficiency formula for electrons
    EfficiencyFormula = "                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.73) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e2) * (0.95) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e2)                * (0.99) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.50) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e2) * (0.83) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e2)                * (0.90) + \
                         (abs(eta) > 2.5)                                                  * (0.00)",
)

##########################
# Muon tracking efficiency
##########################

MuonTrackingEfficiency = delphes.Module("Efficiency",
    InputArray = "ParticlePropagator/muons",
    OutputArray = "muons",

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    # tracking efficiency formula for muons
    EfficiencyFormula = "                                                    (pt <= 0.1)   * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 0.1   && pt <= 1.0)   * (0.75) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0   && pt <= 1.0e3) * (0.99) + \
                                           (abs(eta) <= 1.5) * (pt > 1.0e3 )               * (0.99 * exp(0.5 - pt*5.0e-4)) + \
                         \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1   && pt <= 1.0)   * (0.70) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0   && pt <= 1.0e3) * (0.98) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 1.0e3)                * (0.98 * exp(0.5 - pt*5.0e-4)) + \
                         (abs(eta) > 2.5)                                                  * (0.00)",
)

########################################
# Momentum resolution for charged tracks
########################################

ChargedHadronMomentumSmearing = delphes.Module("MomentumSmearing",
    InputArray = "ChargedHadronTrackingEfficiency/chargedHadrons",
    OutputArray = "chargedHadrons",

    # set ResolutionFormula {resolution formula as a function of eta and pt}

    # resolution formula for charged hadrons
    # based on arXiv:1405.6569
    ResolutionFormula = "                  (abs(eta) <= 0.5) * (pt > 0.1) * sqrt(0.06^2 + pt^2*1.3e-3^2) + \
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 0.1) * sqrt(0.10^2 + pt^2*1.7e-3^2) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1) * sqrt(0.25^2 + pt^2*3.1e-3^2)",
)

###################################
# Momentum resolution for electrons
###################################

ElectronMomentumSmearing = delphes.Module("MomentumSmearing",
  InputArray = "ElectronTrackingEfficiency/electrons",
  OutputArray = "electrons",

  # set ResolutionFormula {resolution formula as a function of eta and energy}

  # resolution formula for electrons
  # based on arXiv:1502.02701
  ResolutionFormula = "                  (abs(eta) <= 0.5) * (pt > 0.1) * sqrt(0.03^2 + pt^2*1.3e-3^2) + \
                       (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 0.1) * sqrt(0.05^2 + pt^2*1.7e-3^2) + \
                       (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1) * sqrt(0.15^2 + pt^2*3.1e-3^2)",
)

###############################
# Momentum resolution for muons
###############################

MuonMomentumSmearing = delphes.Module("MomentumSmearing",
    InputArray = "MuonTrackingEfficiency/muons",
    OutputArray = "muons",

    # set ResolutionFormula {resolution formula as a function of eta and pt}

    # resolution formula for muons
    ResolutionFormula = "                  (abs(eta) <= 0.5) * (pt > 0.1) * sqrt(0.01^2  + pt^2*1.0e-4^2) + \
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 0.1) * sqrt(0.015^2 + pt^2*1.5e-4^2) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 0.1) * sqrt(0.025^2 + pt^2*3.5e-4^2)",
)

##############
# Track merger
##############

TrackMerger = delphes.Module("Merger",
    InputArray = [
        "ChargedHadronMomentumSmearing/chargedHadrons",
        "ElectronMomentumSmearing/electrons",
        "MuonMomentumSmearing/muons",
    ],
    OutputArray = "tracks",
)

#############
#   ECAL
#############

ECal = delphes.Module("SimpleCalorimeter",
    ParticleInputArray = "ParticlePropagator/stableParticles",
    TrackInputArray = "TrackMerger/tracks",

    TowerOutputArray = "ecalTowers",
    EFlowTrackOutputArray = "eflowTracks",
    EFlowTowerOutputArray = "eflowPhotons",

    IsEcal = True,

    EnergyMin = 0.5,
    EnergySignificanceMin = 2.0,

    SmearTowerCenter = True,

    # lists of the edges of each tower in eta and phi
    # each list starts with the lower edge of the first tower
    # the list ends with the higher edged of the last tower
    EtaPhiBins = [],

    EnergyFraction = [
        0, [0.0],
        11, [1.0],  # energy fractions for e, gamma and pi0
        22, [1.0],
        111, [1.0],
        12, [0.0],  # energy fractions for muon, neutrinos and neutralinos
        13, [0.0],
        14, [0.0],
        16, [0.0],
        1000022, [0.0],
        1000023, [0.0],
        1000025, [0.0],
        1000035, [0.0],
        1000045, [0.0],
        310, [0.3],  # energy fractions for K0short and Lambda
        3122, [0.3],
    ],

    # set ResolutionFormula {resolution formula as a function of eta and energy}

    # for the ECAL barrel (|eta| < 1.5), see hep-ex/1306.2016 and 1502.02701

    # set ECalResolutionFormula {resolution formula as a function of eta and energy}
    # Eta shape from arXiv:1306.2016, Energy shape from arXiv:1502.02701
    ResolutionFormula = "                      (abs(eta) <= 1.5) * (1+0.64*eta^2) * sqrt(energy^2*0.008^2 + energy*0.11^2 + 0.40^2) + \
                             (abs(eta) > 1.5 && abs(eta) <= 2.5) * (2.16 + 5.6*(abs(eta)-2)^2) * sqrt(energy^2*0.008^2 + energy*0.11^2 + 0.40^2) + \
                             (abs(eta) > 2.5 && abs(eta) <= 5.0) * sqrt(energy^2*0.107^2 + energy*2.08^2)",

)
# assume 0.02 x 0.02 resolution in eta,phi in the barrel |eta| < 1.5
# 0.02 unit in eta up to eta = 1.5 (barrel)
for eta in [float(i * 0.0174) for i in range(-85, 87)]:
    ECal["EtaPhiBins"].append([eta,])
    ECal["EtaPhiBins"].append([i * math.pi / 180. for i in range(-180, 181)])

# assume 0.02 x 0.02 resolution in eta,phi in the endcaps 1.5 < |eta| < 3.0 (HGCAL- ECAL)
# 0.02 unit in eta up to eta = 3
for eta in [float(-2.958 + i * 0.0174) for i in range(1, 85)]:
    ECal["EtaPhiBins"].append([eta,])
    ECal["EtaPhiBins"].append([i * math.pi / 180. for i in range(-180, 181)])
for eta in [float(1.4964 + i * 0.0174) for i in range(1, 85)]:
    ECal["EtaPhiBins"].append([eta,])
    ECal["EtaPhiBins"].append([i * math.pi / 180. for i in range(-180, 181)])

# take present CMS granularity for HF
# 0.175 x (0.175 - 0.35) resolution in eta,phi in the HF 3.0 < |eta| < 5.0
for eta in [-5., -4.7, -4.525, -4.35, -4.175, -4., -3.825, -3.65, -3.475, -3.3, -3.125, -2.958,
            3.125, 3.3, 3.475, 3.65, 3.825, 4., 4.175, 4.35, 4.525, 4.7, 5.]:
    ECal["EtaPhiBins"].append([eta,])
    ECal["EtaPhiBins"].append([i * math.pi / 18. for i in range(-18, 19)])

#############
#   HCAL
#############

HCal = delphes.Module("SimpleCalorimeter",
    ParticleInputArray = "ParticlePropagator/stableParticles",
    TrackInputArray = "ECal/eflowTracks",

    TowerOutputArray = "hcalTowers",
    EFlowTrackOutputArray = "eflowTracks",
    EFlowTowerOutputArray = "eflowNeutralHadrons",

    IsEcal = False,

    EnergyMin = 1.0,
    EnergySignificanceMin = 1.0,

    SmearTowerCenter = True,

    pi = math.pi,

    # lists of the edges of each tower in eta and phi
    # each list starts with the lower edge of the first tower
    # the list ends with the higher edged of the last tower

    EtaPhiBins = [],

    EnergyFraction = [
        0, [1.0],  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
        11, [0.0],  # energy fractions for e, gamma and pi0
        22, [0.0],
        111, [0.0],
        12, [0.0],  # energy fractions for muon, neutrinos and neutralinos
        13, [0.0],
        14, [0.0],
        16, [0.0],
        1000022, [0.0],
        1000023, [0.0],
        1000025, [0.0],
        1000035, [0.0],
        1000045, [0.0],
        310, [0.7],  # energy fractions for K0short and Lambda
        3122, [0.7],
    ],

    # set HCalResolutionFormula {resolution formula as a function of eta and energy}
    ResolutionFormula = "                  (abs(eta) <= 3.0) * sqrt(energy^2*0.050^2 + energy*1.50^2) + \
                         (abs(eta) > 3.0 && abs(eta) <= 5.0) * sqrt(energy^2*0.130^2 + energy*2.70^2)",
)

# 5 degrees towers
for eta in [-1.566, -1.479, -1.392, -1.305, -1.218, -1.131, -1.044, -0.957, -0.87, -0.783, -0.696, -0.609, -0.522, -0.435, -0.348, -0.261, -0.174, -0.087,
            0.,
            0.087, 0.174, 0.261, 0.348, 0.435, 0.522, 0.609, 0.696, 0.783, 0.87, 0.957, 1.044, 1.131, 1.218, 1.305, 1.392, 1.479, 1.566, 1.653]:
    HCal["EtaPhiBins"].append([eta,])
    HCal["EtaPhiBins"].append([i * math.pi / 36. for i in range(-36, 37)])
# 10 degrees towers
for eta in [-4.35, -4.175, -4., -3.825, -3.65, -3.475, -3.3, -3.125, -2.95, -2.868, -2.65, -2.5, -2.322, -2.172, -2.043, -1.93, -1.83, -1.74, -1.653,
            1.74, 1.83, 1.93, 2.043, 2.172, 2.322, 2.5, 2.65, 2.868, 2.95, 3.125, 3.3, 3.475, 3.65, 3.825, 4., 4.175, 4.35, 4.525]:
    HCal["EtaPhiBins"].append([eta,])
    HCal["EtaPhiBins"].append([i * math.pi / 18. for i in range(-18, 19)])
# 20 degrees towers
for eta in [-5., -4.7, -4.525, 4.7, 5.]:
    HCal["EtaPhiBins"].append([eta,])
    HCal["EtaPhiBins"].append([i * math.pi / 9. for i in range(-9, 10)])

#################
# Electron filter
#################

ElectronFilter = delphes.Module("PdgCodeFilter",
  InputArray = "HCal/eflowTracks",
  OutputArray = "electrons",
  Invert = True,
  PdgCode = [11, -11],
)

######################
# ChargedHadronFilter
######################

ChargedHadronFilter = delphes.Module("PdgCodeFilter",
  InputArray = "HCal/eflowTracks",
  OutputArray = "chargedHadrons",

  PdgCode = [11, -11, 13, -13],
)

###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

Calorimeter = delphes.Module("Merger",
    InputArray = [
        "ECal/ecalTowers",
        "HCal/hcalTowers",
    ],
    OutputArray = "towers",
)

####################
# Energy flow merger
####################

EFlowMerger = delphes.Module("Merger",
    InputArray = [
        "HCal/eflowTracks",
        "ECal/eflowPhotons",
        "HCal/eflowNeutralHadrons",
    ],
    OutputArray = "eflow",
)

######################
# EFlowFilter
######################

EFlowFilter = delphes.Module("PdgCodeFilter",
    InputArray = "EFlowMerger/eflow",
    OutputArray = "eflow",

    PdgCode = [11, -11, 13, -13],
)

###################
# Photon efficiency
###################

PhotonEfficiency = delphes.Module("Efficiency",
    InputArray = "ECal/eflowPhotons",
    OutputArray = "photons",

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    # efficiency formula for photons
    EfficiencyFormula = "                                      (pt <= 10.0) * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 10.0)  * (0.85) + \
                         (abs(eta) > 2.5)                                   * (0.00)",
)

##################
# Photon isolation
##################

PhotonIsolation = delphes.Module("Isolation",
    CandidateInputArray = "PhotonEfficiency/photons",
    IsolationInputArray = "EFlowFilter/eflow",

    OutputArray = "photons",

    DeltaRMax = 0.5,

    PTMin = 0.5,

    PTRatioMax = 0.12,
)

#####################
# Electron efficiency
#####################

ElectronEfficiency = delphes.Module("Efficiency",
    InputArray = "ElectronFilter/electrons",
    OutputArray = "electrons",

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    # efficiency formula for electrons
    EfficiencyFormula = "                                      (pt <= 10.0) * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 10.0)  * (0.95) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.5) * (pt > 10.0)  * (0.85) + \
                         (abs(eta) > 2.5)                                   * (0.00)",
)

####################
# Electron isolation
####################

ElectronIsolation = delphes.Module("Isolation",
    CandidateInputArray = "ElectronEfficiency/electrons",
    IsolationInputArray = "EFlowFilter/eflow",

    OutputArray = "electrons",

    DeltaRMax = 0.5,

    PTMin = 0.5,

    PTRatioMax = 0.12,
)

#################
# Muon efficiency
#################

MuonEfficiency = delphes.Module("Efficiency",
    InputArray = "MuonMomentumSmearing/muons",
    OutputArray = "muons",

    # set EfficiencyFormula {efficiency as a function of eta and pt}

    # efficiency formula for muons
    EfficiencyFormula = "                                     (pt <= 10.0)                * (0.00) + \
                                           (abs(eta) <= 1.5) * (pt > 10.0)                * (0.95) + \
                         (abs(eta) > 1.5 && abs(eta) <= 2.4) * (pt > 10.0)                * (0.95) + \
                         (abs(eta) > 2.4)                                                 * (0.00)",
)

################
# Muon isolation
################

MuonIsolation = delphes.Module("Isolation",
    CandidateInputArray = "MuonEfficiency/muons",
    IsolationInputArray = "EFlowFilter/eflow",

    OutputArray = "muons",

    DeltaRMax = 0.5,

    PTMin = 0.5,

    PTRatioMax = 0.25,
)

###################
# Missing ET merger
###################

MissingET = delphes.Module("Merger",
    InputArray = ["EFlowMerger/eflow"],
    MomentumOutputArray = "momentum",
)

##################
# Scalar HT merger
##################

ScalarHT = delphes.Module("Merger",
    InputArray = [
        "UniqueObjectFinder/jets",
        "UniqueObjectFinder/electrons",
        "UniqueObjectFinder/photons",
        "UniqueObjectFinder/muons",
    ],
    EnergyOutputArray = "energy",
)

#####################
# Neutrino Filter
#####################

NeutrinoFilter = delphes.Module("PdgCodeFilter",
    InputArray = "Delphes/stableParticles",
    OutputArray = "filteredParticles",

    PTMin = 0.0,

    PdgCode = [12, 14, 16, -12, -14, -16],
)

#####################
# MC truth jet finder
#####################

GenJetFinder = delphes.Module("FastJetFinder",
    InputArray = "NeutrinoFilter/filteredParticles",

    OutputArray = "jets",

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
    JetAlgorithm = 6,
    ParameterR = 0.5,

    JetPTMin = 20.0,
)

#########################
# Gen Missing ET merger
########################

GenMissingET = delphes.Module("Merger",
    InputArray = ["NeutrinoFilter/filteredParticles"],
    MomentumOutputArray = "momentum",
)

############
# Jet finder
############

FastJetFinder = delphes.Module("FastJetFinder",
    #InputArray = "Calorimeter/towers",
    InputArray = "EFlowMerger/eflow",

    OutputArray = "jets",

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
    JetAlgorithm = 6,
    ParameterR = 0.5,

    JetPTMin = 20.0,
)

##################
# Fat Jet finder
##################

FatJetFinder = delphes.Module("FastJetFinder",
    InputArray = "EFlowMerger/eflow",

    OutputArray = "jets",

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
    JetAlgorithm = 6,
    ParameterR = 0.8,

    ComputeNsubjettiness = True,
    Beta = 1.0,
    AxisMode = 4,

    ComputeTrimming = True,
    RTrim = 0.2,
    PtFracTrim = 0.05,

    ComputePruning = True,
    ZcutPrun = 0.1,
    RcutPrun = 0.5,
    RPrun = 0.8,

    ComputeSoftDrop = True,
    BetaSoftDrop = 0.0,
    SymmetryCutSoftDrop = 0.1,
    R0SoftDrop = 0.8,

    JetPTMin = 200.0,
)

##################
# Jet Energy Scale
##################

JetEnergyScale = delphes.Module("EnergyScale",
    InputArray = "FastJetFinder/jets",
    OutputArray = "jets",

    # scale formula for jets
    ScaleFormula = "sqrt( (2.5 - 0.15*(abs(eta)))^2 / pt + 1.0 )",
)

########################
# Jet Flavor Association
########################

JetFlavorAssociation = delphes.Module("JetFlavorAssociation",
    PartonInputArray = "Delphes/partons",
    ParticleInputArray = "Delphes/allParticles",
    ParticleLHEFInputArray = "Delphes/allParticlesLHEF",
    JetInputArray = "JetEnergyScale/jets",

    DeltaR = 0.5,
    PartonPTMin = 1.0,
    PartonEtaMax = 2.5,
)

###########
# b-tagging
###########

BTagging = delphes.Module("BTagging",
    JetInputArray = "JetEnergyScale/jets",

    BitNumber = 0,

    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
    # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
    # gluon's PDG code has the lowest priority

    # based on arXiv:1211.4462

    EfficiencyFormula = [
        0, "0.01+0.000038*pt",                          # default efficiency formula (misidentification rate)
        4, "0.25*tanh(0.018*pt)*(1/(1+ 0.0013*pt))",    # efficiency formula for c-jets (misidentification rate)
        5, "0.85*tanh(0.0025*pt)*(25.0/(1+0.063*pt))",  # efficiency formula for b-jets
    ],
)

#############
# tau-tagging
#############

TauTagging = delphes.Module("TauTagging",
    ParticleInputArray = "Delphes/allParticles",
    PartonInputArray = "Delphes/partons",
    JetInputArray = "JetEnergyScale/jets",

    DeltaR = 0.5,

    TauPTMin = 1.0,

    TauEtaMax = 2.5,

    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

    EfficiencyFormula = [
        0, "0.01",  # default efficiency formula (misidentification rate)
        15, "0.6",  # efficiency formula for tau-jets
    ],
)

#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

UniqueObjectFinder = delphes.Module("UniqueObjectFinder",
    # earlier arrays take precedence over later ones
    InputArray = [  # InputArray OutputArray
        "PhotonIsolation/photons", "photons",
        "ElectronIsolation/electrons", "electrons",
        "MuonIsolation/muons", "muons",
        "JetEnergyScale/jets", "jets",
    ],
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

        "FatJetFinder/jets", "FatJet", "Jet",

        "MissingET/momentum", "MissingET", "MissingET",
        "ScalarHT/energy", "ScalarHT", "ScalarHT",
    ],
)
