#######################################
# Order of execution of various modules
# 15.9.2014 try a fine granulated calorimeter for LHeC
# 0.025 in eta and phi
# 19.12.16 add paramterisations according to Peter Kostka's design
# muons not considered here - take CMS/FCC values
# without CTagging, with electrons, muons and b-tagging issues considered
# 30.1.2018 hammad add Muon Filter module
#######################################

import Config.Core as delphes
import math


#################################
# Propagate particles in cylinder
#################################

ParticlePropagator = delphes.Module("ParticlePropagator",
    InputArray = "Delphes/stableParticles",
    OutputArray = "stableParticles",
    ChargedHadronOutputArray = "chargedHadrons",
    ElectronOutputArray = "electrons",
    MuonOutputArray = "muons",

    # LHeC
    # radius of the magnetic field coverage, in m
    Radius = 1.14,
    # half-length of the magnetic field coverage, in m
    HalfLength = 2.9,

    # magnetic field
    Bz = 3.8,
)

####################################
# Charged hadron tracking efficiency
####################################

ChargedHadronTrackingEfficiency = delphes.Module("Efficiency",
    InputArray = "ParticlePropagator/chargedHadrons",
    OutputArray = "chargedHadrons",

    # add EfficiencyFormula {efficiency formula as a function of eta and pt}

    # tracking efficiency formula for charged hadrons
    # set to 1 for full range (-0.2 eta already subtracted from edges)
    EfficiencyFormula = "( eta <= 4.9 && eta >= -4.3)  * (1.0) + \
				(eta > 4.9 || eta < -4.3 ) * (0.00)",
)

##############################
# Electron tracking efficiency
##############################

ElectronTrackingEfficiency = delphes.Module("Efficiency",
    InputArray = "ParticlePropagator/electrons",
    OutputArray = "electrons",

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    # tracking efficiency formula for electrons
    # set to 1
    EfficiencyFormula = "( (eta <=4.9 && eta >= -4.3) ) * (1.0) + \
			       ((eta > 4.9 || eta < -4.3)) * (0.00)",
)

##########################
# Muon tracking efficiency
##########################

MuonTrackingEfficiency = delphes.Module("Efficiency",
    InputArray = "ParticlePropagator/muons",
    OutputArray = "muons",

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    # tracking efficiency formula for muons
    # set to 1
    EfficiencyFormula = "(abs(eta) <= 4)  * (1.0) + (abs(eta) > 4) * (0.00)",
)

########################################
# Momentum resolution for charged tracks
########################################

ChargedHadronMomentumSmearing = delphes.Module("MomentumSmearing",
    InputArray = "ChargedHadronTrackingEfficiency/chargedHadrons",
    OutputArray = "chargedHadrons",

    # set ResolutionFormula {resolution formula as a function of eta and pt}

    # resolution formula for charged hadrons
    # taken from FCC but adjust to ep
    #the abs(eta) <= 2.0 should be the "central best" and abs(eta) >
    #2.0  and  <= 4.9 the fwd    /    <= 4.3 the bwd    (LHeC)
    #        and abs(eta) > 2.0 and <= 5.2 fwd   /   <= 5.0   the bwd    (FCC)
    ResolutionFormula = "                  (abs(eta) <= 2.0) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                                           (abs(eta) <= 2.0) * (pt > 1.0   && pt <= 1.0e1) * (0.01) + \
                                           (abs(eta) <= 2.0) * (pt > 1.0e1 && pt <= 2.0e2) * (0.03) + \
                                           (abs(eta) <= 2.0) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 2.0 && ( eta <= 4.9 && eta >= -4.3 ) ) * (pt > 0.1   && pt <= 1.0)   * (0.03) + \
                         (abs(eta) > 2.0 && ( eta <= 4.9 && eta >= -4.3 ) ) * (pt > 1.0   && pt <= 1.0e1) * (0.02) + \
                         (abs(eta) > 2.0 && ( eta <= 4.9 && eta >= -4.3 ) ) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                        (abs(eta) > 2.0 && ( eta <= 4.9 && eta >= -4.3 ) ) * (pt > 2.0e2)                * (0.05)",
)

#################################
# Energy resolution for electrons
#################################

ElectronMomentumSmearing = delphes.Module("MomentumSmearing",
    InputArray = "ElectronTrackingEfficiency/electrons",
    OutputArray = "electrons",

    # set ResolutionFormula {resolution formula as a function of eta and pt}

    # same as resolution formula for charged hadrons
    # taken from FCC but adjust to ep
    #the abs(eta) <= 2.0 should be the "central best" and abs(eta) >
    #2.0  and  <= 4.9 the fwd    /    <= 4.3 the bwd    (LHeC)
    #        and abs(eta) > 2.0 and <= 5.2 fwd   /   <= 5.0   the bwd    (FCC)

    ResolutionFormula = "                  (abs(eta) <= 2.0) * (pt > 0.1   && pt <= 1.0)   * (0.02) + \
                                           (abs(eta) <= 2.0) * (pt > 1.0   && pt <= 1.0e1) * (0.01) + \
                                           (abs(eta) <= 2.0) * (pt > 1.0e1 && pt <= 2.0e2) * (0.03) + \
                                           (abs(eta) <= 2.0) * (pt > 2.0e2)                * (0.05) + \
                         (abs(eta) > 2.0 && ( eta <= 4.9 && eta >= -4.3 ) ) * (pt > 0.1   && pt <= 1.0)   * (0.03) + \
                         (abs(eta) > 2.0 && ( eta <= 4.9 && eta >= -4.3 ) ) * (pt > 1.0   && pt <= 1.0e1) * (0.02) + \
                         (abs(eta) > 2.0 && ( eta <= 4.9 && eta >= -4.3 ) ) * (pt > 1.0e1 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 2.0 && ( eta <= 4.9 && eta >= -4.3 ) ) * (pt > 2.0e2)                * (0.05)",
)

###############################
# Momentum resolution for muons
###############################

MuonMomentumSmearing = delphes.Module("MomentumSmearing",
    InputArray = "MuonTrackingEfficiency/muons",
    OutputArray = "muons",

    # set ResolutionFormula {resolution formula as a function of eta and pt}

    # resolution formula for muons
    # taken from FCC
    ResolutionFormula = "                  (abs(eta) <= 0.5) * (pt > 0.1   && pt <= 5.0)   * (0.02) + \
                                           (abs(eta) <= 0.5) * (pt > 5.0   && pt <= 1.0e2) * (0.015) + \
                                           (abs(eta) <= 0.5) * (pt > 1.0e2 && pt <= 2.0e2) * (0.03) + \
                                           (abs(eta) <= 0.5) * (pt > 2.0e2)                * (0.05 + pt*1.e-4) + \
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 0.1   && pt <= 5.0)   * (0.03) + \
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 5.0   && pt <= 1.0e2) * (0.02) + \
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 1.0e2 && pt <= 2.0e2) * (0.04) + \
                         (abs(eta) > 0.5 && abs(eta) <= 1.5) * (pt > 2.0e2)                * (0.05 + pt*1.e-4) + \
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 0.1   && pt <= 5.0)   * (0.04) + \
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 5.0   && pt <= 1.0e2) * (0.035) + \
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 1.0e2 && pt <= 2.0e2) * (0.05) + \
                         (abs(eta) > 1.5 && abs(eta) <= 4.0) * (pt > 2.0e2)                * (0.05 + pt*1.e-4)",
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

##################
# LHeC Calorimeter
##################

Calorimeter = delphes.Module("Calorimeter",
    ParticleInputArray = "ParticlePropagator/stableParticles",
    TrackInputArray = "TrackMerger/tracks",

    TowerOutputArray = "towers",
    PhotonOutputArray = "photons",

    EFlowTrackOutputArray = "eflowTracks",
    EFlowPhotonOutputArray = "eflowPhotons",
    EFlowNeutralHadronOutputArray = "eflowNeutralHadrons",

    # lists of the edges of each tower in eta and phi
    # each list starts with the lower edge of the first tower
    # the list ends with the higher edged of the last tower

    # use a delta_phi=0.025 for 2*phi range == 252 cells
    # use a delta_eta-0.025 for eta = -5 to 5 == 400 cells
    # Uta 22.12.16 use eta range from Peter Kostka -4.7 - 5.2
    PhiBins = [i * math.pi / 126. for i in range(-126, 127)],
    EtaPhiBins = [],

    EnergyFraction = [
        0, [0.0, 1.0],  # default energy fractions {abs(PDG code)} {Fecal Fhcal}
        11, [1.0, 0.0], # energy fractions for e, gamma and pi0
        22, [1.0, 0.0],
        111, [1.0, 0.0],
        12, [0.0, 0.0],  # energy fractions for muon, neutrinos and neutralinos
        13, [0.0, 0.0],
        14, [0.0, 0.0],
        16, [0.0, 0.0],
        1000022, [0.0, 0.0],
        1000023, [0.0, 0.0],
        1000025, [0.0, 0.0],
        1000035, [0.0, 0.0],
        1000045, [0.0, 0.0],
        310, [0.3, 0.7],  # energy fractions for K0short and Lambda
        3122, [0.3, 0.7],
    ],

    ## no minimum p range given (Uta 22.12.16)
    ## sigma_E = E* (a/sqrt(E) + b) = sqrt ( E*a^2 E +E^2*b^2 )
    ## a = sampling term, b = constant term
    # ECAL eta : 5.0 - 2.7 - -2.1 - -4.4 (subtract -0.1 in eta as safety distance to edges)
    ECalResolutionFormula = "( eta <= 5.0 && eta >   2.7 ) * sqrt(energy*0.1^2  + energy^2*0.01^2) + \
                             ( eta <= 2.7 && eta >= -2.1 ) * sqrt(energy*0.09^2 + energy^2*0.02^2) + \
                             ( eta < -2.1 && eta >= -4.4 ) * sqrt(energy*0.1^2  + energy^2*0.01^2)",

    #### HCAL eta : 5.1 - 2.1 - -1.7 - -4.6 (subtract -0.1 in eta as safety distance to edges)
    #HCalResolutionFormula = resolution formula as a function of eta and energy

    HCalResolutionFormula = "( eta <= 5.1 && eta >   2.1 ) * sqrt(energy*0.4^2 + energy^2*0.02^2) + \
                             ( eta <= 2.1 && eta >= -1.7 ) * sqrt(energy*0.4^2 + energy^2*0.02^2) + \
                             ( eta < -1.7 && eta >= -4.7 ) * sqrt(energy*0.4^2 + energy^2*0.04^2)",
)

for eta in [-4.700, -4.675, -4.650, -4.625, -4.600, -4.575, -4.550, -4.525,
            -4.500, -4.475, -4.450, -4.425, -4.400, -4.375, -4.350, -4.325, -4.300, -4.275, -4.250, -4.225, -4.200, -4.175, -4.150, -4.125, -4.100, -4.075, -4.050, -4.025,
            -4.000, -3.975, -3.950, -3.925, -3.900, -3.875, -3.850, -3.825, -3.800, -3.775, -3.750, -3.725, -3.700, -3.675, -3.650, -3.625, -3.600, -3.575, -3.550, -3.525,
            -3.500, -3.475, -3.450, -3.425, -3.400, -3.375, -3.350, -3.325, -3.300, -3.275, -3.250, -3.225, -3.200, -3.175, -3.150, -3.125, -3.100, -3.075, -3.050, -3.025,
            -3.000, -2.975, -2.950, -2.925, -2.900, -2.875, -2.850, -2.825, -2.800, -2.775, -2.750, -2.725, -2.700, -2.675, -2.650, -2.625, -2.600, -2.575, -2.550, -2.525,
            -2.500, -2.475, -2.450, -2.425, -2.400, -2.375, -2.350, -2.325, -2.300, -2.275, -2.250, -2.225, -2.200, -2.175, -2.150, -2.125, -2.100, -2.075, -2.050, -2.025,
            -2.000, -1.975, -1.950, -1.925, -1.900, -1.875, -1.850, -1.825, -1.800, -1.775, -1.750, -1.725, -1.700, -1.675, -1.650, -1.625, -1.600, -1.575, -1.550, -1.525,
            -1.500, -1.475, -1.450, -1.425, -1.400, -1.375, -1.350, -1.325, -1.300, -1.275, -1.250, -1.225, -1.200, -1.175, -1.150, -1.125, -1.100, -1.075, -1.050, -1.025,
            -1.000, -0.975, -0.950, -0.925, -0.900, -0.875, -0.850, -0.825, -0.800, -0.775, -0.750, -0.725, -0.700, -0.675, -0.650, -0.625, -0.600, -0.575, -0.550, -0.525,
            -0.500, -0.475, -0.450, -0.425, -0.400, -0.375, -0.350, -0.325, -0.300, -0.275, -0.250, -0.225, -0.200, -0.175, -0.150, -0.125, -0.100, -0.075, -0.050, -0.025,
             0.000,  0.025,  0.050,  0.075,  0.100,  0.125,  0.150,  0.175,  0.200,  0.225,  0.250,  0.275,  0.300,  0.325,  0.350,  0.375,  0.400,  0.425,  0.450,  0.475,
             0.500,  0.525,  0.550,  0.575,  0.600,  0.625,  0.650,  0.675,  0.700,  0.725,  0.750,  0.775,  0.800,  0.825,  0.850,  0.875,  0.900,  0.925,  0.950,  0.975,
             1.000,  1.025,  1.050,  1.075,  1.100,  1.125,  1.150,  1.175,  1.200,  1.225,  1.250,  1.275,  1.300,  1.325,  1.350,  1.375,  1.400,  1.425,  1.450,  1.475,
             1.500,  1.525,  1.550,  1.575,  1.600,  1.625,  1.650,  1.675,  1.700,  1.725,  1.750,  1.775,  1.800,  1.825,  1.850,  1.875,  1.900,  1.925,  1.950,  1.975,
             2.000,  2.025,  2.050,  2.075,  2.100,  2.125,  2.150,  2.175,  2.200,  2.225,  2.250,  2.275,  2.300,  2.325,  2.350,  2.375,  2.400,  2.425,  2.450,  2.475,
             2.500,  2.525,  2.550,  2.575,  2.600,  2.625,  2.650,  2.675,  2.700,  2.725,  2.750,  2.775,  2.800,  2.825,  2.850,  2.875,  2.900,  2.925,  2.950,  2.975,
             3.000,  3.025,  3.050,  3.075,  3.100,  3.125,  3.150,  3.175,  3.200,  3.225,  3.250,  3.275,  3.300,  3.325,  3.350,  3.375,  3.400,  3.425,  3.450,  3.475,
             3.500,  3.525,  3.550,  3.575,  3.600,  3.625,  3.650,  3.675,  3.700,  3.725,  3.750,  3.775,  3.800,  3.825,  3.850,  3.875,  3.900,  3.925,  3.950,  3.975,
             4.000,  4.025,  4.050,  4.075,  4.100,  4.125,  4.150,  4.175,  4.200,  4.225,  4.250,  4.275,  4.300,  4.325,  4.350,  4.375,  4.400,  4.425,  4.450,  4.475,
             4.500,  4.525,  4.550,  4.575,  4.600,  4.625,  4.650,  4.675,  4.700,  4.725,  4.750,  4.775,  4.800,  4.825,  4.850,  4.875,  4.900,  4.925,  4.950,  4.975,
             5.000,  5.025,  5.050,  5.075,  5.100,  5.125,  5.150,  5.175,  5.200]:
    Calorimeter['EtaPhiBins'].append([eta,])
    Calorimeter['EtaPhiBins'].append(Calorimeter['PhiBins'])

####################
# Energy flow merger
####################

EFlowMerger = delphes.Module("Merger",
    # add InputArray InputArray
    # exchange UK:  25.7.16 change back to all eflowtracks!
    #  add InputArray ImpactParameterSmearing/tracks
    ## eflowTracks kills electrons for Delphes 3.3.3 and ep??
    InputArray = [
        "Calorimeter/eflowTracks",
        "Calorimeter/eflowPhotons",
        "Calorimeter/eflowNeutralHadrons",
    ],
    OutputArray = "eflow",
)

###################
# Photon efficiency
###################

PhotonEfficiency = delphes.Module("Efficiency",
    #  set InputArray Calorimeter/eflowPhotons
    # 21.7.16 use photons
    InputArray = "Calorimeter/photons",
    OutputArray = "photons",

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    # efficiency formula for photons - same as for electrons
    # expand ad set flat
    ## OLD EFFICIENCY
    ## EfficiencyFormula = "( (eta <=4.9 && eta >= -4.3)) * (1.0) + \
    ##	                     ((eta > 4.9 || eta < -4.3))  * (0.00)",
    #########################################################
    # OF, 20.03.24: NEW EFFICIENCY TO AVOID SOFT PHOTONS   ##
    #########################################################
    EfficiencyFormula = "   (pt <= 0.5) * (0.00) + \
                            (pt > 0.5) *( (eta <=4.9 && eta >= -4.3) ) * (1.0) + \
                            ((eta > 4.9 || eta < -4.3)) * (0.00)",
)

##################
# Photon isolation
##################

PhotonIsolation = delphes.Module("Isolation",
    CandidateInputArray = "PhotonEfficiency/photons",
    IsolationInputArray = "EFlowMerger/eflow",
    OutputArray = "photons",
    DeltaRMax = 0.4,
    PTMin = 0.5,
    PTRatioMax = 0.1,
)

#################
# Electron filter
#################

ElectronFilter = delphes.Module("PdgCodeFilter",
    InputArray = "Calorimeter/eflowTracks",
    OutputArray = "electrons",
    Invert = True,
    PdgCode = [11, -11],
)

#################
# Muon filter
#################

MuonFilter = delphes.Module("PdgCodeFilter",
    InputArray = "Calorimeter/eflowTracks",
    OutputArray = "muons",
    Invert = True,
    PdgCode = [13, -13],
)

#####################
# Electron efficiency
#####################

ElectronEfficiency = delphes.Module("Efficiency",
    InputArray = "ElectronFilter/electrons",
    OutputArray = "electrons",

    # set EfficiencyFormula {efficiency formula as a function of eta and pt}

    # efficiency formula for electrons
    #  1 degree =  eta<4.7
    ## OLD EFFICIENCY
    ## EfficiencyFormula = "(eta <=4.9 && eta >= -4.3)) * (1.0) + \
    ##                     ((eta > 4.9 || eta < -4.3))  * (0.00)",
    #########################################################
    # OF, 20.03.24: NEW EFFICIENCY TO AVOID SOFT ELECTRONS ##
    #########################################################
    EfficiencyFormula = "(pt <= 0.5) * (0.00) + \
                         (pt >  0.5) * ((eta <=4.9 && eta >= -4.3) * (1.0) + \
                                        (eta > 4.9 || eta <  -4.3) * (0.0))",
)

####################
# Electron isolation
####################

ElectronIsolation = delphes.Module("Isolation",
    CandidateInputArray = "ElectronEfficiency/electrons",
    IsolationInputArray = "EFlowMerger/eflow",
    OutputArray = "electrons",
    DeltaRMax = 0.4,
    PTMin = 0.5,
    PTRatioMax = 0.1,
)

#################
# Muon efficiency
#################

MuonEfficiency = delphes.Module("Efficiency",
    #InputArray = "MuonMomentumSmearing/muons",  #FIXME correct?
    InputArray = "MuonFilter/muons",
    OutputArray = "muons",

    # set EfficiencyFormula {efficiency as a function of eta and pt}

    # efficiency formula for muons
    # set to 1 for |eta|<4
    EfficiencyFormula = "(abs(eta) <= 4)  * (1.0) + \
                         (abs(eta) > 4)   * (0.00)",
)

################
# Muon isolation
################

MuonIsolation = delphes.Module("Isolation",
    CandidateInputArray = "MuonEfficiency/muons",
    IsolationInputArray = "EFlowMerger/eflow",
    OutputArray = "muons",
    DeltaRMax = 0.4,
    PTMin = 0.5,
    PTRatioMax = 0.1,
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
# MC truth jet finder
#####################

GenJetFinder = delphes.Module("FastJetFinder",
    InputArray = "Delphes/stableParticles",
    OutputArray = "jets",

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
    JetAlgorithm = 6,
    ParameterR = 0.7,
    JetPTMin = 3.0,
)

############
# Jet finder
############

FastJetFinder = dict(
    type ="FastJetFinder",
    # ATLAS uses towers
    # InputArray = "Calorimeter/towers",
    # CMS uses eflow?
    InputArray = "EFlowMerger/eflow",
    OutputArray = "jets",

    # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
    JetAlgorithm = 6,
    ParameterR = 0.4,
    JetPTMin = 1.0,
)

##################
# Jet Energy Scale
##################

JetEnergyScale = delphes.Module("EnergyScale",
    InputArray = "FastJetFinder/jets",
    OutputArray = "jets",
    ScaleFormula = "1.00",  # scale formula for jets
)

########################
# GenJet Flavor Association
########################

GenJetFlavorAssociation = delphes.Module("JetFlavorAssociation",
    PartonInputArray = "Delphes/partons",
    ParticleInputArray = "Delphes/allParticles",
    ParticleLHEFInputArray = "Delphes/allParticlesLHEF",
    JetInputArray = "GenJetFinder/jets",

    DeltaR = 0.4,
    PartonPTMin = 0.5,
    PartonEtaMax = 3.0,
)

########################
# Jet Flavor Association
########################

JetFlavorAssociation = delphes.Module("JetFlavorAssociation",
    PartonInputArray = "Delphes/partons",
    ParticleInputArray = "Delphes/allParticles",
    ParticleLHEFInputArray = "Delphes/allParticlesLHEF",
    JetInputArray = "JetEnergyScale/jets",

    DeltaR = 0.4,
    PartonPTMin = 0.5,
    PartonEtaMax = 3.0,
)

###########
# b-tagging
###########
GenBTagging = delphes.Module("BTagging",
    PartonInputArray = "Delphes/partons",
    JetInputArray = "GenJetFinder/jets",
    BitNumber = 0,
    DeltaR = 0.4,
    PartonPTMin = 0.5,
    PartonEtaMax = 6.0,

    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
    # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
    # gluon's PDG code has the lowest priority

    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsBTV
    # default efficiency formula (misidentification rate)
    EfficiencyFormula = [
        0, "0.001",
        4, "0.05",  # efficiency formula for c-jets (misidentification rate)
        5, "0.75",  # efficiency formula for b-jets
    ],
)

####################
# standard B-tagging
## ATTENTION : idealised values - use for cross check only
####################
BTagging = delphes.Module("BTagging",
    PartonInputArray = "Delphes/partons",
    JetInputArray = "JetEnergyScale/jets",
    BitNumber = 0,
    DeltaR = 0.4,
    PartonPTMin = 0.5,
    PartonEtaMax = 6.0,

    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}
    # PDG code = the highest PDG code of a quark or gluon inside DeltaR cone around jet axis
    # gluon's PDG code has the lowest priority

    # https://twiki.cern.ch/twiki/bin/view/CMSPublic/PhysicsResultsBTV
    # default efficiency formula (misidentification rate)
    EfficiencyFormula = [
        0, "0.001",
        4, "0.05",  # efficiency formula for c-jets (misidentification rate)
        5, "0.75",  # efficiency formula for b-jets
    ],
)

TauTagging = delphes.Module("TauTagging",
    ParticleInputArray = "Delphes/allParticles",
    PartonInputArray = "Delphes/partons",
    JetInputArray = "JetEnergyScale/jets",

    DeltaR = 0.4,
    TauPTMin = 0.5,
    TauEtaMax = 6.0,

    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

    # default efficiency formula (misidentification rate)
    EfficiencyFormula = [
        0, "0.001",
        15, "0.4",  # efficiency formula for tau-jets
    ],
)

GenTauTagging = delphes.Module("TauTagging",
    ParticleInputArray = "Delphes/allParticles",
    PartonInputArray = "Delphes/partons",
    JetInputArray = "GenJetFinder/jets",

    DeltaR = 0.4,
    TauPTMin = 0.5,
    TauEtaMax = 6.0,

    # add EfficiencyFormula {abs(PDG code)} {efficiency formula as a function of eta and pt}

    # default efficiency formula (misidentification rate)
    EfficiencyFormula = [
        0, "0.001",
        15, "0.4",  # efficiency formula for tau-jets
    ],
)

#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

UniqueObjectFinder = delphes.Module("UniqueObjectFinder",
    # earlier arrays take precedence over later ones
    # add InputArray InputArray OutputArray
    InputArray = [
        "PhotonIsolation/photons", "photons",
        "ElectronIsolation/electrons", "electrons",
        "MuonIsolation/muons", "muons",
        "JetEnergyScale/jets", "jets",
    ]
)

##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

TreeWriter = delphes.Module("TreeWriter",
    # add Branch InputArray BranchName BranchClass
    Branch = [
        "Delphes/allParticles", "Particle", "GenParticle",
        #"TrackMerger/tracks", ["Track", "Track"],
        #"Calorimeter/towers", ["Tower", "Tower"],
        #"Calorimeter/eflowTracks", ["EFlowTrack", "Track"],
        #"Calorimeter/eflowPhotons", [EFlowPhoton", "Tower"],
        #"Calorimeter/eflowNeutralHadrons", ["EFlowNeutralHadron", "Tower"],
        "GenJetFinder/jets", "GenJet", "Jet",
        "UniqueObjectFinder/jets", "Jet", "Jet",
        "UniqueObjectFinder/electrons", "Electron", "Electron",
        "UniqueObjectFinder/photons", "Photon", "Photon",
        "UniqueObjectFinder/muons", "Muon", "Muon",
        "MissingET/momentum", "MissingET", "MissingET",
        "ScalarHT/energy", "ScalarHT", "ScalarHT",
    ]
)

ExecutionPath = [
    "ParticlePropagator",

    "ChargedHadronTrackingEfficiency",
    "ElectronTrackingEfficiency",
    "MuonTrackingEfficiency",

    "ChargedHadronMomentumSmearing",
    "ElectronMomentumSmearing",
    "MuonMomentumSmearing",

    "TrackMerger",
    "Calorimeter",
    "EFlowMerger",

    "PhotonEfficiency",
    "PhotonIsolation",

    "ElectronFilter",
    "MuonFilter",
    "ElectronEfficiency",
    "ElectronIsolation",

    "MuonEfficiency",
    "MuonIsolation",

    "MissingET",

    "GenJetFinder",
    "FastJetFinder",

    "JetEnergyScale",

    "GenJetFlavorAssociation",
    "JetFlavorAssociation",

    "GenBTagging",
    "BTagging",

    "GenTauTagging",
    "TauTagging",

    "UniqueObjectFinder",

    "ScalarHT",

    "TreeWriter",
    #"TrackCountingBTagging",  # exclude since it adds also an entry to Jet->BTag
]
