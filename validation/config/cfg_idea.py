import os
from collections import OrderedDict
from src.utils import (
    compute_bins,
    Observable,
    Particle,
    Slices1D,
    EfficiencyParticleBlock,
    EfficiencyTaggingBlock,
    ResolutionBlock,
    ResolutionHisto,
    ResolutionPlot,
    Text,
)

## provide relative path to the directory
card = "validation_delphes_card_IDEA.tcl"
name = os.path.basename(card.replace(".tcl", ""))

## run parameters:
njobs = 100
nevts_per_job = 25000

#
pmin = 0.1
pmax = 200.0
log = True
etamin = -3.0
etamax = 3.0
ecm = 400.0

# plotting parameters
nbins = 100

nbins_plin = nbins
plin_min = 0.0
plin_min = 0.0
plin_max = pmax

nbins_plog = nbins
plog_min = pmin
plog_max = pmax

nbins_elin = nbins
elin_min = 0.0
elin_max = pmax

nbins_elog = nbins
elog_min = pmin
elog_max = pmax

nbins_theta = nbins
theta_min = 0.1
theta_max = 3.14

nbins_eta = nbins
eta_min = -3.0
eta_max = 3.0

bins_thetalinf = compute_bins(theta_min, theta_max, nbins_theta, 2, "lin")

bins_plinf = compute_bins(plin_min, plin_max, nbins_plin, 2, "lin")
bins_plogf = compute_bins(plog_min, plog_max, nbins_plog, 2, "log")

bins_jetplinf = compute_bins(20.0, plin_max, nbins_plin, 2, "lin")
bins_jetplogf = compute_bins(20.0, plog_max, nbins_plog, 2, "log")

bins_elinf = compute_bins(elin_min, elin_max, nbins_elin, 2, "lin")
bins_elogf = compute_bins(elog_min, elog_max, nbins_elog, 2, "log")

bins_jetelinf = compute_bins(20.0, elin_max, nbins_elin, 2, "lin")
bins_jetelogf = compute_bins(20.0, elog_max, nbins_elog, 2, "log")

bins_thetalinc = [(1.4, 1.8), (0.6, 1.0), (0.1, 0.5)]
bins_plogc = [
    (0.5, 1.5),
    (5, 15),
    (20, 50),
]

bins_jetplogc = [
    (20, 50),
    (50, 100),
    (100, 500),
]

bins_elogc = bins_plogc
bins_jetelogc = bins_jetplogc

# observables (name, varname, label, nbins, xmin, xmax, opt)
# last 4 params only used if obs in ResolutionHisto
obs_pt = Observable("pt", "PT", "p_{T}", "[GeV]", 500, 0.8, 1.2, "rel")
obs_p = Observable("p", "P", "p", "[GeV]", 500, 0.9, 1.1, "rel")
obs_e = Observable("e", "E", "E", "[GeV]", 500, 0.5, 1.5, "rel")
obs_d0 = Observable("d0", "D0", "d_0", "[mm]", 500, -0.2, 0.2, "abs")
obs_dz = Observable("dz", "DZ", "d_z", "[mm]", 500, -0.2, 0.2, "abs")
obs_phi = Observable("phi", "Phi", "\phi", " [rad]", 500, -0.001, 0.001, "abs")
obs_theta = Observable("theta", "Theta", "\\theta", "[rad]", 500, -0.001, 0.001, "abs")

## covariance matrix
obs_corr_d0phi = Observable("corr_d0phi", "CorrD0Phi", "\\rho(d_0, \phi)", "", 100, -1.0, 1.0, "mean")
obs_corr_d0c = Observable("corr_d0c", "CorrD0C", "\\rho(d_0, C)", "", 100, -1.0, 1.0, "mean")
obs_corr_d0dz = Observable("corr_d0dz", "CorrD0DZ", "\\rho(d_0, d_z)", "", 100, -1.0, 1.0, "mean")
obs_corr_d0ctgtheta = Observable(
    "corr_d0ctgtheta", "CorrD0CtgTheta", "\\rho(d_0, ctg(\\theta))", "", 100, -1.0, 1.0, "mean"
)
obs_corr_phic = Observable("corr_phic", "CorrPhiC", "\\rho(\phi, C)", "", 100, -1.0, 1.0, "mean")
obs_corr_phidz = Observable("corr_phidz", "CorrPhiDZ", "\\rho(\phi, d_z)", "", 100, -1.0, 1.0, "mean")
obs_corr_phictgtheta = Observable(
    "corr_phictgtheta", "CorrPhiCtgTheta", "\\rho(\phi, ctg(\\theta))", "", 100, -1.0, 1.0, "mean"
)
obs_corr_cdz = Observable("corr_cdz", "CorrCDZ", "\\rho(C, d_z)", "", 100, -1.0, 1.0, "mean")
obs_corr_cctgtheta = Observable("corr_cctgtheta", "CorrCCtgTheta", "\\rho(C, ctg(\\theta))", "", 100, -1.0, 1.0, "mean")
obs_corr_dzctgtheta = Observable(
    "corr_dzctgtheta", "CorrDZCtgTheta", "\\rho(d_z, ctg(\\theta))", "", 100, -1.0, 1.0, "mean"
)


# particles
pion = Particle("pi", "$\pi^{\pm}$", 211)
electron = Particle("ele", "$e^{\pm}$", 11)
muon = Particle("mu", "$\mu^{\pm}$", 13)
tau = Particle("tau", "$\tau$", 15)
photon = Particle("gam", "$\gamma$", 22)
kl = Particle("kl", "$K_{L}$", 130)
up = Particle("up", "u", 1)
down = Particle("down", "d", 2)
strange = Particle("strange", "s", 3)
charm = Particle("charm", "c", 4)
bottom = Particle("bottom", "b", 5)
gluon = Particle("gluon", "g", 21)


plogf = Slices1D(obs_p, bins_plogf, "log")
pjetlogf = Slices1D(obs_p, bins_jetplogf, "log")

elogf = Slices1D(obs_e, bins_elogf, "log")
ejetlogf = Slices1D(obs_e, bins_jetelogf, "log")
thetalinf = Slices1D(obs_theta, bins_thetalinf, "linear")

plogc = Slices1D(obs_p, bins_plogc, "log")
pjetlogc = Slices1D(obs_p, bins_jetplogc, "log")
elogc = Slices1D(obs_e, bins_elogc, "log")
thetalinc = Slices1D(obs_theta, bins_thetalinc, "linear")

## plots ##
reso_plots = []
eff_plots = []  ## all eff plots (for step1)
eff_tag_plots = []  ## all eff plots (for step1)
report = OrderedDict()


########################
# Tracking
########################


report["tracks"] = OrderedDict()

# resolution
mom = (plogf, plogc)
theta = (thetalinf, thetalinc)

obs = [
    obs_p,
    obs_d0,
    obs_dz,
    obs_phi,
    obs_theta,
    obs_corr_d0phi,
    obs_corr_d0c,
    obs_corr_d0dz,
    obs_corr_d0ctgtheta,
    obs_corr_phic,
    obs_corr_phidz,
    obs_corr_phictgtheta,
    obs_corr_cdz,
    obs_corr_cctgtheta,
    obs_corr_dzctgtheta,
]

tracks = [pion, muon, electron]
# tracks = [pion]
reso_track_plots = ResolutionBlock(tracks, "Track", obs, mom, theta, reso_plots)

report["tracks"]["resolution"] = OrderedDict()
for p in tracks:
    for o in obs:
        title = "{} track resolution: ${}$".format(p.label, o.label)
        report["tracks"]["resolution"][title] = [
            reso_track_plots[((p, o, "pt"))],
            reso_track_plots[((p, o, "eta"))],
        ]


# efficiency
mom = (plogf, plogc)
theta = (thetalinf, thetalinc)
tracks = [pion, muon, electron]
# tracks = [pion]
eff_track_plots = EfficiencyParticleBlock(tracks, "Track", mom, theta, eff_plots)

report["tracks"]["efficiency"] = OrderedDict()
for p in tracks:
    title = "{} track efficiency".format(p.label)
    report["tracks"]["efficiency"][title] = [
        eff_track_plots[((p, "pt"))],
        eff_track_plots[((p, "eta"))],
    ]

########################
# Calorimeter particles
########################

report["calorimeter"] = OrderedDict()
caloparticles = [electron, photon, pion, kl]

# resolution
mom = (elogf, elogc)
theta = (thetalinf, thetalinc)
obs = [obs_e, obs_theta, obs_phi]

caloparticles = [electron, photon, pion, kl]
reso_calo_plots = ResolutionBlock(caloparticles, "Tower", obs, mom, theta, reso_plots)

report["calorimeter"]["resolution"] = OrderedDict()
for p in caloparticles:
    for o in obs:
        title = "{} calorimeter resolution: ${}$".format(p.label, o.label)
        report["calorimeter"]["resolution"][title] = [
            reso_calo_plots[((p, o, "pt"))],
            reso_calo_plots[((p, o, "eta"))],
        ]

# efficiency
mom = (plogf, plogc)
theta = (thetalinf, thetalinc)
caloparticles = [electron, photon, pion, kl]
eff_calo_plots = EfficiencyParticleBlock(caloparticles, "Tower", mom, theta, eff_plots)

report["calorimeter"]["efficiency"] = OrderedDict()
for p in caloparticles:
    title = "{} calorimeter efficiency".format(p.label)
    report["calorimeter"]["efficiency"][title] = [
        eff_calo_plots[((p, "pt"))],
        eff_calo_plots[((p, "eta"))],
    ]


########################
# PF particles
########################


report["particle-flow"] = OrderedDict()
report["particle-flow"]["resolution"] = OrderedDict()


pftracks = [electron, pion]
pfgamma = [photon]
pfnh = [kl]

# resolution
mom = (elogf, elogc)
theta = (thetalinf, thetalinc)
obs = [obs_e, obs_theta, obs_phi]

reso_pftrack_plots = ResolutionBlock(pftracks, "EFlowTrack", obs, mom, theta, reso_plots)
reso_pfgamma_plots = ResolutionBlock(pfgamma, "EFlowPhoton", obs, mom, theta, reso_plots)
reso_pfnh_plots = ResolutionBlock(pfnh, "EFlowNeutralHadron", obs, mom, theta, reso_plots)


for p in pftracks:
    for o in obs:
        title = "{} particle-flow resolution: ${}$".format(p.label, o.label)
        report["particle-flow"]["resolution"][title] = [
            reso_pftrack_plots[((p, o, "pt"))],
            reso_pftrack_plots[((p, o, "eta"))],
        ]
for p in pfgamma:
    for o in obs:
        title = "{} particle-flow resolution: ${}$".format(p.label, o.label)
        report["particle-flow"]["resolution"][title] = [
            reso_pfgamma_plots[((p, o, "pt"))],
            reso_pfgamma_plots[((p, o, "eta"))],
        ]
for p in pfnh:
    for o in obs:
        title = "{} particle-flow resolution: ${}$".format(p.label, o.label)
        report["particle-flow"]["resolution"][title] = [
            reso_pfnh_plots[((p, o, "pt"))],
            reso_pfnh_plots[((p, o, "eta"))],
        ]


pftest_parts = [electron, pion]
reso_pftest_dict = OrderedDict()
collection_name_dict = {"Track": "track", "Tower": "calo", "EFlowTrack": "particle-flow"}

reso_eletrk_e = ResolutionHisto(
    electron,
    "Track",
    obs_e,
    elogf,
    {obs_theta: (0.73, 1.57)},
    "track, 0 < $\\theta$ < 1.57",
)

reso_elecalo_e = ResolutionHisto(
    electron,
    "Tower",
    obs_e,
    elogf,
    {obs_theta: (0.73, 1.57)},
    "calo, 0 < $\\theta$ < 1.57",
)

reso_elepf_e = ResolutionHisto(
    electron,
    "EFlowTrack",
    obs_e,
    elogf,
    {obs_theta: (0.73, 1.57)},
    "particle-flow, 0 < $\\theta$ < 1.57",
)

reso_ele_pftest_histos = [reso_eletrk_e, reso_elecalo_e, reso_elepf_e]
reso_ele_pftest_plot = ResolutionPlot("reso_ele_pftest_e", reso_ele_pftest_histos, Text("electrons", (0.5, 0.5)))

reso_pitrk_e = ResolutionHisto(
    pion,
    "Track",
    obs_e,
    elogf,
    {obs_theta: (0.73, 1.57)},
    "track, 0.73 < $\\theta$ < 1.57",
)

reso_picalo_e = ResolutionHisto(
    pion,
    "Tower",
    obs_e,
    elogf,
    {obs_theta: (0.73, 1.57)},
    "calo, 0.73 < $\\theta$ < 1.57",
)

reso_pipf_e = ResolutionHisto(
    pion,
    "EFlowTrack",
    obs_e,
    elogf,
    {obs_theta: (0.73, 1.57)},
    "particle-flow, 0.73 < $\\theta$ < 1.57",
)

reso_pi_pftest_histos = [reso_pitrk_e, reso_picalo_e, reso_pipf_e]
reso_pi_pftest_plot = ResolutionPlot("reso_pi_pftest_e", reso_pi_pftest_histos, Text("charged pions", (0.5, 0.5)))


reso_jetcalo_e = ResolutionHisto(
    gluon,
    "CaloJet",
    obs_e,
    ejetlogf,
    {obs_theta: (0.73, 1.57)},
    "calo, 0.73 < $\\theta$ < 1.57",
)

reso_jetpf_e = ResolutionHisto(
    gluon,
    "PFJet",
    obs_e,
    ejetlogf,
    {obs_theta: (0.73, 1.57)},
    "particle-flow, 0.73 < $\\theta$ < 1.57",
)

reso_jet_pftest_histos = [reso_jetcalo_e, reso_jetpf_e]
reso_jet_pftest_plot = ResolutionPlot("reso_jet_pftest_e", reso_jet_pftest_histos, Text("PF-jets", (0.5, 0.5)))

title = "particle-flow/calo/track resolution: ${}$".format(obs_e.label)
report["particle-flow"]["resolution"][title] = [
    reso_ele_pftest_plot,
    reso_pi_pftest_plot,
    reso_jet_pftest_plot,
]

reso_plots.append(reso_ele_pftest_plot)
reso_plots.append(reso_pi_pftest_plot)
reso_plots.append(reso_jet_pftest_plot)


# efficiency
mom = (elogf, elogc)
theta = (thetalinf, thetalinc)

eff_pftrack_plots = EfficiencyParticleBlock(pftracks, "EFlowTrack", mom, theta, eff_plots)
eff_pfgamma_plots = EfficiencyParticleBlock(pfgamma, "EFlowPhoton", mom, theta, eff_plots)
eff_pfnh_plots = EfficiencyParticleBlock(pfnh, "EFlowNeutralHadron", mom, theta, eff_plots)

report["particle-flow"]["efficiency"] = OrderedDict()
for p in pftracks:
    title = "{} particle-flow efficiency: ${}$".format(p.label, o.label)
    report["particle-flow"]["efficiency"][title] = [
        eff_pftrack_plots[((p, "pt"))],
        eff_pftrack_plots[((p, "eta"))],
    ]
for p in pfgamma:
    title = "{} particle-flow efficiency: ${}$".format(p.label, o.label)
    report["particle-flow"]["efficiency"][title] = [
        eff_pfgamma_plots[((p, "pt"))],
        eff_pfgamma_plots[((p, "eta"))],
    ]
for p in pfnh:
    title = "{} particle-flow efficiency: ${}$".format(p.label, o.label)
    report["particle-flow"]["efficiency"][title] = [
        eff_pfnh_plots[((p, "pt"))],
        eff_pfnh_plots[((p, "eta"))],
    ]

########################
# Electron
########################

report["electron"] = OrderedDict()

# efficiency
mom = (plogf, plogc)
theta = (thetalinf, thetalinc)
coll = [electron]
eff_electron_plots = EfficiencyParticleBlock(coll, "Electron", mom, theta, eff_plots)

report["electron"]["efficiency"] = OrderedDict()
for p in coll:
    title = "{} electron efficiency".format(p.label)
    report["electron"]["efficiency"][title] = [
        eff_electron_plots[((p, "pt"))],
        eff_electron_plots[((p, "eta"))],
    ]


########################
# Muon
########################

report["muon"] = OrderedDict()

# efficiency
mom = (plogf, plogc)
theta = (thetalinf, thetalinc)
coll = [muon]
eff_muon_plots = EfficiencyParticleBlock(coll, "Muon", mom, theta, eff_plots)

report["muon"]["efficiency"] = OrderedDict()
for p in coll:
    title = "{} muon efficiency".format(p.label)
    report["muon"]["efficiency"][title] = [
        eff_muon_plots[((p, "pt"))],
        eff_muon_plots[((p, "eta"))],
    ]

########################
# Photon
########################

report["photon"] = OrderedDict()

# efficiency
mom = (plogf, plogc)
theta = (thetalinf, thetalinc)
coll = [photon]
eff_photon_plots = EfficiencyParticleBlock(coll, "Photon", mom, theta, eff_plots)

report["photon"]["efficiency"] = OrderedDict()
for p in coll:
    title = "{} photon efficiency".format(p.label)
    report["photon"]["efficiency"][title] = [
        eff_photon_plots[((p, "pt"))],
        eff_photon_plots[((p, "eta"))],
    ]

########################
# Flavor Tagging
########################

mom = (pjetlogf, pjetlogc)
theta = (thetalinf, thetalinc)

gen_cond_btag = dict()
gen_cond_btag[bottom] = "p.Flavor == 5"
gen_cond_btag[charm] = "p.Flavor == 4"
gen_cond_btag[up] = "p.Flavor < 4"

gen_cond_ctag = gen_cond_btag

gen_cond_tautag = dict()
gen_cond_tautag[tau] = "p.TauFlavor == 15"
gen_cond_tautag[up] = "p.TauFlavor == 0"


eff_btag_plots = EfficiencyTaggingBlock("Jet", "b-tag", mom, theta, gen_cond_btag, [0], "BTag", eff_tag_plots)

eff_ctag_plots = EfficiencyTaggingBlock("Jet", "c-tag", mom, theta, gen_cond_ctag, [1], "BTag", eff_tag_plots)

eff_tautag_plots = EfficiencyTaggingBlock("Jet", "tau-tag", mom, theta, gen_cond_tautag, [0], "TauTag", eff_tag_plots)

eff_tag_dict = dict()
eff_tag_dict["b-tag"] = eff_btag_plots
eff_tag_dict["c-tag"] = eff_ctag_plots
eff_tag_dict["tau-tag"] = eff_tautag_plots

algos = ["b-tag", "c-tag", "tau-tag"]
efftype = dict()
efftype["b-tag"] = ["b-tag", "c-mistag", "light-mistag"]
efftype["c-tag"] = ["c-tag", "b-mistag", "light-mistag"]
efftype["tau-tag"] = ["tau-tag", "light-mistag"]

type_part = dict()
type_part["b-tag"] = bottom
type_part["c-mistag"] = charm
type_part["light-mistag"] = up
type_part["c-tag"] = charm
type_part["b-mistag"] = bottom
type_part["tau-tag"] = tau

wp_dict = dict()
wp_dict[("b-tag")] = 0
wp_dict[("c-tag")] = 1
wp_dict[("tau-tag")] = 0

for algo in algos:
    report[algo] = OrderedDict()
    report[algo]["efficiency"] = OrderedDict()
    for type in efftype[algo]:
        title = "{}".format(type)
        # print(algo, title, type_part[type].name, wp, wp_dict[(algo, wp)])
        report[algo]["efficiency"][title] = [
            eff_tag_dict[algo][(type_part[type], wp_dict[(algo)], "pt")],
            eff_tag_dict[algo][(type_part[type], wp_dict[(algo)], "eta")],
        ]
