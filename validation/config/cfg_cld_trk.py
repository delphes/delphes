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
card = "validation_delphes_card_CLD.tcl"
name = os.path.basename(card.replace(".tcl", ""))

## run parameters:
njobs = 500
nevts_per_job = 25000

#
pmin = 1
pmax = 150.0
log = True
etamin = -3.0
etamax = 3.0
ecm = 400.0

# plotting parameters
nbins = 50

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
theta_min = 0.
theta_max = 90

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

bins_thetalinc = [(85, 95), (40, 50), (0, 10)]
bins_plogc = [
    (0.8, 1.2),
    (8, 12),
    (80, 120),
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
obs_theta = Observable("theta", "Theta", "\\theta", "[°]", 500, -0.01, 0.01, "abs")

## covariance matrix
obs_corr_d0phi = Observable(
    "corr_d0phi", "CorrD0Phi", "\\rho(d_0, \phi)", "", 100, -1.0, 1.0, "mean"
)
obs_corr_d0c = Observable(
    "corr_d0c", "CorrD0C", "\\rho(d_0, C)", "", 100, -1.0, 1.0, "mean"
)
obs_corr_d0dz = Observable(
    "corr_d0dz", "CorrD0DZ", "\\rho(d_0, d_z)", "", 100, -1.0, 1.0, "mean"
)
obs_corr_d0ctgtheta = Observable(
    "corr_d0ctgtheta",
    "CorrD0CtgTheta",
    "\\rho(d_0, ctg(\\theta))",
    "",
    100,
    -1.0,
    1.0,
    "mean",
)
obs_corr_phic = Observable(
    "corr_phic", "CorrPhiC", "\\rho(\phi, C)", "", 100, -1.0, 1.0, "mean"
)
obs_corr_phidz = Observable(
    "corr_phidz", "CorrPhiDZ", "\\rho(\phi, d_z)", "", 100, -1.0, 1.0, "mean"
)
obs_corr_phictgtheta = Observable(
    "corr_phictgtheta",
    "CorrPhiCtgTheta",
    "\\rho(\phi, ctg(\\theta))",
    "",
    100,
    -1.0,
    1.0,
    "mean",
)
obs_corr_cdz = Observable(
    "corr_cdz", "CorrCDZ", "\\rho(C, d_z)", "", 100, -1.0, 1.0, "mean"
)
obs_corr_cctgtheta = Observable(
    "corr_cctgtheta",
    "CorrCCtgTheta",
    "\\rho(C, ctg(\\theta))",
    "",
    100,
    -1.0,
    1.0,
    "mean",
)
obs_corr_dzctgtheta = Observable(
    "corr_dzctgtheta",
    "CorrDZCtgTheta",
    "\\rho(d_z, ctg(\\theta))",
    "",
    100,
    -1.0,
    1.0,
    "mean",
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

# tracks = [pion, muon, electron]
tracks = [pion]
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
