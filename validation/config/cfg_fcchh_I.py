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
scenario = "I"
card = "FCC/scenarios/validation_delphes_card_FCChh_{}.tcl".format(scenario)
name = os.path.basename(card.replace(".tcl", ""))

## run parameters:
njobs = 100
nevts_per_job = 25000

pmin = 0.1
pmax = 50000.0
log = True
etamin = -6.0
etamax = 6.0
ecm = 100000.0

## plotting parameters
nbins = 100

nbins_ptlin = nbins
ptlin_min = 0.0
ptlin_min = 0.0
ptlin_max = pmax

nbins_ptlog = nbins
ptlog_min = pmin
ptlog_max = pmax

nbins_elin = nbins
elin_min = 0.0
elin_max = pmax

nbins_elog = nbins
elog_min = pmin
elog_max = pmax

nbins_eta = nbins
eta_min = etamin
eta_max = etamax

bins_etalinf = compute_bins(eta_min, eta_max, nbins_eta, 2, "lin")

bins_ptlinf = compute_bins(ptlin_min, ptlin_max, nbins_ptlin, 2, "lin")
bins_ptlogf = compute_bins(ptlog_min, ptlog_max, nbins_ptlog, 2, "log")

bins_jetptlinf = compute_bins(20.0, ptlin_max, nbins_ptlin, 2, "lin")
bins_jetptlogf = compute_bins(20.0, ptlog_max, nbins_ptlog, 2, "log")

bins_elinf = compute_bins(elin_min, elin_max, nbins_elin, 2, "lin")
bins_elogf = compute_bins(elog_min, elog_max, nbins_elog, 2, "log")

bins_jetelinf = compute_bins(20.0, elin_max, nbins_elin, 2, "lin")
bins_jetelogf = compute_bins(20.0, elog_max, nbins_elog, 2, "log")

bins_etalinc = [(0.0, 2.5), (2.5, 4.0), (4.0, 6.0)]
bins_ptlogc = [
    (0.5, 1.5),
    (5, 15),
    (20, 50),
    (50, 100),
    (100, 500),
    (500, 1500),
]

bins_jetptlogc = [
    (20, 50),
    (50, 100),
    (100, 500),
    (500, 1500),
]

bins_elogc = bins_ptlogc
bins_jetelogc = bins_jetptlogc

# observables (name, varname, label, nbins, xmin, xmax, opt)
# last 4 params only used if obs in ResolutionHisto
obs_pt = Observable("pt", "PT", "p_{T}", "[GeV]", 500, 0.8, 1.2, "rel")
obs_p = Observable("p", "P", "p", "[GeV]", 500, 0.9, 1.1, "rel")
obs_e = Observable("e", "E", "E", "[GeV]", 500, 0.5, 1.5, "rel")
obs_d0 = Observable("d0", "D0", "d_0", "[mm]", 500, -0.2, 0.2, "abs")
obs_dz = Observable("dz", "DZ", "d_z", "[mm]", 500, -0.2, 0.2, "abs")
obs_phi = Observable("phi", "Phi", "\phi", " [rad]", 500, -0.001, 0.001, "abs")
obs_theta = Observable("theta", "Theta", "\\theta", "[rad]", 500, -0.001, 0.001, "abs")

# only used in slicing, last 4 params are dummy
obs_eta = Observable("eta", "Eta", "|\eta|", "", 500, -0.2, 0.2, "abs")

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


ptlogf = Slices1D(obs_pt, bins_ptlogf, "log")
ptjetlogf = Slices1D(obs_pt, bins_jetptlogf, "log")

elogf = Slices1D(obs_e, bins_elogf, "log")
ejetlogf = Slices1D(obs_e, bins_jetelogf, "log")
etalinf = Slices1D(obs_eta, bins_etalinf, "linear")

ptlogc = Slices1D(obs_pt, bins_ptlogc, "log")
ptjetlogc = Slices1D(obs_pt, bins_jetptlogc, "log")
elogc = Slices1D(obs_e, bins_elogc, "log")
etalinc = Slices1D(obs_eta, bins_etalinc, "linear")

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
mom = (ptlogf, ptlogc)
eta = (etalinf, etalinc)
obs = [obs_p]

tracks = [pion, muon]
reso_track_plots = ResolutionBlock(tracks, "Track", obs, mom, eta, reso_plots)

report["tracks"]["resolution"] = OrderedDict()
for p in tracks:
    for o in obs:
        title = "{} track resolution: ${}$".format(p.label, o.label)
        report["tracks"]["resolution"][title] = [
            reso_track_plots[((p, o, "pt"))],
            reso_track_plots[((p, o, "eta"))],
        ]

# efficiency
mom = (ptlogf, ptlogc)
eta = (etalinf, etalinc)
tracks = [pion, muon, electron]
eff_track_plots = EfficiencyParticleBlock(tracks, "Track", mom, eta, eff_plots)

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
eta = (etalinf, etalinc)
obs = [obs_e, obs_eta, obs_phi]

caloparticles = [electron, photon, pion, kl]
reso_calo_plots = ResolutionBlock(caloparticles, "Tower", obs, mom, eta, reso_plots)

report["calorimeter"]["resolution"] = OrderedDict()
for p in caloparticles:
    for o in obs:
        title = "{} calorimeter resolution: ${}$".format(p.label, o.label)
        report["calorimeter"]["resolution"][title] = [
            reso_calo_plots[((p, o, "pt"))],
            reso_calo_plots[((p, o, "eta"))],
        ]

# efficiency
mom = (ptlogf, ptlogc)
eta = (etalinf, etalinc)
caloparticles = [electron, photon, pion, kl]
eff_calo_plots = EfficiencyParticleBlock(caloparticles, "Tower", mom, eta, eff_plots)

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
eta = (etalinf, etalinc)
obs = [obs_e, obs_eta, obs_phi]

reso_pftrack_plots = ResolutionBlock(pftracks, "EFlowTrack", obs, mom, eta, reso_plots)
reso_pfgamma_plots = ResolutionBlock(pfgamma, "EFlowPhoton", obs, mom, eta, reso_plots)
reso_pfnh_plots = ResolutionBlock(pfnh, "EFlowNeutralHadron", obs, mom, eta, reso_plots)


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
    {obs_eta: (0.0, 2.5)},
    "track, 0 < $|\eta|$ < 2.5",
)

reso_elecalo_e = ResolutionHisto(
    electron,
    "Tower",
    obs_e,
    elogf,
    {obs_eta: (0.0, 2.5)},
    "calo, 0 < $|\eta|$ < 2.5",
)

reso_elepf_e = ResolutionHisto(
    electron,
    "EFlowTrack",
    obs_e,
    elogf,
    {obs_eta: (0.0, 2.5)},
    "particle-flow, 0 < $|\eta|$ < 2.5",
)

reso_ele_pftest_histos = [reso_eletrk_e, reso_elecalo_e, reso_elepf_e]
reso_ele_pftest_plot = ResolutionPlot("reso_ele_pftest_e", reso_ele_pftest_histos, Text("electrons", (0.5, 0.5)))

reso_pitrk_e = ResolutionHisto(
    pion,
    "Track",
    obs_e,
    elogf,
    {obs_eta: (0.0, 2.5)},
    "track, 0 < $|\eta|$ < 2.5",
)

reso_picalo_e = ResolutionHisto(
    pion,
    "Tower",
    obs_e,
    elogf,
    {obs_eta: (0.0, 2.5)},
    "calo, 0 < $|\eta|$ < 2.5",
)

reso_pipf_e = ResolutionHisto(
    pion,
    "EFlowTrack",
    obs_e,
    elogf,
    {obs_eta: (0.0, 2.5)},
    "particle-flow, 0 < $|\eta|$ < 2.5",
)

reso_pi_pftest_histos = [reso_pitrk_e, reso_picalo_e, reso_pipf_e]
reso_pi_pftest_plot = ResolutionPlot("reso_pi_pftest_e", reso_pi_pftest_histos, Text("charged pions", (0.5, 0.5)))


reso_jetcalo_e = ResolutionHisto(
    up,
    "CaloJet",
    obs_e,
    ejetlogf,
    {obs_eta: (0.0, 2.5)},
    "calo, 0 < $|\eta|$ < 2.5",
)

reso_jetpf_e = ResolutionHisto(
    up,
    "PFJet",
    obs_e,
    ejetlogf,
    {obs_eta: (0.0, 2.5)},
    "particle-flow, 0 < $|\eta|$ < 2.5",
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
eta = (etalinf, etalinc)

eff_pftrack_plots = EfficiencyParticleBlock(pftracks, "EFlowTrack", mom, eta, eff_plots)
eff_pfgamma_plots = EfficiencyParticleBlock(pfgamma, "EFlowPhoton", mom, eta, eff_plots)
eff_pfnh_plots = EfficiencyParticleBlock(pfnh, "EFlowNeutralHadron", mom, eta, eff_plots)

report["particle-flow"]["efficiency"] = OrderedDict()
for p in pftracks:
    title = "{} particle-flow efficiency".format(p.label)
    report["particle-flow"]["efficiency"][title] = [
        eff_pftrack_plots[((p, "pt"))],
        eff_pftrack_plots[((p, "eta"))],
    ]
for p in pfgamma:
    title = "{} particle-flow efficiency".format(p.label)
    report["particle-flow"]["efficiency"][title] = [
        eff_pfgamma_plots[((p, "pt"))],
        eff_pfgamma_plots[((p, "eta"))],
    ]
for p in pfnh:
    title = "{} particle-flow efficiency".format(p.label)
    report["particle-flow"]["efficiency"][title] = [
        eff_pfnh_plots[((p, "pt"))],
        eff_pfnh_plots[((p, "eta"))],
    ]

########################
# Electron
########################

report["electron"] = OrderedDict()

# efficiency
mom = (ptlogf, ptlogc)
eta = (etalinf, etalinc)
coll = [electron]
eff_electron_plots = EfficiencyParticleBlock(coll, "Electron", mom, eta, eff_plots)

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
mom = (ptlogf, ptlogc)
eta = (etalinf, etalinc)
coll = [muon]
eff_muon_plots = EfficiencyParticleBlock(coll, "Muon", mom, eta, eff_plots)

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
mom = (ptlogf, ptlogc)
eta = (etalinf, etalinc)
coll = [photon]
eff_photon_plots = EfficiencyParticleBlock(coll, "Photon", mom, eta, eff_plots)

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

mom = (ptjetlogf, ptjetlogc)
eta = (etalinf, etalinc)

gen_cond_btag = dict()
gen_cond_btag[bottom] = "p.Flavor == 5"
gen_cond_btag[charm] = "p.Flavor == 4"
gen_cond_btag[up] = "p.Flavor < 4"

gen_cond_ctag = gen_cond_btag

gen_cond_tautag = dict()
gen_cond_tautag[tau] = "p.TauFlavor == 15"
gen_cond_tautag[up] = "p.TauFlavor == 0"


eff_btag_plots = EfficiencyTaggingBlock("Jet", "b-tag", mom, eta, gen_cond_btag, [0, 1, 2], "BTag", eff_tag_plots)

eff_ctag_plots = EfficiencyTaggingBlock("Jet", "c-tag", mom, eta, gen_cond_ctag, [4, 5, 6], "BTag", eff_tag_plots)

eff_tautag_plots = EfficiencyTaggingBlock(
    "Jet", "tau-tag", mom, eta, gen_cond_tautag, [0, 1, 2], "TauTag", eff_tag_plots
)

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

wps = ["Loose", "Medium", "Tight"]

wp_dict = dict()
wp_dict[("b-tag", "Loose")] = 0
wp_dict[("b-tag", "Medium")] = 1
wp_dict[("b-tag", "Tight")] = 2
wp_dict[("c-tag", "Loose")] = 4
wp_dict[("c-tag", "Medium")] = 5
wp_dict[("c-tag", "Tight")] = 6
wp_dict[("tau-tag", "Loose")] = 0
wp_dict[("tau-tag", "Medium")] = 1
wp_dict[("tau-tag", "Tight")] = 2

for algo in algos:
    report[algo] = OrderedDict()
    report[algo]["efficiency"] = OrderedDict()
    for type in efftype[algo]:
        for wp in wps:
            title = "{} - {}".format(type, wp)
            # print(algo, title, type_part[type].name, wp, wp_dict[(algo, wp)])
            report[algo]["efficiency"][title] = [
                eff_tag_dict[algo][(type_part[type], wp_dict[(algo, wp)], "pt")],
                eff_tag_dict[algo][(type_part[type], wp_dict[(algo, wp)], "eta")],
            ]
