import os
from src.utils import compute_bins, Observable, Particle, Slices1D, EfficiencyBlock


"""
1) prepare validation card:

- validation_delphes_card_FCChh_I.tcl

2) make sure all the other needed cards are in the same directory:

- muonMomentumResolution_I.tcl
- trackMomentumResolution_I.tcl
- FCChh_I.tcl

"""

### global parameters for event generation
delphes_path = "/Users/michele/Programs/delphes/"

## provide relative path to the directory
scenario = "I"
card = "FCC/scenarios/validation_delphes_card_FCChh_{}.tcl".format(scenario)
name = os.path.basename(card.replace(".tcl", ""))

## run parameters:
njobs = 2
nevts_per_job = 1000

pmin = 0.1
pmax = 50000.0
log = True
etamin = -6.0
etamax = 6.0
ecm = 100000.0

## plotting parameters
nbins = 50

nbins_ptlin = nbins
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

bins_etalinf = compute_bins(eta_min, eta_max, nbins_eta, 1, "lin")

bins_ptlinf = compute_bins(ptlin_min, ptlin_max, nbins_ptlin, 1, "lin")
bins_ptlogf = compute_bins(ptlog_min, ptlog_max, nbins_ptlog, 2, "log")

bins_elinf = compute_bins(elin_min, elin_max, nbins_elin, 1, "lin")
bins_elogf = compute_bins(elog_min, elog_max, nbins_elog, 2, "log")

bins_etalinc = [(0.0, 2.5), (2.5, 4.0), (4.0, 6.0)]
bins_ptlogc = [
    (0.5, 1.5),
    (5, 15),
    (50, 150),
    (500, 1500),
]
bins_elogc = bins_ptlogc

# observables (name, varname, label, nbins, xmin, xmax, opt)
# last 4 params only used if obs in ResolutionHisto
obs_pt = Observable("pt", "PT", "p_{T}", "[GeV]", 500, 0.8, 1.2, "rel")
obs_p = Observable("p", "P", "p", "[GeV]", 500, 0.9, 1.1, "rel")
obs_e = Observable("e", "E", "E", "[GeV]", 500, 0.5, 1.5, "rel")
obs_d0 = Observable("d0", "D0", "d_0'", "[mm]", 500, -0.2, 0.2, "abs")
obs_dz = Observable("dz", "DZ", "d_z", "[mm]", 500, -0.2, 0.2, "abs")
obs_phi = Observable("phi", "Phi", "\phi", " [rad]", 500, -0.001, 0.001, "abs")
obs_theta = Observable("theta", "Theta", "\theta", "[rad]", 500, -0.001, 0.001, "abs")

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

# jet particles flavour
bjetl = Particle("bjetl", "b-jet (L)", 5)
bjetm = Particle("bjetm", "b-jet (M)", 5)
bjett = Particle("bjett", "b-jet (T)", 5)

taujetl = Particle("taujetl", "tau-jet (L)", 15)
taujetm = Particle("taujetm", "tau-jet (M)", 15)
taujett = Particle("taujett", "tau-jet (T)", 15)


ptlogf = Slices1D(obs_pt, bins_ptlogf, "log")
elogf = Slices1D(obs_e, bins_elogf, "log")
etalinf = Slices1D(obs_eta, bins_etalinf, "linear")

ptlogc = Slices1D(obs_pt, bins_ptlogc, "log")
elogc = Slices1D(obs_e, bins_elogc, "log")
etalinc = Slices1D(obs_eta, bins_etalinc, "linear")


### definition of efficiency plots
eff_plots = []  ## all eff plots (for step1)

mom = (ptlogf, ptlogc)
eta = (etalinf, etalinc)

tracks = [pion, muon, electron]
eff_track_plots = EfficiencyBlock(tracks, "Track", mom, eta, eff_plots)

mom = (elogf, elogc)
eta = (etalinf, etalinc)

eftracks = [pion, electron]
eff_eftrack_plots = EfficiencyBlock(eftracks, "EFlowTrack", mom, eta, eff_plots)

efphotons = [photon]
eff_efphoton_plots = EfficiencyBlock(efphotons, "EFlowPhoton", mom, eta, eff_plots)

efnhs = [kl]
eff_efnh_plots = EfficiencyBlock(efnhs, "EFlowNeutralHadron", mom, eta, eff_plots)


report = {
    "Tracking": {
        "Efficiency": {
            "pion": [
                eff_track_plots[0][pion],  ## pt
                eff_track_plots[1][pion],  ## eta
            ],
            "muon": [
                eff_track_plots[0][muon],
                eff_track_plots[1][muon],
            ],
            "electron": [
                eff_track_plots[0][electron],
                eff_track_plots[1][electron],
            ],
        },
    },
    "Particle Flow": {
        "Efficiency": {
            "pion": [eff_eftrack_plots[0][pion], eff_eftrack_plots[1][pion]],
            "electron": [
                eff_eftrack_plots[0][electron],
                eff_eftrack_plots[1][electron],
            ],
            "photon": [
                eff_efphoton_plots[0][photon],
                eff_efphoton_plots[1][photon],
            ],
            "KL": [eff_efnh_plots[0][kl], eff_efnh_plots[1][kl]],
        }
    },
}


"""
Structure of report:

Tracks:
    Efficiencies:
        pion pt, eta
        muon pt, eta
        ele, pt, eta
    Resolution:
        pt:
         pion vs pt, eta
         muon vs pt, eta
         ele vs pt, eta
        p:
         pion vs pt, eta
         muon vs pt, eta
         ele vs pt, eta
        d0:
         pion vs pt, eta
         muon vs pt, eta
         ele vs pt, eta
         ...

Calo:
    Resolution:
       e:
        ele e, eta
        gamma, e, eta
        pion e, eta
        kl e, eta

PFlow:
    Efficiency:
        ele e, eta
        gamma, e, eta
        pion e, eta
        kl e, eta
    Resolution:
        e:
          ele vs e (trk, pftrack, calo) in central bin / pion vs e (trk, pftrack, calo) /jet vs e (calo, pf) (gluon/light)

Jets:

    Resolution:
        e:
          pf jet e eta (g)
          pf jet e, eta (c)
          pf jet e, eta (b)

    Tagging Efficiency:

          btag: pt (3WPs) (eta central) / eta (3WP), pt > 30 - 60 GeV
          btag: pt (3WPs) (eta fwd) / eta (3WP), pt > 60 - 200 GeV

          ctag: pt (3WPs) (eta central) / eta (3WP), pt > 30 - 60 GeV
          ctag: pt (3WPs) (eta fwd) / eta (3WP), pt > 60 - 200 GeV

          tautag: pt (3WPs) (eta central) / eta (3WP), pt > 30 - 60 GeV
          tautag: pt (3WPs) (eta fwd) / eta (3WP), pt > 60 - 200 GeV

electron
   pt, eta
mu
   pt, eta
gamma
   pt, eta


Efficiency plot (section)
"""
