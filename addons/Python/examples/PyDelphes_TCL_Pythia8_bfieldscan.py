import DelphesPython as delphes
import matplotlib.pyplot as plt


delphes.load("libDelphesPythia.so")

processor = delphes.Delphes()
#processor.reader = delphes.Reader("LHEF", inputFiles = ["/home/laurent/work/dev/cepgen/build/output.lhe"])
processor.reader = delphes.Reader("Pythia8",
    pythiaConfig = [
        'Beams:idA = 2212               ! proton beam 1',
        'Beams:idB = 2212               ! proton beam 2',
        'Beams:eCM = 13600.             ! 13.6 TeV center-of-mass energy',
        'HardQCD:all = on               ! Enable all 2-to-2 QCD processes (qg, gg, qq, etc.)',
        'PhaseSpace:pThatMin = 50.      ! Minimum pT of hard scattering (GeV)',
    ],
)
processor.loadTCL("../cards/delphes_card_LHeC.tcl")

num_events = 500

bfields = [0., 1.2, 2.4, 3.8, 5.2]

ele_pts_vs_bfields = []
ele_etas_vs_bfields = []
muon_pts_vs_bfields = []
muon_etas_vs_bfields = []
jet_pts_vs_bfields = []
jet_etas_vs_bfields = []

for bfield in bfields:
    processor.modules['ParticlePropagator']['Bz'] = bfield
    print(f"Bz={processor.modules['ParticlePropagator']['Bz']}T")

    ele_pts = []
    ele_etas = []
    muon_pts = []
    muon_etas = []
    jet_pts = []
    jet_etas = []
    i = 0
    while i < num_events:
        event = processor.next()
        ele_pts += [ele.pt for ele in event["UniqueObjectFinder/electrons"]]
        ele_etas += [ele.eta for ele in event["UniqueObjectFinder/electrons"]]
        muon_pts += [mu.pt for mu in event["UniqueObjectFinder/muons"]]
        muon_etas += [mu.eta for mu in event["UniqueObjectFinder/muons"]]
        jet_pts += [jet.pt for jet in event["UniqueObjectFinder/jets"]]
        jet_etas += [jet.eta for jet in event["UniqueObjectFinder/jets"]]
        i += 1
    ele_pts_vs_bfields.append(ele_pts)
    ele_etas_vs_bfields.append(ele_etas)
    muon_pts_vs_bfields.append(muon_pts)
    muon_etas_vs_bfields.append(muon_etas)
    jet_pts_vs_bfields.append(jet_pts)
    jet_etas_vs_bfields.append(jet_etas)

labels = [f"B = {bfield}T" for bfield in bfields]

fig, ax = plt.subplots(figsize=(9,5))
ax.hist(ele_pts_vs_bfields, bins=25, histtype='step', linewidth=2, alpha=0.7, label=labels)
ax.legend()
plt.xlabel('Electron $p_T$')
plt.ylabel('Frequency')
plt.show()

fig, ax = plt.subplots(figsize=(9,5))
ax.hist(ele_etas_vs_bfields, bins=25, histtype='step', linewidth=2, alpha=0.7, label=labels)
ax.legend()
plt.xlabel('Electron $\\eta$')
plt.ylabel('Frequency')
plt.show()

fig, ax = plt.subplots(figsize=(9,5))
ax.hist(muon_pts_vs_bfields, bins=25, histtype='step', linewidth=2, alpha=0.7, label=labels, log=True)
ax.legend()
plt.xlabel('Muon $p_T$')
plt.ylabel('Frequency')
plt.show()

fig, ax = plt.subplots(figsize=(9,5))
ax.hist(muon_etas_vs_bfields, bins=25, histtype='step', linewidth=2, alpha=0.7, label=labels)
ax.legend()
plt.xlabel('Muon $\\eta$')
plt.ylabel('Frequency')
plt.show()

fig, ax = plt.subplots(figsize=(9,5))
ax.hist(jet_pts_vs_bfields, bins=50, histtype='step', linewidth=2, alpha=0.7, label=labels, log=True)
ax.legend()
plt.xlabel('Jet $p_T$')
plt.ylabel('Frequency')
plt.show()

fig, ax = plt.subplots(figsize=(9,5))
ax.hist(jet_etas_vs_bfields, bins=50, histtype='step', linewidth=2, alpha=0.7, label=labels)
ax.legend()
plt.xlabel('Jet $\\eta$')
plt.ylabel('Frequency')
plt.show()

