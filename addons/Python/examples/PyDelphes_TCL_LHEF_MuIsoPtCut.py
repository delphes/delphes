import DelphesPython as delphes
import matplotlib.pyplot as plt


#delphes.load("libDelphesPythia.so")

processor = delphes.Delphes()
processor.reader = delphes.Reader("LHEF", inputFiles = ["test.lhe"])
'''processor.reader = delphes.Reader("Pythia8",
    pythiaConfig = [
        'Beams:idA = 2212               ! proton beam 1',
        'Beams:idB = 2212               ! proton beam 2',
        'Beams:eCM = 13600.             ! 13.6 TeV center-of-mass energy',
        'HardQCD:all = on               ! Enable all 2-to-2 QCD processes (qg, gg, qq, etc.)',
        'PhaseSpace:pThatMin = 50.      ! Minimum pT of hard scattering (GeV)',
    ],
)'''
processor.loadTCL("../cards/delphes_card_CMS.tcl")
print(f"Bz={processor.modules['ParticlePropagator']['Bz']}T")

num_events = 15000

ptcuts = [5., 10., 20., 50., 100.]

muon_pts_vs_ptcuts = []
muon_etas_vs_ptcuts = []

for ptmin in ptcuts:
    #processor.reset()  # ensure to always use the same events
    processor.modules['MuonIsolation']['PTMin'] = ptmin
    print(f"Min muon pt={processor.modules['MuonIsolation']['PTMin']} GeV/c")

    muon_pts = []
    muon_etas = []
    i = 0
    while i < num_events:
        event = processor.next()
        muon_pts += [muon.pt for muon in event["UniqueObjectFinder/muons"]]
        muon_etas += [muon.eta for muon in event["UniqueObjectFinder/muons"]]
        i += 1
    muon_pts_vs_ptcuts.append(muon_pts)
    muon_etas_vs_ptcuts.append(muon_etas)

labels = [f"$p_T^{{\\mu}}$ > {ptmin} GeV/c" for ptmin in ptcuts]

fig, ax = plt.subplots(figsize=(9,5))
ax.hist(muon_pts_vs_ptcuts, bins=100, histtype='step', linewidth=2, alpha=0.7, label=labels, log=True)
ax.legend()
plt.xlabel('$p_T^\\mu$')
plt.ylabel('Frequency')
plt.show()

fig, ax = plt.subplots(figsize=(9,5))
ax.hist(muon_etas_vs_ptcuts, bins=50, histtype='step', linewidth=2, alpha=0.7, label=labels)
ax.legend()
plt.xlabel('$\\eta^\\mu$')
plt.ylabel('Frequency')
plt.show()
