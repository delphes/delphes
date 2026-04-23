import DelphesPython as delphes
import matplotlib.pyplot as plt


processor = delphes.Delphes()
processor.reader = delphes.Reader("LHEF",
    inputFiles = ["/home/laurent/work/dev/cepgen/build/test.lhe"],
)
processor.modules = dict(
    ParticlePropagator = delphes.Module("ParticlePropagator",
        InputArray = "Delphes/stableParticles",
        OutputArray = "stableParticles",
        ChargedHadronOutputArray = "chargedHadrons",
        ElectronOutputArray = "electrons",
        MuonOutputArray = "muons",
        Radius = 1.14,     # radius of the magnetic field coverage, in m
        HalfLength = 2.9,  # half-length of the magnetic field coverage, in m
        Bz = 3.8,          # magnetic field
    ),
    ChargedHadronTrackingEfficiency = delphes.Module("Efficiency",
        InputArray = "ParticlePropagator/chargedHadrons",
        OutputArray = "chargedHadrons",
        # add EfficiencyFormula {efficiency formula as a function of eta and pt}
        # tracking efficiency formula for charged hadrons
        # set to 1 for full range (-0.2 eta already subtracted from edges)
        EfficiencyFormula = "(eta <= 4.9 && eta >= -4.3) * (1.0) + \
                             (eta  > 4.9 || eta <  -4.3) * (0.0)",
    ),
)

num_events = 10000

muon_pts = []
muon_etas = []
i = 0
while i < num_events:
    event = processor.next()
    muon_pts += [muon.pt for muon in event["ParticlePropagator/muons"]]
    muon_etas += [muon.eta for muon in event["ParticlePropagator/muons"]]
    i += 1

fig, ax = plt.subplots(figsize=(9,5))
ax.hist([muon_pts], bins=15, histtype='step', linewidth=2, alpha=0.7)
plt.xlabel('Muon $p_T$ (GeV)')
plt.ylabel('Frequency')
plt.show()

fig, ax = plt.subplots(figsize=(9,5))
ax.hist([muon_etas], bins=15, histtype='step', linewidth=2, alpha=0.7)
plt.xlabel('Muon $\\eta$')
plt.ylabel('Frequency')
plt.show()

