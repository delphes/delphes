from src.utils import *
from src.init import *

debug = False
# debug = False

for plot in cfg.eff_plots:
    plot.construct()
    # print(eff_plot_name)
    for eff_hist in plot.eff_histos:
        print(eff_hist.histname)

# Loop over all events
if debug:
    numberOfEntries = 100

for entry in range(0, numberOfEntries):
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)

    stable_particles = [p for p in branch["Particle"] if p.Status == 1]

    for eff_plot in cfg.eff_plots:
        for eff_hist in eff_plot.eff_histos:
            eff_hist.fill(stable_particles, branch[eff_hist.collection])


# write tree
out_root = ROOT.TFile(outputFile, "RECREATE")

for plot in cfg.eff_plots:
    plot.write_histos()
