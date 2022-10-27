# from src.utils import *
from src.init import *

debug = False
# debug = False

# reads from input file which gun we are running
gun_pid = int(inputFile.split("_")[1])

plots = cfg.eff_plots + cfg.eff_tag_plots + cfg.reso_plots

for plot in plots:
    plot.construct(gun_pid)


# Loop over all events
if debug:
    numberOfEntries = 100

for entry in range(0, numberOfEntries):

    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)

    if (entry + 1) % 100 == 0:
        print(" ... processed {} events ...".format(entry + 1))

    stable_particles = [p for p in branch["Particle"] if p.Status == 1]

    for eff_plot in cfg.eff_plots:
        for eff_hist in eff_plot.eff_histos:
            if eff_hist.particle.pid == gun_pid:
                # print("eff part hist {} {}".format(eff_hist.particle.name, inputFile))
                eff_hist.fill(stable_particles, branch[eff_hist.collection])

    for eff_plot in cfg.eff_tag_plots:
        for eff_hist in eff_plot.eff_histos:
            if eff_hist.particle.pid == gun_pid:
                # print("eff tag hist {} {}".format(eff_hist.particle.name, inputFile))
                eff_hist.fill(branch[eff_hist.collection])

    for reso_plot in cfg.reso_plots:
        for reso_hist in reso_plot.res_histos:
            if reso_hist.particle.pid == gun_pid:
                if reso_hist.particle.pid in [1, 2, 3, 4, 5, 21]:
                    # print("reso jet hist {} {}".format(reso_hist.particle.name, inputFile))
                    reso_hist.fill(branch["GenJet"], branch[reso_hist.collection])
                else:
                    # print("reso hist {} {}".format(reso_hist.particle.name, inputFile))
                    reso_hist.fill(stable_particles, branch[reso_hist.collection])

# write tree
out_root = ROOT.TFile(outputFile, "RECREATE")

for plot in plots:
    plot.write_histos(gun_pid)
