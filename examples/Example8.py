#!/usr/bin/env python
import sys
import ROOT, math
from collections import OrderedDict
from ROOT import TLorentzVector
import array


# _______________________________________________________________________________
def print_part(p, all_particles, basespace):
    print(
        "{}i: {}, PID: {}, Q: {}, Status: {}, E: {:.2e}, Eta: {:.2e}, Phi: {:.2e}, X: {:.2e}, Y: {:.2e}, Z: {:.2e}, D1: {}, D2: {}".format(
            basespace,
            all_particles.index(p),
            p.PID,
            p.Charge,
            p.Status,
            p.E,
            p.Eta,
            p.Phi,
            p.X,
            p.Y,
            p.Z,
            p.D1,
            p.D2,
        )
    )


# _______________________________________________________________________________
def print_descendence(particle, all_particles, basespace, list_indexed):

    print_part(particle, all_particles, basespace)
    d1 = particle.D1
    d2 = particle.D2

    position = all_particles.index(particle)

    # exit recursion if particle already treated
    if position in list_indexed:
        # print("{} particle descendence already studied".format(position))
        return

    list_indexed.append(position)

    basespace += "-- "

    drange = []
    maxd = max(d1, d2)
    mind = min(d1, d2)

    if maxd < 0:
        # stable particle
        drange = []
    elif mind < 0:
        drange = [maxd]
    elif d1 > d2:
        drange = [d1, d2]
    else:
        drange = range(mind, maxd + 1)

    # exit recursion if particle is stable
    if len(drange) == 0:
        # print("{} reached stable particle".format(basespace))
        return

    # this particle was not treated and is unstable
    else:
        for i in drange:
            part = all_particles[i]
            print_descendence(part, all_particles, basespace, list_indexed)


# _______________________________________________________________________________
if len(sys.argv) < 2:
    print(" Usage: Example1.py input_file")
    sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
    ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
    ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
    pass

inputFile = sys.argv[1]
outputFile = sys.argv[2]

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointers to branches used in this analysis

branchParticle = treeReader.UseBranch("Particle")

debug = True
if debug:
    numberOfEntries = 10


# Loop over all events
for entry in range(0, numberOfEntries):
    # for entry in range(0, 1000):
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)

    if (entry + 1) % 100 == 0:
        print(" ... processed {} events ...".format(entry + 1))
    i = -1

    if debug:
        print(
            "---------------------------------------------------------------------------------------"
        )
        print("")

    ###  first produce 3 lists: higgs daughters stable, z daughters stables, all stable particles
    list_indexed = []
    for gen in branchParticle:
        i += 1

        """
        if debug:
            print(
                "N: ",
                i,
                ", St: ",
                gen.Status,
                ", PID: ",
                gen.PID,
                ", E: ",
                gen.E,
                ", PT: ",
                gen.PT,
                ", Eta: ",
                gen.Eta,
                ", M: ",
                gen.Mass,
                ", M1: ",
                gen.M1,
                ", M2: ",
                gen.M2,
                ", D1: ",
                gen.D1,
                ", D2: ",
                gen.D2,
            )
        """
        basespace = ""
        if i not in list_indexed:
            print_descendence(gen, branchParticle, basespace, list_indexed)

# Show resulting histograms
out_root = ROOT.TFile(outputFile, "RECREATE")
