#!/usr/bin/env python

import sys

import ROOT

try:
  input = raw_input
except:
  pass

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

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointers to branches used in this analysis
branchJet = treeReader.UseBranch("Jet")
branchElectron = treeReader.UseBranch("Electron")

# Book histograms
histJetPT = ROOT.TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0)
histMass = ROOT.TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0)

# Loop over all events
for entry in range(0, numberOfEntries):
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)

  # If event contains at least 1 jet
  if branchJet.GetEntries() > 0:
    # Take first jet
    jet = branchJet.At(0)

    # Plot jet transverse momentum
    histJetPT.Fill(jet.PT)

    # Print jet transverse momentum
    print(jet.PT)

  # If event contains at least 2 electrons
  if branchElectron.GetEntries() > 1:
    # Take first two electrons
    elec1 = branchElectron.At(0)
    elec2 = branchElectron.At(1)

    # Plot their invariant mass
    histMass.Fill(((elec1.P4()) + (elec2.P4())).M())

# Show resulting histograms
histJetPT.Draw()
histMass.Draw()

input("Press Enter to continue...")
