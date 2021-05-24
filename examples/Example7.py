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
branchJet = treeReader.UseBranch("JetPUPPITight")
branchElectron = treeReader.UseBranch("ElectronMedium")
branchWeight = treeReader.UseBranch("Weight")

# Book histograms
histJetPT = ROOT.TH1F("jet_pt", "jet P_{T}", 100, 0.0, 1000.0)
histElectronPT = ROOT.TH1F("Electron_pt", "electron P_{T}", 100, 0.0, 1000.0)

# Loop over all events
for entry in range(0, numberOfEntries):
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)

  w = 1.0  
  ## read MC event weight
  if branchWeight.GetEntries() > 0:  
    weight = branchWeight.At(0).Weight

  # If event contains at least 1 jet
  if branchJet.GetEntries() > 0:
    # Take first jet
    jet = branchJet.At(0)

    ## 0 - Loose , 1 - Medium, 2 - Tight
    wp = 1

    BtagOk = ( jet.BTag & (1 << wp) )
    pt = jet.PT
    eta = abs(jet.Eta)

    # Plot jet transverse momentum
    if (BtagOk and pt > 30. and eta < 5.):
        histJetPT.Fill(jet.PT, w)


  # If event contains at least 1 electron
  if branchElectron.GetEntries() > 0:
    # Take first electron
    electron = branchElectron.At(0)

    pt = electron.PT
    eta = abs(electron.Eta)

    IsoCut = 0.2
    IsoOk = electron.IsolationVar < IsoCut

    # Plot electron transverse momentum
    if (IsoOk and pt > 10. and eta < 5.):
        histElectronPT.Fill(electron.PT, w)


# Show resulting histograms
cnv = ROOT.TCanvas("cnv", "cnv", 50, 50, 800, 500)
cnv.Divide(2, 1)
cnv.cd(1)
ROOT.gStyle.SetOptStat(0)

histJetPT.Draw()

cnv.cd(2)
ROOT.gStyle.SetOptStat(0)
histElectronPT.Draw()

input("Press Enter to continue...")
