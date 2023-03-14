"""
This script computes and prints the signal efficiency when reinterpreting the CMS analysis searching for LLPs that decay in the endcap muon detectors (https://arxiv.org/abs/2107.04838)
The event-level and cluster-level selections follow the exact selections applied in the paper and the recasting instructions provided in the HEPData entry (https://www.hepdata.net/record/104408)
The user would need to normalize to the correct signal cross section and luminosity to get the expected signal yield.
Need input ROOT file from running Delphes
"""

import sys
import ROOT
import math
def deltaR(eta1, phi1, eta2, phi2):
    return (dPhi(phi1,phi2)**2+(eta1-eta2)**2)**0.5

def dPhi(phi1, phi2):
    delta = phi1-phi2
    while delta > math.pi: delta -= 2* math.pi
    while delta < -math.pi: delta += 2* math.pi
    return delta

if __name__ == '__main__':
    try:
      input = raw_input
    except:
      pass

    if len(sys.argv) < 2:
      print(" Usage: ExampleCscCluster.py input_root_file")
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
    branchCluster      = treeReader.UseBranch("CscCluster")
    branchMET      = treeReader.UseBranch("MissingET")
    branchElectron = treeReader.UseBranch("Electron")
    branchMuon = treeReader.UseBranch("Muon")
    branchJet = treeReader.UseBranch("Jet")
    branchWeight   = treeReader.UseBranch("Weight")
    branchEvent    = treeReader.UseBranch("Event")


    signal_yield = 0
    total_weight = 0

    # Loop over all events
    for entry in range(0, numberOfEntries):
      # Load selected branches with data from specified event
      treeReader.ReadEntry(entry)

      ## main MC event weight
      w =  branchWeight[0].Weight
      total_weight += w

      ################################
      # Event-level selections
      ################################
      # Require MET > 200 GeV
      if branchMET.At(0).MET < 200: continue

      # Require at least 1 jet with pT > 50 GeV and abs(eta) < 2.4
      nJet = 0
      for i in range(branchJet.GetEntries()):
          jet = branchJet.At(i)
          if jet.PT > 50 and abs(jet.Eta)< 2.4: nJet+=1
      if nJet == 0: continue

      # Require 0 lepton
      nLeptons = 0
      for i in range(branchElectron.GetEntries()):
          if branchElectron.At(i).PT > 35 and abs(branchElectron.At(i).Eta)< 2.5: nLeptons+=1
      for i in range(branchElectron.GetEntries()):
          if branchElectron.At(i).PT > 25 and abs(branchElectron.At(i).Eta)< 2.4: nLeptons+=1
      if nLeptons > 0:continue
      ################################
      # Cluster-level selections
      ################################
      nCscCluster = 0
      for i in range(branchCluster.GetEntries()):
        cluster = branchCluster.At(i)

        # check for jet veto
        maxJetVetoPt = 0
        for j in range(branchJet.GetEntries()):
            jet = branchJet.At(j)
            if deltaR(cluster.Eta, cluster.Phi, jet.Eta, jet.Phi) < 0.4:
                maxJetVetoPt = max(maxJetVetoPt, jet.PT)
        nCscCluster+=1
        if (maxJetVetoPt<10 and abs(dPhi(cluster.Phi, branchMET.At(0).Phi)) < 0.75 and cluster.T < 12.5 and cluster.T > -5): nCscCluster+=1
      if nCscCluster == 0:continue

      signal_yield+= w

    print("final signal efficiency is:" + str(signal_yield/total_weight))
