#!/usr/bin/env python
import sys
import ROOT, math
from collections import OrderedDict
from ROOT import TLorentzVector
import array

#_________________________________________________________________________________________


if len(sys.argv) < 2:
  print " Usage: Example1.py input_file"
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

branchParticle   = treeReader.UseBranch("Particle")

## fill histogram with particle status
hStatus            = ROOT.TH1F("hStatus", "hStatus",  101, -0.5,  100.5)

debug = True
if debug: numberOfEntries = 10


# Loop over all events
for entry in range(0, numberOfEntries):
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)

    if (entry+1)%100 == 0:
       print ' ... processed {} events ...'.format(entry+1)
    i=0

    print '---------------------------------------------------------------------------------------'
    print ''


    for gen in branchParticle:
        i +=1
        print "N: ", i, ", St: ", gen.Status,", PID: ",gen.PID,", E: ",gen.E,", PT: ",gen.PT,", Eta: ",gen.Eta,", M: ",gen.Mass,", M1: ",gen.M1,", M2: ",gen.M2,", D1: ",gen.D1,", D2: ",gen.D2        

        hStatus.Fill(gen.Status)

# Show resulting histograms
out_root = ROOT.TFile(outputFile,"RECREATE")
 
hStatus.Write()
