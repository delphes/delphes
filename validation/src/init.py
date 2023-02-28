import ROOT, os, sys
import importlib

ROOT.gSystem.Load("libDelphes")
try:
    ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
    ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
    pass

inputFile = sys.argv[1]
outputFile = sys.argv[2]
config = sys.argv[3]

config_module = os.path.basename(config).strip(".py")
cfg = importlib.import_module(config_module)

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# get magnetic field
Bz = treeReader.GetInfo("Bz")
if Bz < -9:
    Bz = 4

branch = dict()
collection_names = [
    "Particle",
    "Track",
    "Tower",
    "EFlowTrack",
    "EFlowPhoton",
    "EFlowNeutralHadron",
    "Electron",
    "Muon",
    "Photon",
    "PFJet",
    "Jet",
    "CaloJet",
    "GenJet",
]

for colname in collection_names:
    branch[colname] = treeReader.UseBranch(colname)
