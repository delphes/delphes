from BaseControlPlots import BaseControlPlots

# Requirements:
# event.muons
# event.electrons

class LeptonControlPlots(BaseControlPlots):
    """A class to create control plots for leptons"""

    def __init__(self, dir=None, dataset=None, mode="plots"):
      # create output file if needed. If no file is given, it means it is delegated
      BaseControlPlots.__init__(self, dir=dir, purpose="leptons", dataset=dataset, mode=mode)

    def beginJob(self):
      # declare histograms
      self.add("ElectronPt","Electron Pt",100,0,200)
      self.add("MuonPt","Muon Pt",100,0,200)
      self.add("ElectronEta","Electron Eta",50,-2.5,2.5)
      self.add("MuonEta","Muon Eta",50,-2.5,2.5)
      self.add("NMuons","Muon multiplicity",10,0,10)
      self.add("NElectrons","Electron multiplicity",10,0,10)

    def process(self, event):
      #get information
      result = { }
      result["ElectronPt"] = [ ]
      result["MuonPt"] = [ ]
      result["ElectronEta"] = [ ]
      result["MuonEta"] = [ ]
      for mu in event.muons:
        result["MuonPt"].append(mu.PT)
        result["MuonEta"].append(mu.Eta)
      for ele in event.electrons:
        result["ElectronPt"].append(ele.PT)
        result["ElectronEta"].append(ele.Eta)
      result["NMuons"] = event.muons.GetEntries()
      result["NElectrons"] = event.electrons.GetEntries()
      return result

if __name__=="__main__":
  import sys
  from DelphesAnalysis.BaseControlPlots import runTest
  runTest(sys.argv[1], LeptonControlPlots())

