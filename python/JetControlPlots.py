from BaseControlPlots import BaseControlPlots

# Requirements:
# event.jets
# event.bjets
# event.MEt

class JetControlPlots(BaseControlPlots):
    """A class to create control plots for jetmet"""

    def __init__(self, dir=None, dataset=None, mode="plots"):
      # create output file if needed. If no file is given, it means it is delegated
      BaseControlPlots.__init__(self, dir=dir, purpose="jetmet", dataset=dataset, mode=mode)

    def beginJob(self):
      # declare histograms
      self.add("JetPt","Jet Pt",100,0,200)
      self.add("JetEta","Jet Eta",100,-5,5)
      self.add("JetPhi","Jet Phi",64,-3.2,3.2)
      self.add("BjetPt","B Jet Pt",100,0,200)
      self.add("BjetEta","B Jet Eta",100,-5,5)
      self.add("BjetPhi","B Jet Phi",64,-3.2,3.2)
      self.add("MET","MET",100,0,200)
      self.add("METphi","MET phi",64,-3.2,3.2)
      self.add("Njets","Jet multiplicity",10,0,10)
      self.add("Nbjets","B Jet multiplicity",10,0,10)

    def process(self, event):
      #get information
      result = { }
      result["JetPt"] = [ ]
      result["JetEta"] = [ ]
      result["JetPhi"] = [ ]
      result["BjetPt"] = [ ]
      result["BjetEta"] = [ ]
      result["BjetPhi"] = [ ]
      nb=0
      for jet in event.selectedJets:
        result["JetPt"].append(jet.PT)
        result["JetEta"].append(jet.Eta)
        result["JetPhi"].append(jet.Phi)
        if jet.BTag:
          nb+=1
          result["BjetPt"].append(jet.PT)
          result["BjetEta"].append(jet.Eta)
          result["BjetPhi"].append(jet.Phi)
      result["Njets"]  = len(event.selectedJets)
      result["Nbjets"] = nb
      result["MET"] = event.MEt[0].MET
      result["METphi"] = event.MEt[0].Phi
      return result

if __name__=="__main__":
  import sys
  from DelphesAnalysis.BaseControlPlots import runTest
  runTest(sys.argv[1], JetControlPlots())

