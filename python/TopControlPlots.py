from BaseControlPlots import BaseControlPlots
import ROOT

# Requirements:
#   event.topCandidates_L
#   event.topCandidates_H
#   event.topPairCandidates

class TopControlPlots(BaseControlPlots):
    """A class to create control plots for jetmet"""

    def __init__(self, dir=None, dataset=None, mode="plots"):
      # create output file if needed. If no file is given, it means it is delegated
      BaseControlPlots.__init__(self, dir=dir, purpose="top", dataset=dataset, mode=mode)

    def beginJob(self):
      # declare histograms
      self.add("Hmass","hadronic top mass",1000,0,1000)
      self.add("Hpt","hadronic top transverse momentum",1000,0,1000)
      self.add("bestHmass","hadronic top mass",1000,0,1000)
      self.add("bestHpt","Hadronic top transverse momentum",1000,0,1000)
      self.add("bestLmass","leptonic top transverse mass",1000,0,1000)
      self.add("bestLpt","leptonic top transverse momentum",1000,0,1000)
      self.add("mass1","hadronic top mass",1000,0,1000)
      self.add("pt1","hadronic top transverse momentum",1000,0,1000)
      self.add("mass2","leptonic top transverse mass",1000,0,1000)
      self.add("pt2","leptonic top transverse momentum",1000,0,1000)

    def process(self, event):
      #get information
      result = { }
      result["Hmass"] = [ ]
      result["Hpt"] = [ ]
      for top in event.topCandidates_H:
        result["Hmass"].append((top[0].P4()+top[1].P4()+top[2].P4()).M())
        result["Hpt"].append((top[0].P4()+top[1].P4()+top[2].P4()).Pt())
      if len(event.topCandidates_H):
        top = event.topCandidates_H[0]
        result["bestHmass"] = (top[0].P4()+top[1].P4()+top[2].P4()).M()
        result["bestHpt"] = (top[0].P4()+top[1].P4()+top[2].P4()).Pt()
      if len(event.topCandidates_L):
        top = event.topCandidates_L[0]
        met = ROOT.TLorentzVector()
        met.SetPtEtaPhiE(top[2].MET,0,top[2].Phi,top[2].MET)
        momentum = top[0].P4()+top[1].P4()+met
        momentum.SetPz(0)
        result["bestLmass"] = momentum.M()
        result["bestLpt"] = momentum.Pt()
      if len(event.topPairCandidates):
        top = event.topPairCandidates[0][:3]
        result["mass1"] = (top[0].P4()+top[1].P4()+top[2].P4()).M()
        result["pt1"] = (top[0].P4()+top[1].P4()+top[2].P4()).Pt()
        top = event.topPairCandidates[0][3:]
        met = ROOT.TLorentzVector()
        met.SetPtEtaPhiE(top[2].MET,0,top[2].Phi,top[2].MET)
        momentum = top[0].P4()+top[1].P4()+met
        momentum.SetPz(0)
        result["mass2"] = momentum.M()
        result["pt2"] =momentum.Pt()
      return result

if __name__=="__main__":
  import sys
  from DelphesAnalysis.BaseControlPlots import runTest
  runTest(sys.argv[1], TopControlPlots())

