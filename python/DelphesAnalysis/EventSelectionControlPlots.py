import ROOT
import sys
import os
from AnalysisEvent import AnalysisEvent
from BaseControlPlots import BaseControlPlots
from CPconfig import configuration
import EventSelection

# Requirements:
# fully implemented EventSelection
# categories stored as event.category

class EventSelectionControlPlots(BaseControlPlots):
    """A class to create control plots for event selection"""

    def __init__(self, dir=None, dataset=None, mode="plots"):
      # create output file if needed. If no file is given, it means it is delegated
      BaseControlPlots.__init__(self, dir=dir, purpose="eventSelection", dataset=dataset, mode=mode)
      self.eventCategories = EventSelection.eventCategories()
 
    def beginJob(self):
      # declare histograms
      self.add("event","Event number",1000,0,100000)
      self.add("category","event category",self.eventCategories+1,0,self.eventCategories+1)

    def process(self, event):
      """EventSelectionControlPlots"""
      result = { }
      ## event category
      categoryData = event.category
      result["category"] = [ ]
      for category in range(self.eventCategories):
        if EventSelection.isInCategory(category, categoryData):
          result["category"].append(category)
      result["event"] = event.event()
      return result

if __name__=="__main__":
  import sys
  from BaseControlPlots import runTest
  runTest(sys.argv[1], EventSelectionControlPlots())

