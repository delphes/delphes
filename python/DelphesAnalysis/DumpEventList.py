#!/usr/bin/env python 
import os
import ROOT
import Delphes
from AnalysisEvent import AnalysisEvent
import EventSelection

def DumpEventList(category, path="./", output="eventlist.txt"):
  """Dump a list of events in a given category"""
  # input
  if os.path.isdir(path):
    dirList=os.listdir(path)
    files=[]
    for fname in dirList:
      files.append(path+fname)
  elif os.path.isfile(path):
    files=[path]
  else:
    files=[]
  # events
  events = AnalysisEvent(files)
  # output
  event_list = open(output,"w")
  # collections and producers used in the analysis
  EventSelection.prepareAnalysisEvent(events)
  for event in events:
    # check category
    if EventSelection.isInCategory(category, event.category):
      # print
      print >> event_list , "Event", event.event()

if __name__=="__main__":
  import sys
  DumpEventList(int(sys.argv[1]), path=sys.argv[2], output="eventlist.txt")

