#!/usr/bin/env python 
import os
import ROOT
import Delphes
from AnalysisEvent import AnalysisEvent
import EventSelection
from CPconfig import configuration, eventDumpConfig

def DumpEventInfo(event=None, eventNumber=None, path="./"):
  """Dump informations about a given event"""
  # in case no event is provided, find it using eventNumber
  if event is None:
    if eventNumber is None:
      print "DumpEventInfo Error: either pass an event or an event number"
      return
    # find event based on run and event
    if os.path.isdir(path):
      dirList=os.listdir(path)
      files=[]
      for fname in dirList:
        files.append(path+fname)
    elif os.path.isfile(path):
      files=[path]
    else:
      files=[]
    events = AnalysisEvent(files)
    EventSelection.prepareAnalysisEvent(events)
    DumpEventInfo(events[eventNumber])
    return
  # run the producers when we want to print the outcome, and mute unneeded collections
  for product in eventDumpConfig.productsToPrint: getattr(event,product)
  for collection in eventDumpConfig.collectionsToHide: event.removeCollection(collection)
  # Now, we can go on with the printing.
  print event

if __name__=="__main__":
  import sys
  DumpEventInfo(eventNumber=int(sys.argv[1]), path=sys.argv[2])

