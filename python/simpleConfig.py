#configuration of the ControlPlot machinery

from collections import namedtuple
controlPlot     = namedtuple("controlPlot",    ["label","module","classname","kwargs"])
eventCollection = namedtuple("eventCollection",["label","collection"])
eventProducer   = namedtuple("eventProducer",  ["label","module","function","kwargs"])
eventWeight     = namedtuple("eventWeight",    ["label","module","classname","kwargs"])

class configuration:
  # default I/O
  defaultFilename = "controlPlots"
  RDSname = "rds_delphes"
  WSname = "workspace_ras"

  # mode: plots or dataset
  runningMode = "plots"
  #runningMode = "dataset"

  # event selection class
  eventSelection = "SimpleEventSelection"

  # control plot classes
  controlPlots = [ controlPlot("selection","EventSelectionControlPlots","EventSelectionControlPlots", { }),
                   controlPlot("leptons","LeptonControlPlots","LeptonControlPlots", { }) ]

  # event content: lists of eventCollection, eventProducer, and eventWeight objects respectively.
  eventCollections = [ eventCollection("genEvent","Event"),
                       eventCollection("muons","EFlowMuon"), 
                       eventCollection("electrons","Electron"), 
                       eventCollection("jets","Jet"),
                       eventCollection("tracks","Track")] 
  eventProducers   = [ eventProducer("category","SimpleEventSelection","eventCategory",{ }) ] 
  eventWeights     = [ ]

class eventDumpConfig:
  # fine-tuning of the event content for display
  productsToPrint   = [ "category" ] # list of product to display (use the producer label)
  collectionsToHide = [ "tracks" ] # collections used in the analysis but not printed (use the collection label) 

