#configuration of the ControlPlot machinery

from collections import namedtuple
controlPlot     = namedtuple("controlPlot",    ["label","module","classname","kwargs"])
eventCollection = namedtuple("eventCollection",["label","collection"])
eventProducer   = namedtuple("eventProducer",  ["label","module","function","kwargs"])
eventWeight     = namedtuple("eventWeight",    ["label","module","classname","kwargs"])

class configuration:
  # default I/O
  defaultFilename = "controlPlots_ttbar"
  RDSname = "rds_delphes_ttbar"
  WSname = "workspace_ras_ttbar"

  # mode: plots or dataset
  runningMode = "plots"
  #runningMode = "dataset"

  # event selection class
  eventSelection = "TtbarEventSelection"

  # control plot classes
  controlPlots = [ controlPlot("selection","EventSelectionControlPlots","EventSelectionControlPlots", { }),
                   controlPlot("leptons","LeptonControlPlots","LeptonControlPlots", { }),
                   controlPlot("jets","JetControlPlots","JetControlPlots", { }),
                   controlPlot("top","TopControlPlots","TopControlPlots", { }) 
                 ]

  # event content: lists of eventCollection, eventProducer, and eventWeight objects respectively.
  eventCollections = [ eventCollection("genEvent","Event"),
                       eventCollection("muons","Muon"), 
                       eventCollection("electrons","Electron"), 
                       eventCollection("jets","Jet"),
                       eventCollection("MEt","MissingET")
                     ] 
  eventProducers   = [ eventProducer("category","TtbarEventSelection","eventCategory",{ }),
                       eventProducer("selectedJets","TopReconstruction","jetSelection",{ "ptcut":20., "etacut":2.4 }),
                       eventProducer("bJets","TopReconstruction","bjets",{}),
                       eventProducer("lJets","TopReconstruction","ljets",{}),
                       eventProducer("topCandidates_L","TopReconstruction","topCandidates",{"leptonic":True,"hadronic":False}),
                       eventProducer("topCandidates_H","TopReconstruction","topCandidates",{"leptonic":False,"hadronic":True}),
                       eventProducer("topPairCandidates","TopReconstruction","topCandidates",{"leptonic":True,"hadronic":True})
                     ] 
  eventWeights     = [ ]

class eventDumpConfig:
  # fine-tuning of the event content for display
  productsToPrint   = [ "category", "topPairCandidates" ] # list of product to display (use the producer label)
  collectionsToHide = [ ] # collections used in the analysis but not printed (use the collection label) 

