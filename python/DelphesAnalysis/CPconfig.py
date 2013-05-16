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

  # event selection class
  eventSelection = ""

  # control plot classes
  controlPlots = [ ]

  # event content: lists of eventCollection, eventProducer, and eventWeight objects respectively.
  eventCollections = [ ]
  eventProducers   = [ ]
  eventWeights     = [ ]

class eventDumpConfig:
  # fine-tuning of the event content for display
  productsToPrint   = [ ] # list of product to display (use the producer label)
  collectionsToHide = [ ] # collections used in the analysis but not printed (use the collection label) 

# import the actual implementation of the configuration
import os
theConfig = os.getenv("DelphesAnalysisCfg")
if theConfig is not None:
  configImplementation = __import__(os.path.splitext(theConfig)[0])
  configuration = configImplementation.configuration
  eventDumpConfig = configImplementation.eventDumpConfig

