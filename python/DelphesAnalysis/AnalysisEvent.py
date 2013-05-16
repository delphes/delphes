from inspect import getargspec
from ROOT import TLorentzVector,TVector3,TChain,TClass,TDatabasePDG
from datetime import datetime
from collections import Iterable
from types import StringTypes
from os import path

class AnalysisEvent(TChain):
   """A class that complements fwlite::Events with analysis facilities.
      The class provides the following additional functionalities:
        1. instrumentation for event weight
             A set of weight classes can be defined, and the event weight
             is computed and cached using those.
        2. list of event products used in the analysis
             It makes the iteration faster by only enabling required branches.
        3. a list of "producers" of analysis high-level quantities
             It allows to run "analysis on demand", by automatically running
             the defined producers to fill the cache, and later use that one.
        4. a volatile dictionary
             It allows to use the event as an heterogenous container for
             any analysis product. The event is properly reset when iterating 
             to the next event.
   """

   def __init__(self, inputFiles = '', maxEvents=0):
     """Initialize the AnalysisEvent like a standard Event, plus additional features."""
     # initialization of base functionalities
     TChain.__init__(self,"Delphes","Delphes")
     if isinstance(inputFiles,Iterable) and not isinstance(inputFiles,StringTypes):
       for thefile in inputFiles:
         if path.isfile(thefile): 
           self.AddFile(thefile)
         else:
           print "Warning: ",thefile," do not exist."
     elif  isinstance(inputFiles,StringTypes):
       thefile = inputFiles
       if path.isfile(thefile): 
         self.AddFile(thefile)
       else:
         print "Warning: ",thefile," do not exist."
     else:
       print "Warning: invalid inputFiles"
     self.BuildIndex("Event[0].Number")
     self.SetBranchStatus("*",0)
     self._eventCounts        = 0
     self._maxEvents          = maxEvents
     # additional features:
     # 1. instrumentation for event weight
     self._weightCache = {}
     self._weightEngines = {}
     # 2. a list of event products used in the analysis
     self._collections = {} 
     self._branches = dict((b,False) for b in map(lambda b:b.GetName(),self.GetListOfBranches()))
     # 3. a list of "producers" of analysis high-level quantities
     self._producers = {}
     # 4. volatile dictionary. User can add any quantity to the event and it will be 
     #    properly erased in the iteration step. 
     self.__dict__["vardict"] = {}
     
   def addWeight(self, name, weightClass):
     """Declare a new class (engine) to compute the weights.
        weightClass must have a weight() method returning a float."""
     if name in self._weightEngines:
       raise KeyError("%s weight engine is already declared" % name)
     self._weightEngines[name] = weightClass
     self._weightCache.clear()

   def delWeight(self, name):
     """Remove one weight engine from the internal list."""
     # just to clean the dictionnary
     del self._weightEngines[name]
     self._weightCache.clear()

   def weight(self, weightList=None, **kwargs):
     """Return the event weight. Arguments:
         * weightList is the list of engines to use, as a list of strings. 
              Default: all defined engines.
         * the other named arguments are forwarded to the engines.
        The output is the product of the selected individual weights."""
     # first check in the cache if the result is there already
     if weightList is None:
       weightList=self._weightEngines.keys()
     kwargs["weightList"] = weightList
     # compute the weight or use the cached value
     myhash = self._dicthash(kwargs)
     if not myhash in self._weightCache : 
       w = 1.
       for weightElement in weightList:
         engine = self._weightEngines[weightElement]
         engineArgs = getargspec(engine.weight).args
         subargs = dict((k,v) for k,v in kwargs.iteritems() if k in engineArgs)
         w *= self._weightCache.setdefault("weightElement:%s # %s" %(weightElement,self._dicthash(subargs)),engine.weight(self,**subargs))
       self._weightCache[myhash] = w
     return self._weightCache[myhash]

   def addCollection(self, name, inputTag):
     """Register an event collection as used by the analysis.
        Example: addCollection("myjets","jets") 
        Note that the direct access to the branch is still possible but unsafe."""
     if name in self._collections:
       raise KeyError("%r collection is already declared", name)
     if name in self._producers:
       raise KeyError("%r is already declared as a producer", name)
     if hasattr(self,name):
       raise AttributeError("%r object already has attribute %r" % (type(self).__name__, name))
     if inputTag not in self._branches:
       raise AttributeError("%r object has no branch %r" % (type(self).__name__, inputTag))
     self._collections[name] = inputTag
     self.SetBranchStatus(inputTag+"*",1)
     self._branches[inputTag] = True

   def removeCollection(self,name):
     """Forget about the named event collection.
        This method will delete both the product from the cache (if any) and the definition.
        To simply clear the cache, use "del event.name" instead. """
     self.SetBranchStatus(self._collections[name]+"*",0)
     self._branches[self._collections[name]] = False
     del self._collections[name]
     if name in self.vardict:
       delattr(self,name)

   def getCollection(self,name):
     """Retrieve the event product or return the cached collection.
        Note that the prefered way to get the collection is instead to access the "event.name" attribute."""
     if not name in self._collections:
       raise AttributeError("%r object has no attribute %r" % (type(self).__name__,name))
     if not name in self.vardict:
       self.vardict[name] = TChain.__getattr__(self,self._collections[name])
     return getattr(self,name)

   def addProducer(self,name,producer,**kwargs):
     """Register a producer to create new high-level analysis objects."""
     # sanity checks
     if name in self._producers:
       raise KeyError("%r producer is already declared", name)
     if name in self._collections:
       raise KeyError("%r is already declared as a collection", name)
     if hasattr(self,name):
       raise AttributeError("%r object already has attribute %r" % (type(self).__name__, attr))
     # remove name and producer from kwargs
     if "name" in kwargs: del kwargs["name"]
     if "producer" in kwargs: del kwargs["producer"]
     # store
     self._producers[name] = (producer,kwargs)

   def removeProducer(self,name):
     """Forget about the producer.
        This method will delete both the product from the cache (if any) and the producer.
        To simply clear the cache, use "del event.name" instead."""
     del self._producers[name]
     if name in self.vardict:
       delattr(self,name)

   def event(self):
     """Event number"""
     if self._branches["Event"]:
       return self.Event.At(0).Number
     else:
       return 0

   def to(self,event):
     """Jump to some event"""
     self.GetEntryWithIndex(event)
 
   def __getitem__(self,index):
     """Jump to some event"""
     self.GetEntryWithIndex(index)
     return self

   def __iter__ (self):
     """Iterator"""
     self._eventCounts = 0
     while self.GetEntry(self._eventCounts):
       self.vardict.clear()
       self._weightCache.clear()
       yield self
       self._eventCounts += 1
       if self._maxEvents > 0 and self._eventCounts >= self._maxEvents:
         break

   def __getattr__(self, attr):
     """Overloaded getter to handle properly:
          - volatile analysis objects
          - event collections
          - data producers"""
     if attr in self.__dict__["vardict"]:
       return self.vardict[attr]
     if attr in self._collections:
       return self.vardict.setdefault(attr, TChain.__getattr__(self,self._collections[attr]))
     if attr in self._producers:
       return self.vardict.setdefault(attr, self._producers[attr][0](self, **self._producers[attr][1]))
     return TChain.__getattr__(self,attr)

   def __setattr__(self, name, value):
     """Overloaded setter that puts any new attribute in the volatile dict."""
     if name in self.__dict__  or not "vardict" in self.__dict__ or name[0]=='_':
       self.__dict__[name] = value
     else:
       if name in self._collections or name in self._producers:
         raise AttributeError("%r object %r attribute is read-only (event collection)" % (type(self).__name__, name))
       self.vardict[name] = value

   def __delattr__(self, name):
     """Overloaded del method to handle the volatile internal dictionary."""
     if name=="vardict":
       raise AttributeError("%r object has no attribute %r" % (type(self).__name__, name))
     if name in self.__dict__:
       del self.__dict__[name]
     elif name in self.vardict:
       del self.vardict[name]
     else:
       raise AttributeError("%r object has no attribute %r" % (type(self).__name__, name))

   def _dicthash(self,dict):
     return (lambda d,j='=',s=';': s.join([j.join((str(k),str(v))) for k,v in d.iteritems()]))(dict)

   def __str__(self):
     """Event text dump."""
     dictjoin = lambda d,j=' => ',s='\n': s.join([j.join((str(k),str(v))) for k,v in d.iteritems()])
     mystring = "=================================================================\n"
     # general information
     if self._branches["Event"]:
       mystring += str(self.Event.At(0))
     else:
       mystring += "Event %d\n" % self.GetReadEvent()
     mystring += "-----------------------------------------------------------------\n"
     # weights
     if len(self._weightCache)==0: 
       mystring += "No weight computed so far. Default weight is %f.\n" % self.weight()
     else:
       mystring += "Weights:\n"
       mystring += dictjoin(self._weightCache)
     mystring += "\n-----------------------------------------------------------------\n"
     # list the collections
     mystring += "Collections:\n"
     for colname in self._collections.keys():
       collection = self.getCollection(colname)
       if collection.GetEntries()>0:
         if collection.At(0).IsA()==TClass.GetClass("HepMCEvent"):
           pass
         else:
           mystring += "*** %s has %d element(s)\n" % (colname,collection.GetEntries())
           mystring += reduce(lambda a,b: a+b,map(str,collection))
     mystring += "\n-----------------------------------------------------------------\n"
     # list the registered producers
     mystring += "Producers:\n"
     mystring += dictjoin(self._producers)
     mystring += "\n-----------------------------------------------------------------\n"
     # list the content of vardict, excluding collections
     mystring += "Content of the cache:\n"
     for k, v in self.vardict.iteritems():
       if k in self._collections.keys() : continue
       if isinstance(v,Iterable) and not isinstance(v,StringTypes):
         try:
           thisstring = "%s => vector of %d objects(s)\n" % (k,len(v))
         except:
           mystring += "%s => %s\n"%(k,str(v))
         else:
           try:
             for it,vec in enumerate(v):
               thisstring += "%s[%d] = %s\n"%(k,it,str(vec))
           except: 
             mystring += "%s => %s\n"%(k,str(v))
           else:
             mystring += thisstring
       else:
         mystring += "%s => %s\n"%(k,str(v))
     return mystring

   def decayTree(self, genparticles):
     db = TDatabasePDG()
     theString = ""
     for part in genparticles:
       if part.M1==-1 and part.M2==-1:
         theString += part.printDecay(db, genparticles)
     return theString

