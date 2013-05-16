from CPconfig import configuration
from importlib import import_module

# A proper implementation of the event selection will at least 
# define the following:
# - categoryNames (list of names of categories)
# - eventCategory (category data producer)
# - isInCategory (category definition, from category data)

# the list of category names
categoryNames = [ ]

def eventCategory(event):
  """Check analysis requirements for various steps
     and return a tuple of data used to decide 
     to what category an event belong """
  raise NotImplementedError
  return []

def isInCategory(category, categoryData):
  """Check if the event enters category X, given the tuple computed by eventCategory."""
  raise NotImplementedError
  return True

# specific implementation of the "virtual methods" above
EventSelectionImplementation = __import__(configuration.eventSelection)
categoryNames  = EventSelectionImplementation.categoryNames
eventCategory  = EventSelectionImplementation.eventCategory
isInCategory   = EventSelectionImplementation.isInCategory

# Functions below should not be touched in any implementation.

def eventCategories():
  """Return the number of catefories """
  return len(categoryNames)

def categoryName(category):
  """Return the name of category x"""
  if category<eventCategories() and category>=0: return categoryNames[category]
  else: return "None"

def categoriesHierarchy():
  """Return a hierarchically organized list of categories"""
  dct = {}
  for i,item in enumerate(categoryNames):
    p = dct
    splitted = item.split('/')
    for x in splitted[:-1]:
      p = p.setdefault(x,{})
    p[splitted[-1]] = i
  return dct

def prepareAnalysisEvent(event):
  """Define collections and producers"""
  for coll in configuration.eventCollections:
    event.addCollection(coll.label,coll.collection)
  for prod in configuration.eventProducers:
    event.addProducer(prod.label,getattr(import_module(prod.module),prod.function),**prod.kwargs)
  for weight in configuration.eventWeights:
    events.addWeight(weight.label,getattr(import_module(weight.module),weight.classname)(**weight.kwargs))

