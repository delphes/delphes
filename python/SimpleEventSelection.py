#very simple EventSelection class aimed at demonstrating the 
#typical implementation of an EventSelection class

#On purpose, the categories are not in any "logical" order, and there is some redundant check when testing the category.
#This is mostly for illustration.

# requirements:
#   event.muons
#   event.electrons
#   event.jets

# the list of category names
categoryNames = [ "MuonChannel/SingleMuon", "ElectronChannel/SingleElectron", "ElectronChannel/Jet", "MuonChannel/DoubleMuon", "ElectronChannel/DoubleElectron" ]

def eventCategory(event):
  """Check analysis requirements for various steps
     and return a tuple of data used to decide 
     to what category an event belong """
  categoryData = [ ]
  # 0: number of muons
  categoryData.append(event.muons.GetEntries())
  # 1: number of electrons
  categoryData.append(event.electrons.GetEntries())
  # 2: number of jets
  categoryData.append(event.jets.GetEntries())
  # 3: Pt of the leading muon > 2GeV
  if event.muons.GetEntries():
    categoryData.append(event.muons[0].PT>1.)
  else:
    categoryData.append(False)
  # 4: Pt of the leading electron > 2 GeV
  if event.electrons.GetEntries():
    categoryData.append(event.electrons[0].PT>1.)
  else:
    categoryData.append(False)
  # DONE
  return categoryData

def isInCategory(category, categoryData):
  """Check if the event enters category X, given the tuple computed by eventCategory."""
  if category==0:
    return categoryData[0]>0 
  elif category==1:
    return categoryData[1]>0
  elif category==2:
    return categoryData[2]>0
  elif category==3:
    return isInCategory(0,categoryData) and categoryData[3]==True
  elif category==4:
    return isInCategory(1,categoryData) and categoryData[4]==True
  else:
    return False

