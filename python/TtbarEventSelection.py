# basic selection for semi-leptonic ttbar events

# requirements:
#   event.muons
#   event.electrons
#   event.jets
#   event.selectedJets
#   event.bJets
#   event.topCandidates_L
#   event.topCandidates_H
#   event.topPairCandidates

# the list of category names
categoryNames = [ "l4j", "withPtCuts", "withBtag", "LeptonicTopCandidate", "HadronicTopCandidate", "TtbarCandidate" ]

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
  # 3: Pt of the leading muon 
  if event.muons.GetEntries():
    categoryData.append(event.muons[0].PT)
  else:
    categoryData.append(0.)
  # 4: Pt of the leading electron 
  if event.electrons.GetEntries():
    categoryData.append(event.electrons[0].PT)
  else:
    categoryData.append(0.)
  # 5: number of selected jets
  categoryData.append(len(event.selectedJets))
  # 6: number of bjets
  categoryData.append(len(event.bJets))
  # 7: leptonic top candidate
  categoryData.append(len(event.topCandidates_L))
  # 8: hadronic top candidate
  categoryData.append(len(event.topCandidates_H))
  # 9: top pairs
  categoryData.append(len(event.topPairCandidates))
  # DONE
  return categoryData

def isInCategory(category, categoryData):
  """Check if the event enters category X, given the tuple computed by eventCategory."""
  if category==0:
    return (categoryData[0]>0 or categoryData[1]>0) and categoryData[2]>=4
  elif category==1:
    return isInCategory(0,categoryData) and (categoryData[3]>15 or categoryData[4]>15) and categoryData[5]>=4
  elif category==2:
    return isInCategory(1,categoryData) and categoryData[6]>=2
  elif category==3:
    return isInCategory(2,categoryData) and categoryData[7]>0
  elif category==4:
    return isInCategory(2,categoryData) and categoryData[8]>0
  elif category==5:
    return isInCategory(2,categoryData) and categoryData[9]>0
  else:
    return False

