from itertools import combinations

def jetSelection(event, ptcut=20., etacut=2.4):
  def jetfilter(jet): return jet.PT>ptcut and abs(jet.Eta)<etacut
  return filter(jetfilter,event.jets)

def bjets(event):
  def jetfilter(jet): return jet.BTag
  return filter(jetfilter,event.selectedJets)

def ljets(event):
  def jetfilter(jet): return jet.BTag==0
  return filter(jetfilter,event.selectedJets)

def topCandidates(event,leptonic=True,hadronic=True):
  electrons = event.electrons
  muons = event.muons
  ljets = event.lJets
  bjets = event.bJets
  met = event.MEt[0]
  output = []
  # leptonic top candidates: lepton + bjet + MET
  if leptonic and not hadronic:
    # build all combinations
    for b in bjets:
      for l in electrons:
        output.append( (l,b,met) )
      for l in muons:
        output.append( (l,b,met) )
    # no specific order
    return output
  # hadronic top candidates: 2 jets + bjet
  if hadronic and not leptonic:
    # build all combinations
    for j in combinations(ljets,2):
      for b in bjets:
        output.append( (j[0], j[1], b) )
    # order by distance to top mass
    def massDistance(top):
      mass = (top[0].P4()+top[1].P4()+top[2].P4()).M()
      return abs(mass-172.9)
    return sorted(output,key=massDistance)
  # full event reconstruction
  if hadronic and leptonic:
    # build all combinations
    for j in combinations(ljets,2):
      for b in combinations(bjets,2):
        for l in electrons:
          output.append( (j[0], j[1], b[0], l, b[1], met) )
        for l in muons:
          output.append( (j[0], j[1], b[0], l, b[1], met) )
    # find the best ones: order via hadronic top mass
    def massDistance(top):
      mass = (top[0].P4()+top[1].P4()+top[2].P4()).M()
      return abs(mass-172.9)
    return sorted(output,key=massDistance)

