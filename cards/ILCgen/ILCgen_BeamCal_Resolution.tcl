# P.Sopicki, corrected by A.F.Zarnecki
# BeamCal resolution

  set ResolutionFormula {
  (abs(eta) > 4.0 && abs(eta) <= 4.8) * sqrt(energy^2*0.02^2 + energy*0.30^2) +
  (abs(eta) > 4.8 && abs(eta) <= 5.8) * sqrt(energy^2*0.03^2 + energy*0.45^2) 
  }
