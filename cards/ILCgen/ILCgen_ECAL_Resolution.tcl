# P.Sopicki, corrected by A.F.Zarnecki: based on email from D.Jeans
# ECAL and LumiCal resolution: same formula used at the moment
#

set ResolutionFormula {
    (abs(eta) <= 3 )                 * sqrt(energy^2*0.01^2 + energy*0.17^2)+
    (abs(eta) > 3 && abs(eta) <= 4 ) * sqrt(energy^2*0.01^2 + energy*0.17^2) 
  }
