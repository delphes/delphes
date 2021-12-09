# P.Sopicki, corrected by A.F.Zarnecki: based on email from D.Jeans
# HCAL and LHCAL resolution: same formula used at the moment
#

set ResolutionFormula {
   (abs(eta) <= 2.8 )                 * sqrt(energy^2*0.017^2 + energy*0.45^2)+
   (abs(eta) > 2.8 && abs(eta)<=3.8 ) * sqrt(energy^2*0.017^2 + energy*0.45^2)
  }
