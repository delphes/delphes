# Efficiency formula for photons/electrons in LumiCal at 3 TeV
# Based on CLICdp-Note-2018-005 (Fig. 65)
# A.F.Zarnecki

set EfficiencyFormula {
    (energy < 10 ) *                                                       (0.00) +
    (energy >= 10 ) * (abs(eta) <= 3.0 ) *                                 (0.00) +
    (energy >= 10 && energy < 50 ) *  (abs(eta) > 3.00 && abs(eta) <= 3.22) *
                             pow(0.99 - 3.60283*(3.22-abs(eta)), 50./energy)      +
    (energy >= 50 ) *                 (abs(eta) > 3.00 && abs(eta) <= 3.22) *
                               (0.99 - 3.60283*(3.22-abs(eta)))                   +
    (energy >= 10 && energy < 60 ) *  (abs(eta) > 3.22 && abs(eta) <= 3.37) *
                                 (0.899781 +  0.0510312 * log10(energy))          +
    (energy >= 60 ) *                 (abs(eta) > 3.22 && abs(eta) <= 3.37) *
                                 (0.99)                                           +
    (energy >= 10 && energy < 98 ) *  (abs(eta) > 3.37 && abs(eta) <= 3.55) *
                                 (0.252111 +  0.37063 * log10(energy))            +
    (energy >= 98 ) *                 (abs(eta) > 3.37 && abs(eta) <= 3.55) *
                                 (0.99)                                           +
    (energy >= 10 && energy < 98 ) *  (abs(eta) > 3.55 && abs(eta) <= 4.00) *
     (1.- pow(1.-pow((4.4-abs(eta))/0.9,exp(2.60527-0.0139161*energy)),4.0)) * (0.252111 +  0.37063  * log10(energy)) +
    (energy >= 98 ) *                 (abs(eta) > 3.55 && abs(eta) <= 4.00) *
     (1.- pow(1.-pow((4.4-abs(eta))/0.9,exp(2.60527-0.0139161*energy)),4.0)) * (0.99) +
    (energy >= 10 ) * (abs(eta) > 4.0 ) *                                  (0.00)
  }
