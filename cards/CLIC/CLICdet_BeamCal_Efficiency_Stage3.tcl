# Efficiency formula for photons/electrons in BeamCal
# Based on CLICdp-Note-2018-005 (Fig. 69)
# A.F.Zarnecki

set EfficiencyFormula {
    (energy  < 400 )                                                         * (0.00) +
    (energy >= 400 )                 * (abs(eta) <= 4.00 )                   * (0.00) +
    (energy >= 400 && energy < 800 ) * (abs(eta) > 4.00 && abs(eta) <= 4.15) * (0.00144 * energy - 0.16) +
    (energy >= 800 )                 * (abs(eta) > 4.00 && abs(eta) <= 4.15) * (0.99) +
    (energy >= 400 && energy < 800 ) * (abs(eta) > 4.15 && abs(eta) <= 4.70) * (0.99 + (0.0016*energy-1.2796290)*(abs(eta)-4.15)) +
    (energy >= 800 )                 * (abs(eta) > 4.15 && abs(eta) <= 4.70) * (0.99) +
    (energy >= 400 && energy < 800 ) * (abs(eta) > 4.70 && abs(eta) <= 5.30) * (0.00) +
    (energy >= 800 )                 * (abs(eta) > 4.70 && abs(eta) <= 5.30) *
    (1.- pow(1.-min(pow((5.35-abs(eta))/0.95,0.55),1.0),max(0.0029*energy-0.13,1.0))) * (0.98) +
    (energy >= 400 )                 * (abs(eta) > 5.30 )                    * (0.00)
  }
