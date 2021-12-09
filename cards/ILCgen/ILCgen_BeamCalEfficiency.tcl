# Efficiency formula for photons in Forward BeamCal
# Based on IDR Figure 8.8a and Moritz Hebermehl PhD Thesis Figure 4.10
# A.F.Zarnecki,  June 2020

set EfficiencyFormula {
    (energy <= 10 ) *                                                      (0.00) +
    (energy > 10 ) * (abs(eta)  > 5.88 )                                 * (0.00) +
    (energy > 10 ) * (abs(eta) <= 4.25 )                     * (0.98 - 1.0/energy) +
    (energy > 10 ) * (abs(eta) > 4.25 && abs(eta) <= 5.88)   *
         (1.- pow(1.-pow(3.6-abs(eta)/1.634,pow(85./energy,1.4)),4.0)) * (0.98 - 1.0/energy) 
  }
