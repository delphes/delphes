
  set DeltaR 0.1

  set TauPTMin 1.0

  set TauEtaMax 4.0

  # default efficiency formula (misidentification rate)
  add EfficiencyFormula {0} {
     (abs(eta) <= 3 ) * (energy <= 15 )                *   0.0000 +
     (abs(eta) <= 3 ) * (energy  > 15 && energy <=25 ) *  (0.0600 - 0.00336 * (energy - 15 ) ) +
     (abs(eta) <= 3 ) * (energy  > 25 && energy <=45 ) *  (0.0264 - 0.00041 * (energy - 25 ) ) +
     (abs(eta) <= 3 ) * (energy  > 45 )                *   0.0182 +
     (abs(eta)  > 3 )                                  *   0.0000 
  }

  # efficiency formula for tau-jets
  add EfficiencyFormula {15} {
     (abs(eta) <= 3 ) * (energy <=15 )                 *   0.000 +
     (abs(eta) <= 3 ) * (energy  > 15 && energy <=25 ) *  (0.118 + 0.0352 * (energy - 15 ) ) +
     (abs(eta) <= 3 ) * (energy  > 25 && energy <=45 ) *  (0.470 + 0.0090 * (energy - 25 ) ) +
     (abs(eta) <= 3 ) * (energy  > 45 )                *   0.650 +
     (abs(eta)  > 3 )                                  *   0.000 
  }
