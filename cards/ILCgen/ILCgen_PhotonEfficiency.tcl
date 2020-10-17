# set EfficiencyFormula (efficiency formula as a function of eta and pt)
# efficiency formula for photons (ECAL + LumiCal)

  set EfficiencyFormula {
        (energy <= 2 )                                      * (0.00) +
        (energy > 2 ) * (abs(eta) <= 3.0)                   * (0.95) +
        (energy > 2 ) * (abs(eta) > 3.0 && abs(eta) <= 4.0) * (0.90)
  }
