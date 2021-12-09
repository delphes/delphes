
# efficiency formula for electrons (ECAL and LumiCal)
# LumiCal added for consistency: not in tracking range

  set EfficiencyFormula {
        (energy <= 2 )                                        * (0.00) +
        (energy > 2 )  *  (abs(eta) <= 3.0)                   * (0.95) +
        (energy > 2 )  *  (abs(eta) > 3.0 && abs(eta) <= 4.0) * (0.90) 
  }
