# efficiency formula for muons
# Identification based on tracking only

  set EfficiencyFormula {
        (energy <= 2 )                  * (0.00) +
        (energy > 2  && energy <= 10 )  * (0.95) +
        (energy > 10 )                  * (0.97)
  }
