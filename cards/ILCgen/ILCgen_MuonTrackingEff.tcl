
# Tracking efficiency as given in the ILD IDR (Figure 8.2)
#  for high angles (FTD coverage): educated guess
#  A.F.Zarnecki, June 12, 2020

set EfficiencyFormula {
                                                                 (pt <= 0.1)   * (0.00)  +
                       (abs(eta) <= 1.83)     *      (pt > 0.1 && pt <= 0.2)   * (0.70)  +
                       (abs(eta) <= 1.83)     *      (pt > 0.2 && pt <= 0.3)   * (0.93)  +
                       (abs(eta) <= 1.83)     *      (pt > 0.3 && pt <= 0.7)   * (0.995) +
                       (abs(eta) <= 1.83)     *      (pt > 0.7)                * (1.00)  +
    (abs(eta) > 1.83 && abs(eta) <= 2.65)     *      (pt > 0.1 && pt <= 0.2)   * (0.697) +
    (abs(eta) > 1.83 && abs(eta) <= 2.65)     *      (pt > 0.2 && pt <= 0.3)   * (0.925) +
    (abs(eta) > 1.83 && abs(eta) <= 2.65)     *      (pt > 0.3 && pt <= 0.7)   * (0.99)  +
    (abs(eta) > 1.83 && abs(eta) <= 2.65)     *      (pt > 0.7)                * (0.995) +
    (abs(eta) > 2.65 && abs(eta) <= 3.00)     *      (pt > 0.1 && pt <= 0.2)   * (0.665) +
    (abs(eta) > 2.65 && abs(eta) <= 3.00)     *      (pt > 0.2 && pt <= 0.3)   * (0.884) +
    (abs(eta) > 2.65 && abs(eta) <= 3.00)     *      (pt > 0.3 && pt <= 0.7)   * (0.945) +
    (abs(eta) > 2.65 && abs(eta) <= 3.00)     *      (pt > 0.7)                * (0.95)  +
    (abs(eta) > 3.00)                         *      (pt > 0.1)                * (0.00)
  }
