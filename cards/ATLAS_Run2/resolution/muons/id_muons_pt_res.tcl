set ResolutionFormula {
    (abs(eta) < 1.05) * sqrt((0.0044^2) + (0.00012*pt)^2) +
    (abs(eta) >= 1.05 && abs(eta) < 2.0) * sqrt((0.0067^2) + (0.00031*pt)^2) +
    (abs(eta) > 2.0) * sqrt((0.0094^2) + (0.00008*pt)^2)
}

