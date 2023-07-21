module IdentificationMap pfRICH {
  set InputArray TrackMerger/tracks
  set OutputArray tracks

    add EfficiencyFormula {211} {211} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.969438) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.999618) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.994424) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.988458) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.965391) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.928925) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.884513) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.860498) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.835291) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.791175) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.994687) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.999612) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.997016) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.992627) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.975363) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.944895) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.903183) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.880329) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.841614) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.812721) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.966849) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.999983) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.998702) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.992562) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.973777) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.943895) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.912760) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.877449) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.835307) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.805130)
    }

    add EfficiencyFormula {211} {321} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.030562) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000382) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.005576) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.011542) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.034609) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.071075) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.115487) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.139502) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.164709) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.208825) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.005313) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000388) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.002984) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.007373) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.024637) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.055105) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.096817) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.119671) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.158386) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.187279) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.033151) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000017) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.001298) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.007438) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.026223) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.056105) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.087240) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.122551) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.164693) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.194870)
    }

    add EfficiencyFormula {321} {211} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000027) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000017) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.005802) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.020117) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.027898) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.085784) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.068343) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.137977) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.161484) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.199754) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000027) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000046) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.003797) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.014908) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.027780) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.046924) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.076374) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.103153) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.168062) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.203182) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000027) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000039) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.002762) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.005290) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.014290) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.053736) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.079947) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.133608) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.163086) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.196454)
    }

    add EfficiencyFormula {321} {321} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.999900) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999285) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.939280) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.975657) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.970480) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.972101) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.914117) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.931599) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.860423) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.827644) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.784389) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.999841) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999963) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.973590) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.984528) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.984522) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.972213) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.953012) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.923305) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.894855) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.826918) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.789004) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.988810) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999963) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.967904) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.994326) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.994648) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.985710) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.946219) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.919361) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.864669) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.832001) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.795820)
    }

    add EfficiencyFormula {321} {2212} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000100) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000688) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.060703) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.018541) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.009403) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000001) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000099) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000058) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.001600) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.010872) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.015856) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000159) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000010) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.026364) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.011675) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.000570) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000007) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000063) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000320) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.001992) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.005020) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.007814) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.011190) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000010) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.032057) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.002911) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.000063) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000046) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000691) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.001723) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.004913) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.007725)
    }

    add EfficiencyFormula {2212} {321} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000007) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000148) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.000757) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.002584) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.009338) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000001) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000063) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.000318) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.001823) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.004863) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000002) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000042) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.000359) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.002267) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.006159)
    }

    add EfficiencyFormula {2212} {2212} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.999993) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.999852) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.999243) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.997416) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.990662) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.999999) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.999937) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.999682) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.998177) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.995137) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.999998) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.999958) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.999641) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.997733) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.993841)
    }

  add EfficiencyFormula {0} {0} { 0.00 }
}
