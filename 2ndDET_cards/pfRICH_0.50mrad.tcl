module IdentificationMap pfRICH {
  set InputArray TrackMerger/tracks
  set OutputArray tracks

    add EfficiencyFormula {211} {211} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.972944) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999999) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.999422) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.992784) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.984919) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.958603) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.919817) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.875145) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.848991) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.822742) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.780192) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.995856) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999999) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.999407) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.995783) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.989594) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.968762) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.935193) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.892044) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.867138) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.828418) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.799323) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.970552) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999999) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.999955) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.997901) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.989495) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.967230) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.934204) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.900894) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.864474) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.822755) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.792580)
    }

    add EfficiencyFormula {211} {321} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.027056) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000001) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000578) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.007216) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.015081) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.041397) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.080183) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.124855) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.151009) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.177258) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.219808) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.004144) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000001) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000593) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.004217) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.010406) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.031238) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.064807) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.107956) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.132862) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.171582) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.200677) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.029448) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000001) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000045) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.002099) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.010505) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.032770) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.065796) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.099106) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.135526) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.177245) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.207420)
    }

    add EfficiencyFormula {321} {211} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000046) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000043) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.007470) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.024189) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.034587) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.094240) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.081035) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.149615) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.174369) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.211734) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000046) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000098) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.005179) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.018689) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.034477) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.056836) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.088863) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.117858) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.180251) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.214743) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000046) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000086) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.003942) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.007979) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.020159) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.063498) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.092283) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.145602) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.175805) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.208822)
    }

    add EfficiencyFormula {321} {321} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.999934) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999479) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.944324) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.975555) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.967251) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.965411) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.905604) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.918842) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.848052) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.812466) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.768801) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.999903) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999948) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.976324) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.984454) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.980779) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.965511) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.943060) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.910624) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.879356) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.812961) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.774542) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.991148) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999948) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.972872) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.993617) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.991955) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.979841) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.936424) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.906731) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.851930) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.817517) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.780553)
    }

    add EfficiencyFormula {321} {2212} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000066) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000475) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.055633) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.016974) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.008560) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000002) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000157) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000123) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.002332) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.013165) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.019465) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000097) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000006) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.023578) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.010367) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.000532) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000012) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000104) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000513) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.002787) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.006788) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.010714) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.008852) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000006) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.027042) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.002442) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.000065) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000079) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000985) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.002468) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.006679) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.010625)
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
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000016) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000268) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.001242) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.003875) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.012422) +
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
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000004) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000130) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.000604) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.002926) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.007275) +
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
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000005) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000092) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.000668) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.003503) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.008820)
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
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.999984) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.999732) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.998758) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.996125) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.987578) +
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
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.999996) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.999870) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.999396) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.997074) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.992725) +
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
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.999995) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.999908) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.999332) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.996497) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.991180)
    }

  add EfficiencyFormula {0} {0} { 0.00 }
}
