module IdentificationMap pfRICH {
  set InputArray TrackMerger/tracks
  set OutputArray tracks

    add EfficiencyFormula {211} {211} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.970369) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999999) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.999575) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.994032) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.987604) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.963700) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.926613) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.882131) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.857523) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.832023) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.788312) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.995014) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999999) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.999567) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.996729) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.991910) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.973736) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.942443) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.900329) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.876915) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.838175) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.809210) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.967831) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999999) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.999978) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.998523) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.991836) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.972163) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.941445) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.909718) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.874091) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.832039) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.801846)
    }

    add EfficiencyFormula {211} {321} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.029631) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000001) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000425) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.005968) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.012396) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.036300) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.073387) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.117869) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.142477) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.167977) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.211688) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.004986) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000001) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000433) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.003271) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.008090) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.026264) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.057557) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.099671) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.123085) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.161825) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.190790) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.032169) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000001) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000022) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.001477) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.008164) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.027837) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.058555) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.090282) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.125909) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.167961) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.198154)
    }

    add EfficiencyFormula {321} {211} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000031) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000022) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.006201) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.021119) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.029553) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.087934) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.071577) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.140986) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.164841) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.202885) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000031) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000056) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.004122) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.015829) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.029438) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.049419) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.079569) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.106961) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.171233) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.206199) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000031) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.000048) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.003036) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.005915) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.015701) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.056202) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.083107) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.136709) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.166399) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.199690)
    }

    add EfficiencyFormula {321} {321} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.999910) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999344) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.940613) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.975679) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.969714) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.970445) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.911954) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.928351) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.857246) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.823731) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.780377) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.999860) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999960) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.974324) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.984561) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.983613) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.970554) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.950508) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.920068) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.890863) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.823329) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.785297) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (1.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.989460) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.999960) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.969270) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.994186) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.994022) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.984299) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.943745) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.916135) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.861396) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.828271) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.791896)
    }

    add EfficiencyFormula {321} {2212} {
      (eta< -3.80 || eta>= -1.50 || pt * cosh(eta) <    0.90 || pt * cosh(eta) >=   15.50) * ( 0.00 ) +
      ( -3.80 <= eta && eta <  -2.80) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -3.80 <= eta && eta <  -2.80) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -3.80 <= eta && eta <  -2.80) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000090) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000625) +
      ( -3.80 <= eta && eta <  -2.80) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.059365) +
      ( -3.80 <= eta && eta <  -2.80) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.018120) +
      ( -3.80 <= eta && eta <  -2.80) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.009168) +
      ( -3.80 <= eta && eta <  -2.80) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000001) +
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000112) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000071) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.001768) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.011428) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.016738) +
      ( -2.80 <= eta && eta <  -1.90) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -2.80 <= eta && eta <  -1.90) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -2.80 <= eta && eta <  -1.90) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.000140) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000009) +
      ( -2.80 <= eta && eta <  -1.90) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.025620) +
      ( -2.80 <= eta && eta <  -1.90) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.011317) +
      ( -2.80 <= eta && eta <  -1.90) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.000558) +
      ( -2.80 <= eta && eta <  -1.90) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000008) +
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000072) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000363) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.002176) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.005438) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.008503) +
      ( -1.90 <= eta && eta <  -1.50) * (   0.90 <= pt * cosh(eta) && pt * cosh(eta) <    2.90) * (0.500000) +
      ( -1.90 <= eta && eta <  -1.50) * (   2.90 <= pt * cosh(eta) && pt * cosh(eta) <    3.10) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.10 <= pt * cosh(eta) && pt * cosh(eta) <    3.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   3.50 <= pt * cosh(eta) && pt * cosh(eta) <    4.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   4.50 <= pt * cosh(eta) && pt * cosh(eta) <    5.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (   5.50 <= pt * cosh(eta) && pt * cosh(eta) <    6.10) * (0.010540) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.10 <= pt * cosh(eta) && pt * cosh(eta) <    6.50) * (0.000009) +
      ( -1.90 <= eta && eta <  -1.50) * (   6.50 <= pt * cosh(eta) && pt * cosh(eta) <    7.50) * (0.030682) +
      ( -1.90 <= eta && eta <  -1.50) * (   7.50 <= pt * cosh(eta) && pt * cosh(eta) <    8.50) * (0.002778) +
      ( -1.90 <= eta && eta <  -1.50) * (   8.50 <= pt * cosh(eta) && pt * cosh(eta) <    9.50) * (0.000063) +
      ( -1.90 <= eta && eta <  -1.50) * (   9.50 <= pt * cosh(eta) && pt * cosh(eta) <   10.50) * (0.000000) +
      ( -1.90 <= eta && eta <  -1.50) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000053) +
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000758) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.001895) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.005330) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.008415)
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
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000009) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000174) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.000865) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.002880) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.010077) +
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
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.000002) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000076) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.000378) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.002071) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.005423) +
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
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.000052) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.000425) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.002548) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.006785)
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
      ( -3.80 <= eta && eta <  -2.80) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.999991) +
      ( -3.80 <= eta && eta <  -2.80) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.999826) +
      ( -3.80 <= eta && eta <  -2.80) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.999135) +
      ( -3.80 <= eta && eta <  -2.80) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.997120) +
      ( -3.80 <= eta && eta <  -2.80) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.989923) +
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
      ( -2.80 <= eta && eta <  -1.90) * (  10.50 <= pt * cosh(eta) && pt * cosh(eta) <   11.50) * (0.999998) +
      ( -2.80 <= eta && eta <  -1.90) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.999924) +
      ( -2.80 <= eta && eta <  -1.90) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.999622) +
      ( -2.80 <= eta && eta <  -1.90) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.997929) +
      ( -2.80 <= eta && eta <  -1.90) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.994577) +
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
      ( -1.90 <= eta && eta <  -1.50) * (  11.50 <= pt * cosh(eta) && pt * cosh(eta) <   12.50) * (0.999948) +
      ( -1.90 <= eta && eta <  -1.50) * (  12.50 <= pt * cosh(eta) && pt * cosh(eta) <   13.50) * (0.999575) +
      ( -1.90 <= eta && eta <  -1.50) * (  13.50 <= pt * cosh(eta) && pt * cosh(eta) <   14.50) * (0.997452) +
      ( -1.90 <= eta && eta <  -1.50) * (  14.50 <= pt * cosh(eta) && pt * cosh(eta) <   15.50) * (0.993215)
    }

  add EfficiencyFormula {0} {0} { 0.00 }
}
