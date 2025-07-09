# /author Jan T. Offermann - Brown University
#
# Based on the "trackResolutionCMS.tcl" card that ships with the official Delphes.
#
# It appears that the track resolution took values from arXiv:1405.6569, specifically
# corresponding with isolated muons. That paper also includes separate results for
# electrons and pions, so I'm including those as additional files.
#
# I've also adjusted the resolution formulae for P, CtgTheta and Phi to not be zero,
# otherwise this totally breaks the VertexFinderDA4D since it computes weighted
# positions that incorporate these errors (and having them be zero will result in
# NaNs).
#
# I've tried to cross-check the resolution formulae with the CMS Tracker TDR
# (https://cds.cern.ch/record/368412/), from my limited checks it seems sensible.
# I've extracted the data using the "WebPlotDigitizer" tool.
# Note that for the PResolutionFormula, I'm using the uncertainties on ctgTheta
# and pt, hence the slightly cumbersome formula.

# For PResolutionFormula, we need to get d(|p|) in terms of d(ctgTheta) and d(pt),
# and properly handle propagation of uncertainty.
#
# deta = d(arcsinh(cot(theta))) = 1 / sqrt(ctgTheta^2 + 1) * d(ctgTheta)
# p = pt * cosh(eta)
# => dp = d(pt * cosh(eta)) = \sqrt{ (pt * sinh(eta))^2 * deta^2 + cosh^2(eta) * dpt^2 }
# => dp^2 = pt^2 * sinh^2(eta) * deta^2 + cosh^2(eta) * dpt^2
#         = pt^2 * sinh^2(eta) / (ctgTheta^2 + 1) * d(ctgTheta)^2 + cosh^2(eta) * dpt^2
#         = pt^2 * sinh^2(eta) / (sinh^2(eta) + 1) * d(ctgTheta)^2 + cosh^2(eta) * dpt^2
#
#    (dp^2)/p^2 = [pt^2 * sinh^2(eta) / (sinh^2(eta) + 1) * d(ctgTheta)^2 + cosh^2(eta) * dpt^2] / [pt^2 cosh^2(eta)]
#               =         sinh^2(eta) / [(sinh^2(eta) + 1) * cosh^2(eta)] * d(ctgTheta)^2 + dpt^2 / pt^2
#
# Keep in mind whether using relative or absolute pt uncertainties! The plots used as
# reference give the relative uncertainties.
# Also, we ultimately want to calculate (dp)/p, since we need the relative |p| uncertainty.
set PResolutionFormula {
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.1823 && pt <=   0.2241 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  9.1530)**2 + (0.0149)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.2241 && pt <=   0.2737 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  7.3764)**2 + (0.0145)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.2737 && pt <=   0.3343 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  5.9395)**2 + (0.0135)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.3343 && pt <=   0.4083 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  5.0446)**2 + (0.0129)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.4083 && pt <=   0.4990 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  4.3076)**2 + (0.0119)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.4990 && pt <=   0.6099 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  3.5773)**2 + (0.0115)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.6099 && pt <=   0.7451 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  2.8155)**2 + (0.0105)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.7451 && pt <=   0.9100 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  2.2775)**2 + (0.0095)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.9100 && pt <=   1.1110 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  2.0115)**2 + (0.0092)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   1.1110 && pt <=   1.3567 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.7465)**2 + (0.0085)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   1.3567 && pt <=   1.6565 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.4459)**2 + (0.0082)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   1.6565 && pt <=   2.0232 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.2471)**2 + (0.0082)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   2.0232 && pt <=   2.4727 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.1132)**2 + (0.0085)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   2.4727 && pt <=   3.0212 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.9798)**2 + (0.0085)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   3.0212 && pt <=   3.6917 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.8467)**2 + (0.0085)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   3.6917 && pt <=   4.5097 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.7471)**2 + (0.0085)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   4.5097 && pt <=   5.5078 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.6465)**2 + (0.0089)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   5.5078 && pt <=   6.7269 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.6144)**2 + (0.0092)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   6.7269 && pt <=   8.2145 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.5477)**2 + (0.0089)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   8.2145 && pt <=  10.0282 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4810)**2 + (0.0088)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  10.0282 && pt <=  12.2447 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4811)**2 + (0.0092)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  12.2447 && pt <=  14.9651 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4148)**2 + (0.0095)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  14.9651 && pt <=  18.2986 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4150)**2 + (0.0099)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  18.2986 && pt <=  22.3600 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4149)**2 + (0.0102)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  22.3600 && pt <=  27.3128 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4150)**2 + (0.0105)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  27.3128 && pt <=  33.3396 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3815)**2 + (0.0112)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  33.3396 && pt <=  40.7156 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3818)**2 + (0.0122)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  40.7156 && pt <=  49.7409 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3819)**2 + (0.0139)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  49.7409 && pt <=  60.7340 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3484)**2 + (0.0149)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  60.7340 && pt <=  74.1628 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3486)**2 + (0.0169)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  74.1628 && pt <=  90.6169 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3482)**2 + (0.0209)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  90.6169 && pt <= 110.8043 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3483)**2 + (0.0236)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 110.8043 && pt <= 135.3741 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3472)**2 + (0.0272)**2 ) +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 135.3741                   ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3474)**2 + (0.0339)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.1823 && pt <=   0.2241 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 17.1477)**2 + (0.0257)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.2241 && pt <=   0.2737 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 14.7170)**2 + (0.0245)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.2737 && pt <=   0.3343 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 11.3702)**2 + (0.0208)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.3343 && pt <=   0.4083 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  8.9394)**2 + (0.0175)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.4083 && pt <=   0.4990 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  7.2088)**2 + (0.0169)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.4990 && pt <=   0.6099 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  5.9035)**2 + (0.0166)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.6099 && pt <=   0.7451 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  5.4346)**2 + (0.0162)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.7451 && pt <=   0.9100 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  4.3655)**2 + (0.0156)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.9100 && pt <=   1.1110 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  3.4808)**2 + (0.0146)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   1.1110 && pt <=   1.3567 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  2.9456)**2 + (0.0146)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   1.3567 && pt <=   1.6565 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  2.6056)**2 + (0.0150)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   1.6565 && pt <=   2.0232 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  2.1170)**2 + (0.0146)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   2.0232 && pt <=   2.4727 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.7932)**2 + (0.0146)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   2.4727 && pt <=   3.0212 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.5352)**2 + (0.0146)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   3.0212 && pt <=   3.6917 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.3427)**2 + (0.0150)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   3.6917 && pt <=   4.5097 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.2179)**2 + (0.0150)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   4.5097 && pt <=   5.5078 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.0504)**2 + (0.0149)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   5.5078 && pt <=   6.7269 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.9230)**2 + (0.0146)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   6.7269 && pt <=   8.2145 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.8241)**2 + (0.0149)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   8.2145 && pt <=  10.0282 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.7586)**2 + (0.0156)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  10.0282 && pt <=  12.2447 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.6287)**2 + (0.0159)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  12.2447 && pt <=  14.9651 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.5326)**2 + (0.0162)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  14.9651 && pt <=  18.2986 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4995)**2 + (0.0169)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  18.2986 && pt <=  22.3600 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.5010)**2 + (0.0172)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  22.3600 && pt <=  27.3128 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4694)**2 + (0.0173)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  27.3128 && pt <=  33.3396 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4365)**2 + (0.0176)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  33.3396 && pt <=  40.7156 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4365)**2 + (0.0182)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  40.7156 && pt <=  49.7409 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4033)**2 + (0.0195)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  49.7409 && pt <=  60.7340 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4051)**2 + (0.0212)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  60.7340 && pt <=  74.1628 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4373)**2 + (0.0235)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  74.1628 && pt <=  90.6169 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4051)**2 + (0.0275)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  90.6169 && pt <= 110.8043 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.4059)**2 + (0.0327)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 110.8043 && pt <= 135.3741 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3703)**2 + (0.0377)**2 ) +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 135.3741                   ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  0.3703)**2 + (0.0415)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.1823 && pt <=   0.2241 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 51.4418)**2 + (0.0353)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.2241 && pt <=   0.2737 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 47.6056)**2 + (0.0343)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.2737 && pt <=   0.3343 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 40.5330)**2 + (0.0330)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.3343 && pt <=   0.4083 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 32.9366)**2 + (0.0292)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.4083 && pt <=   0.4990 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 27.3163)**2 + (0.0276)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.4990 && pt <=   0.6099 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 22.4358)**2 + (0.0259)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.6099 && pt <=   0.7451 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 18.2310)**2 + (0.0229)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.7451 && pt <=   0.9100 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 15.2093)**2 + (0.0216)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.9100 && pt <=   1.1110 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 12.6790)**2 + (0.0206)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   1.1110 && pt <=   1.3567 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * ( 10.0390)**2 + (0.0199)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   1.3567 && pt <=   1.6565 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  8.4620)**2 + (0.0202)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   1.6565 && pt <=   2.0232 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  7.1652)**2 + (0.0199)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   2.0232 && pt <=   2.4727 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  6.0065)**2 + (0.0195)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   2.4727 && pt <=   3.0212 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  5.3033)**2 + (0.0199)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   3.0212 && pt <=   3.6917 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  4.3998)**2 + (0.0192)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   3.6917 && pt <=   4.5097 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  3.6882)**2 + (0.0188)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   4.5097 && pt <=   5.5078 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  3.2904)**2 + (0.0209)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   5.5078 && pt <=   6.7269 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  2.8455)**2 + (0.0212)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   6.7269 && pt <=   8.2145 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  2.4361)**2 + (0.0212)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   8.2145 && pt <=  10.0282 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  2.0422)**2 + (0.0219)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  10.0282 && pt <=  12.2447 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.8119)**2 + (0.0236)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  12.2447 && pt <=  14.9651 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.7472)**2 + (0.0246)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  14.9651 && pt <=  18.2986 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.5427)**2 + (0.0262)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  18.2986 && pt <=  22.3600 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.4496)**2 + (0.0276)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  22.3600 && pt <=  27.3128 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.4124)**2 + (0.0279)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  27.3128 && pt <=  33.3396 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.2799)**2 + (0.0312)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  33.3396 && pt <=  40.7156 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.2799)**2 + (0.0359)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  40.7156 && pt <=  49.7409 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.3135)**2 + (0.0400)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  49.7409 && pt <=  60.7340 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.2799)**2 + (0.0449)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  60.7340 && pt <=  74.1628 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.2466)**2 + (0.0540)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  74.1628 && pt <=  90.6169 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.2799)**2 + (0.0623)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  90.6169 && pt <= 110.8043 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.2799)**2 + (0.0705)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 110.8043 && pt <= 135.3741 ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.3135)**2 + (0.0883)**2 ) +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 135.3741                   ) * sqrt( TMath::SinH(eta)**2 / (TMath::CosH(eta)**2 * (TMath::SinH(eta)**2 + 1) ) * (  1.3480)**2 + (0.1049)**2 )
}

# For cross-check, consider TDR Figure 9.12, pg 484.
set CtgThetaResolutionFormula {
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.1823 && pt <=   0.2257 ) *   0.0092 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.2257 && pt <=   0.2754 ) *   0.0074 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.2754 && pt <=   0.3364 ) *   0.0059 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.3364 && pt <=   0.4110 ) *   0.0050 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.4110 && pt <=   0.5023 ) *   0.0043 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.5023 && pt <=   0.6135 ) *   0.0036 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.6135 && pt <=   0.7492 ) *   0.0028 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.7492 && pt <=   0.9149 ) *   0.0023 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.9149 && pt <=   1.1170 ) *   0.0020 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   1.1170 && pt <=   1.3642 ) *   0.0017 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   1.3642 && pt <=   1.6674 ) *   0.0014 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   1.6674 && pt <=   2.0380 ) *   0.0012 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   2.0380 && pt <=   2.4895 ) *   0.0011 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   2.4895 && pt <=   3.0415 ) *   0.0010 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   3.0415 && pt <=   3.7150 ) *   0.0008 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   3.7150 && pt <=   4.5363 ) *   0.0007 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   4.5363 && pt <=   5.5400 ) *   0.0006 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   5.5400 && pt <=   6.7654 ) *   0.0006 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   6.7654 && pt <=   8.2616 ) *   0.0005 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   8.2616 && pt <=  10.0935 ) *   0.0005 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  10.0935 && pt <=  12.3380 ) *   0.0005 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  12.3380 && pt <=  15.0746 ) *   0.0004 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  15.0746 && pt <=  18.4141 ) *   0.0004 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  18.4141 && pt <=  22.4780 ) *   0.0004 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  22.4780 && pt <=  27.4584 ) *   0.0004 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  27.4584 && pt <=  33.5437 ) *   0.0004 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  33.5437 && pt <=  40.9466 ) *   0.0004 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  40.9466 && pt <=  50.0500 ) *   0.0004 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  50.0500 && pt <=  61.1531 ) *   0.0003 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  61.1531 && pt <=  74.6915 ) *   0.0003 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  74.6915 && pt <=  91.2665 ) *   0.0003 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  91.2665 && pt <= 111.4703 ) *   0.0003 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 111.4703 && pt <= 136.1333 ) *   0.0003 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 136.1333                   ) *   0.0003 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.1823 && pt <=   0.2257 ) *   0.0171 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.2257 && pt <=   0.2754 ) *   0.0147 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.2754 && pt <=   0.3364 ) *   0.0114 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.3364 && pt <=   0.4110 ) *   0.0089 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.4110 && pt <=   0.5023 ) *   0.0072 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.5023 && pt <=   0.6135 ) *   0.0059 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.6135 && pt <=   0.7492 ) *   0.0054 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.7492 && pt <=   0.9149 ) *   0.0044 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.9149 && pt <=   1.1170 ) *   0.0035 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   1.1170 && pt <=   1.3642 ) *   0.0029 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   1.3642 && pt <=   1.6674 ) *   0.0026 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   1.6674 && pt <=   2.0380 ) *   0.0021 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   2.0380 && pt <=   2.4895 ) *   0.0018 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   2.4895 && pt <=   3.0415 ) *   0.0015 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   3.0415 && pt <=   3.7150 ) *   0.0013 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   3.7150 && pt <=   4.5363 ) *   0.0012 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   4.5363 && pt <=   5.5400 ) *   0.0011 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   5.5400 && pt <=   6.7654 ) *   0.0009 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   6.7654 && pt <=   8.2616 ) *   0.0008 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   8.2616 && pt <=  10.0935 ) *   0.0008 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  10.0935 && pt <=  12.3380 ) *   0.0006 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  12.3380 && pt <=  15.0746 ) *   0.0005 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  15.0746 && pt <=  18.4141 ) *   0.0005 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  18.4141 && pt <=  22.4780 ) *   0.0005 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  22.4780 && pt <=  27.4584 ) *   0.0005 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  27.4584 && pt <=  33.5437 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  33.5437 && pt <=  40.9466 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  40.9466 && pt <=  50.0500 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  50.0500 && pt <=  61.1531 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  61.1531 && pt <=  74.6915 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  74.6915 && pt <=  91.2665 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  91.2665 && pt <= 111.4703 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 111.4703 && pt <= 136.1333 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 136.1333                   ) *   0.0004 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.1823 && pt <=   0.2257 ) *   0.0514 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.2257 && pt <=   0.2754 ) *   0.0476 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.2754 && pt <=   0.3364 ) *   0.0405 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.3364 && pt <=   0.4110 ) *   0.0329 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.4110 && pt <=   0.5023 ) *   0.0273 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.5023 && pt <=   0.6135 ) *   0.0224 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.6135 && pt <=   0.7492 ) *   0.0182 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.7492 && pt <=   0.9149 ) *   0.0152 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.9149 && pt <=   1.1170 ) *   0.0127 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   1.1170 && pt <=   1.3642 ) *   0.0100 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   1.3642 && pt <=   1.6674 ) *   0.0085 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   1.6674 && pt <=   2.0380 ) *   0.0072 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   2.0380 && pt <=   2.4895 ) *   0.0060 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   2.4895 && pt <=   3.0415 ) *   0.0053 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   3.0415 && pt <=   3.7150 ) *   0.0044 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   3.7150 && pt <=   4.5363 ) *   0.0037 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   4.5363 && pt <=   5.5400 ) *   0.0033 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   5.5400 && pt <=   6.7654 ) *   0.0028 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   6.7654 && pt <=   8.2616 ) *   0.0024 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   8.2616 && pt <=  10.0935 ) *   0.0020 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  10.0935 && pt <=  12.3380 ) *   0.0018 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  12.3380 && pt <=  15.0746 ) *   0.0017 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  15.0746 && pt <=  18.4141 ) *   0.0015 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  18.4141 && pt <=  22.4780 ) *   0.0014 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  22.4780 && pt <=  27.4584 ) *   0.0014 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  27.4584 && pt <=  33.5437 ) *   0.0013 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  33.5437 && pt <=  40.9466 ) *   0.0013 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  40.9466 && pt <=  50.0500 ) *   0.0013 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  50.0500 && pt <=  61.1531 ) *   0.0013 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  61.1531 && pt <=  74.6915 ) *   0.0012 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  74.6915 && pt <=  91.2665 ) *   0.0013 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  91.2665 && pt <= 111.4703 ) *   0.0013 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 111.4703 && pt <= 136.1333 ) *   0.0013 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 136.1333                   ) *   0.0013
}

# For cross-check, consider TDR Figure 9.14, pg 485.
set PhiResolutionFormula {
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.1823 && pt <=   0.2265 ) *   0.0095 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.2265 && pt <=   0.2731 ) *   0.0072 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.2731 && pt <=   0.3229 ) *   0.0058 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.3229 && pt <=   0.3873 ) *   0.0049 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.3873 && pt <=   0.4730 ) *   0.0041 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.4730 && pt <=   0.5777 ) *   0.0034 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.5777 && pt <=   0.7057 ) *   0.0026 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.7057 && pt <=   0.8617 ) *   0.0020 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   0.8617 && pt <=   1.0532 ) *   0.0017 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   1.0532 && pt <=   1.2862 ) *   0.0015 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   1.2862 && pt <=   1.5700 ) *   0.0013 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   1.5700 && pt <=   1.9174 ) *   0.0011 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   1.9174 && pt <=   2.3423 ) *   0.0009 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   2.3423 && pt <=   2.8607 ) *   0.0007 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   2.8607 && pt <=   3.4925 ) *   0.0007 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   3.4925 && pt <=   4.2672 ) *   0.0006 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   4.2672 && pt <=   5.2150 ) *   0.0005 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   5.2150 && pt <=   6.3725 ) *   0.0004 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   6.3725 && pt <=   7.7827 ) *   0.0004 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   7.7827 && pt <=   9.5017 ) *   0.0003 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >   9.5017 && pt <=  11.6071 ) *   0.0003 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  11.6071 && pt <=  14.1733 ) *   0.0002 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  14.1733 && pt <=  17.3070 ) *   0.0002 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  17.3070 && pt <=  21.1331 ) *   0.0002 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  21.1331 && pt <=  25.8090 ) *   0.0001 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  25.8090 && pt <=  31.5393 ) *   0.0001 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  31.5393 && pt <=  38.5454 ) *   0.0001 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  38.5454 && pt <=  47.1046 ) *   0.0001 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  47.1046 && pt <=  57.5351 ) *   0.0001 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  57.5351 && pt <=  70.2519 ) *   0.0001 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  70.2519 && pt <=  85.8030 ) *   0.0001 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt >  85.8030 && pt <= 104.5942 ) *   0.0001 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 104.5942 && pt <= 127.7932 ) *   0.0001 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 127.7932                   ) *   0.0001 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.1823 && pt <=   0.2265 ) *   0.0116 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.2265 && pt <=   0.2731 ) *   0.0096 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.2731 && pt <=   0.3229 ) *   0.0080 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.3229 && pt <=   0.3873 ) *   0.0058 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.3873 && pt <=   0.4730 ) *   0.0050 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.4730 && pt <=   0.5777 ) *   0.0037 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.5777 && pt <=   0.7057 ) *   0.0033 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.7057 && pt <=   0.8617 ) *   0.0029 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   0.8617 && pt <=   1.0532 ) *   0.0023 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   1.0532 && pt <=   1.2862 ) *   0.0020 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   1.2862 && pt <=   1.5700 ) *   0.0015 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   1.5700 && pt <=   1.9174 ) *   0.0013 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   1.9174 && pt <=   2.3423 ) *   0.0011 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   2.3423 && pt <=   2.8607 ) *   0.0009 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   2.8607 && pt <=   3.4925 ) *   0.0008 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   3.4925 && pt <=   4.2672 ) *   0.0007 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   4.2672 && pt <=   5.2150 ) *   0.0006 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   5.2150 && pt <=   6.3725 ) *   0.0005 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   6.3725 && pt <=   7.7827 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   7.7827 && pt <=   9.5017 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >   9.5017 && pt <=  11.6071 ) *   0.0004 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  11.6071 && pt <=  14.1733 ) *   0.0003 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  14.1733 && pt <=  17.3070 ) *   0.0002 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  17.3070 && pt <=  21.1331 ) *   0.0002 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  21.1331 && pt <=  25.8090 ) *   0.0002 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  25.8090 && pt <=  31.5393 ) *   0.0002 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  31.5393 && pt <=  38.5454 ) *   0.0001 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  38.5454 && pt <=  47.1046 ) *   0.0001 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  47.1046 && pt <=  57.5351 ) *   0.0001 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  57.5351 && pt <=  70.2519 ) *   0.0001 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  70.2519 && pt <=  85.8030 ) *   0.0001 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt >  85.8030 && pt <= 104.5942 ) *   0.0001 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 104.5942 && pt <= 127.7932 ) *   0.0001 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 127.7932                   ) *   0.0001 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.1823 && pt <=   0.2265 ) *   0.0180 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.2265 && pt <=   0.2731 ) *   0.0152 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.2731 && pt <=   0.3229 ) *   0.0237 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.3229 && pt <=   0.3873 ) *   0.0132 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.3873 && pt <=   0.4730 ) *   0.0108 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.4730 && pt <=   0.5777 ) *   0.0084 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.5777 && pt <=   0.7057 ) *   0.0069 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.7057 && pt <=   0.8617 ) *   0.0057 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   0.8617 && pt <=   1.0532 ) *   0.0044 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   1.0532 && pt <=   1.2862 ) *   0.0038 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   1.2862 && pt <=   1.5700 ) *   0.0029 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   1.5700 && pt <=   1.9174 ) *   0.0024 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   1.9174 && pt <=   2.3423 ) *   0.0021 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   2.3423 && pt <=   2.8607 ) *   0.0018 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   2.8607 && pt <=   3.4925 ) *   0.0014 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   3.4925 && pt <=   4.2672 ) *   0.0012 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   4.2672 && pt <=   5.2150 ) *   0.0010 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   5.2150 && pt <=   6.3725 ) *   0.0009 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   6.3725 && pt <=   7.7827 ) *   0.0008 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   7.7827 && pt <=   9.5017 ) *   0.0006 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >   9.5017 && pt <=  11.6071 ) *   0.0005 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  11.6071 && pt <=  14.1733 ) *   0.0005 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  14.1733 && pt <=  17.3070 ) *   0.0004 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  17.3070 && pt <=  21.1331 ) *   0.0004 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  21.1331 && pt <=  25.8090 ) *   0.0003 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  25.8090 && pt <=  31.5393 ) *   0.0003 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  31.5393 && pt <=  38.5454 ) *   0.0003 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  38.5454 && pt <=  47.1046 ) *   0.0002 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  47.1046 && pt <=  57.5351 ) *   0.0002 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  57.5351 && pt <=  70.2519 ) *   0.0002 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  70.2519 && pt <=  85.8030 ) *   0.0002 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt >  85.8030 && pt <= 104.5942 ) *   0.0002 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 104.5942 && pt <= 127.7932 ) *   0.0002 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 127.7932                   ) *   0.0002
}

# taken from arXiv:1405.6569 fig. 15
set D0ResolutionFormula {
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 0.3543 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 0.2809 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 0.2304 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 0.1917 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 0.1737 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 0.1434 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.1060 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.0893 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.0753 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.0670 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.0577 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.0524 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.0452 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.0376 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.0350 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.0324 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.0283 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.0258 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.0237 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.0211 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.0191 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.0164 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.0150 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0143 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0130 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0130 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0116 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0116 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0110 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0110 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0110 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0104 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0109 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 134.7130 ) * 0.0110 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 0.4564 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 0.3580 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 0.3010 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 0.2353 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 0.2026 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 0.1595 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.1383 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.1119 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.0926 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.0816 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.0663 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.0553 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.0488 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.0431 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.0399 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.0357 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.0313 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.0277 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.0233 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.0221 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.0214 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.0180 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.0155 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0141 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0128 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0134 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0121 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0108 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0101 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0101 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0101 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0102 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0088 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 134.7130 ) * 0.0095 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 0.6970 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 0.6046 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 0.5315 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 0.4306 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 0.3398 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 0.2788 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.2387 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.1814 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.1557 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.1230 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.1009 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.0914 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.0767 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.0638 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.0544 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.0468 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.0425 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.0385 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.0331 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.0278 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.0256 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.0236 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.0217 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0196 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0176 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0165 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0157 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0150 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0144 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0144 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0137 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0130 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0137 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 134.7130 ) * 0.0137
}
set DZResolutionFormula {
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 0.3693 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 0.3135 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 0.3125 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 0.2578 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 0.2221 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 0.1936 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.1686 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.1351 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.1113 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.0983 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.0882 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.0786 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.0684 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.0615 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.0551 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.0516 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.0484 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.0450 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.0416 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.0416 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.0382 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.0350 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.0317 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0316 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0316 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0316 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0348 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0316 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0316 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0316 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0284 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0283 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0315 +\
    ( abs(eta) > 0.0 && abs(eta) <= 0.9 ) * ( pt > 134.7130 ) * 0.0318 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 0.7062 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 0.6010 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 0.5992 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 0.4959 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 0.3877 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 0.3199 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.2649 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.2518 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.1982 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.1587 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.1399 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.1199 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.1031 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.0967 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.0805 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.0736 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.0707 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.0603 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.0609 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.0541 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.0511 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.0443 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.0409 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0408 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0409 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0377 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0375 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0377 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0342 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0342 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0343 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0343 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0309 +\
    ( abs(eta) > 0.9 && abs(eta) <= 1.4 ) * ( pt > 134.7130 ) * 0.0310 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.1823 && pt <= 0.2227 ) * 2.1717 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.2227 && pt <= 0.2720 ) * 2.0715 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.2720 && pt <= 0.3323 ) * 2.0679 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.3323 && pt <= 0.4060 ) * 1.7679 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.4060 && pt <= 0.4959 ) * 1.4393 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.4959 && pt <= 0.6058 ) * 1.1997 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.6058 && pt <= 0.7401 ) * 0.9800 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.7401 && pt <= 0.9041 ) * 0.8251 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 0.9041 && pt <= 1.1044 ) * 0.6695 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.1044 && pt <= 1.3492 ) * 0.5545 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.3492 && pt <= 1.6481 ) * 0.4366 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 1.6481 && pt <= 2.0134 ) * 0.3711 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 2.0134 && pt <= 2.4595 ) * 0.3319 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 2.4595 && pt <= 3.0045 ) * 0.2721 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 3.0045 && pt <= 3.6703 ) * 0.2443 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 3.6703 && pt <= 4.4837 ) * 0.2085 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 4.4837 && pt <= 5.4772 ) * 0.1816 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 5.4772 && pt <= 6.6910 ) * 0.1641 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 6.6910 && pt <= 8.1737 ) * 0.1451 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 8.1737 && pt <= 9.9849 ) * 0.1317 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 9.9849 && pt <= 12.1976 ) * 0.1117 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 12.1976 && pt <= 14.9005 ) * 0.1020 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 14.9005 && pt <= 18.2024 ) * 0.1017 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 18.2024 && pt <= 22.2360 ) * 0.0983 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 22.2360 && pt <= 27.1635 ) * 0.0882 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 27.1635 && pt <= 33.1828 ) * 0.0847 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 33.1828 && pt <= 40.5360 ) * 0.0814 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 40.5360 && pt <= 49.5187 ) * 0.0784 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 49.5187 && pt <= 60.4919 ) * 0.0817 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 60.4919 && pt <= 73.8967 ) * 0.0750 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 73.8967 && pt <= 90.2720 ) * 0.0816 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 90.2720 && pt <= 110.2760 ) * 0.0820 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 110.2760 && pt <= 134.7130 ) * 0.0814 +\
    ( abs(eta) > 1.4 && abs(eta) <= 2.5 ) * ( pt > 134.7130 ) * 0.0850
}
