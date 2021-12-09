#P.Sopicki: based on plots from D.Jeans:
# corrected by A.F.Zarnecki

set pi [expr {acos(-1)} ]

# ECAL barrel
# 0.003 segm in eta and phi(rad) - towers (5x5 mm^2) for |eta| <= 1.1

# to get dphi = 0.003rad/0.17deg -> 360/0.17 /cosh(0.55) =~ 1813 -> 1800 
set PhiBins {}
  for {set i -900} {$i <= 900} {incr i} {
    add PhiBins [expr {$i * $pi/900} ]
  }

#
# |eta| = 1.1, delta = 0.003 -> 733,(3) =~ 720 
  for {set i 1} {$i <= 720} {incr i} {
    set eta [expr {-1.1 + $i * 2.2/720} ]
    add EtaPhiBins $eta $PhiBins
  }

# ECAL endcaps 
# eta 1.1 - 2   : dphi(rad) and deta = 0.006
# eta 2 - 2.5   : d = 0.011
# eta 2.5 - 3   : d = 0.016

#360/0.34/cosh(1.65) : ~387.9 -> 380 divisions
#0.9/0.006: 150
set PhiBins {}
  for {set i -190} {$i <= 190} {incr i} {
    add PhiBins [expr {$i * $pi/190} ]
  }

for {set i 1} {$i <= 150} {incr i} {
    set eta [expr {-2.0 + $i * 0.9/150.0} ]
    add EtaPhiBins $eta $PhiBins
  }

for {set i 1} {$i <= 150} {incr i} {
    set eta [expr {1.1 + $i * 0.9/150.0} ]
    add EtaPhiBins $eta $PhiBins
  }

#360/0.63/cosh(2.25) : ~119 -> 120 divisions
#0.5/0.011: ~45.5 -> 48 
set PhiBins {}
  for {set i -60} {$i <= 60} {incr i} {
    add PhiBins [expr {$i * $pi/60} ]
  }

for {set i 1} {$i <= 48} {incr i} {
    set eta [expr {-2.5 + $i * 0.5/48.0} ]
    add EtaPhiBins $eta $PhiBins
  }

for {set i 1} {$i <= 48} {incr i} {
    set eta [expr {2.0 + $i * 0.5/48.0} ]
    add EtaPhiBins $eta $PhiBins
  }

#360/0.63/cosh(2.75) : ~50 divisions
#0.5/0.016: ~31 -> 30
set PhiBins {}
  for {set i -25} {$i <= 25} {incr i} {
    add PhiBins [expr {$i * $pi/25} ]
  }

for {set i 0} {$i <= 30} {incr i} {
    set eta [expr {-3.0 + $i * 0.5/30.0} ]
    add EtaPhiBins $eta $PhiBins
  }

for {set i 1} {$i <= 30} {incr i} {
    set eta [expr {2.5 + $i * 0.5/30.0} ]
    add EtaPhiBins $eta $PhiBins
  }
