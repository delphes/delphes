#P.Sopicki: based on plots from D.Jeans:
# corrected by A.F.Zarnecki

set pi [expr {acos(-1)} ]

# HCAL barrel
# 0.015 segm in eta and phi(rad) - towers (5x5 mm^2) for |eta| <= 1.1

# to get dphi = 0.015rad/0.86deg -> 360/0.86 /cosh(0.55) = ~362.7 -> 360
# to get deta = 0.015 2.2/0.015=146.(6)
set PhiBins {}
  for {set i -180} {$i <= 180} {incr i} {
    add PhiBins [expr {$i * $pi/180.0} ]
  }

for {set i 1} {$i <= 147} {incr i} {
    set eta [expr {-1.1 + $i * 2.2/147.0} ]
    add EtaPhiBins $eta $PhiBins
  }

# HCAL endcaps 
# eta 1.1 - 2   : dphi(rad) and deta = 0.03
# eta 2 - 2.5   : d = 0.065
# eta 2.5 - 2.8   : d = 0.08

#360/1.72 deg/cosh(1.65) : ~77.6 -> 76 divisions
#0.9/0.03: 30
set PhiBins {}
  for {set i -38} {$i <= 38} {incr i} {
    add PhiBins [expr {$i * $pi/38} ]
  }

for {set i 1} {$i <= 30} {incr i} {
    set eta [expr {-2.0 + $i * 0.9/30.0} ]
    add EtaPhiBins $eta $PhiBins
  }

for {set i 1} {$i <= 30} {incr i} {
    set eta [expr {1.1 + $i * 0.9/30.0} ]
    add EtaPhiBins $eta $PhiBins
  }

#360/3.72 deg/cosh(2.25) : ~20  divisions
#0.5/0.065: ~7.65 -> 8  
set PhiBins {}
  for {set i -20} {$i <= 20} {incr i} {
    add PhiBins [expr {$i * $pi/20} ]
  }

for {set i 1} {$i <= 8} {incr i} {
    set eta [expr {-2.5 + $i * 0.5/8.0} ]
    add EtaPhiBins $eta $PhiBins
  }

for {set i 1} {$i <= 8} {incr i} {
    set eta [expr {2.0 + $i * 0.5/8.0} ]
    add EtaPhiBins $eta $PhiBins
}

#360/4.58 deg/cosh(2.65) : ~11 -> 10 divisions
#0.3/0.08: ~4
set PhiBins {}
  for {set i -5} {$i <= 5} {incr i} {
    add PhiBins [expr {$i * $pi/5} ]
  }

for {set i 0} {$i <= 4} {incr i} {
    set eta [expr {-2.8 + $i * 0.3/4.0} ]
    add EtaPhiBins $eta $PhiBins
  }

for {set i 1} {$i <= 4} {incr i} {
    set eta [expr {2.5 + $i * 0.3/4.0} ]
    add EtaPhiBins $eta $PhiBins
  }
