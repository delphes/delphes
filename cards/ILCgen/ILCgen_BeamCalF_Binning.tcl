#P.Sopicki: based on plots from D.Jeans:
# corrected by A.F.Zarnecki

set pi [expr {acos(-1)} ]

# BeamCal eta range 4.0 - 5.8 -> 2.099-0.347 =1.752deg
#
# Inner part (4.6-5.8):  1.152-0.347 = 0.805 deg
# 360/0.0973/cosh(5.3) = ~36 => 32 bins (taking key hole into account!!!)
#
# Outer part (4.0-4.6):  2.099-1.152 = 0.947 deg
# 360/0.0973/cosh(4.4) = ~90 => 90 bins
#Front part

  set PhiBins {}
  for {set i -45} {$i <= 45} {incr i} {
    add PhiBins [expr {$i * $pi/45.0}]
  }
  for {set i 0} {$i <= 12} {incr i} {
    set eta [expr {4.0 + $i * 0.5999/12.0}]
    add EtaPhiBins $eta $PhiBins
  }
  set PhiBins {}
  for {set i -16} {$i <= 16} {incr i} {
    add PhiBins [expr {$i * $pi/18.0}]
  }
  for {set i 0} {$i <= 8} {incr i} {
    set eta [expr {4.6 + $i * 1.2/8.0}]
    add EtaPhiBins $eta $PhiBins
  }
