# A.F.Zarnecki, based on ILCDelphes implementation

set pi [expr {acos(-1)} ]

# BeamCal eta range 4.0 - 5.3 -> 2.099-0.572 =1.527deg
#
# Inner part (4.5-5.3):  1.273-0.347 = 0.926 deg
# 360/0.0973/cosh(5.3) = ~36 => 32 bins (taking key hole into account!!!)
#
# Outer part (4.0-4.5):  2.099-1.273 = 0.826 deg
# 360/0.0973/cosh(4.4) = ~90 => 90 bins
# Rear part

  set PhiBins {}
  for {set i -45} {$i <= 45} {incr i} {
    add PhiBins [expr {$i * $pi/45.0}]
  }
  for {set i 0} {$i <= 10} {incr i} {
    set eta [expr {-4.0 - $i * 0.4999/10.0}]
    add EtaPhiBins $eta $PhiBins
  }
  set PhiBins {}
  for {set i -16} {$i <= 16} {incr i} {
    add PhiBins [expr {$i * $pi/18.0}]
  }
  for {set i 0} {$i <= 10} {incr i} {
    set eta [expr {-4.5 - $i * 0.8/10.0}]
    add EtaPhiBins $eta $PhiBins
  }
