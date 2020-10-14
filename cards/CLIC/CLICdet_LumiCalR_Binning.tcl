# Copy of ILCDelphes implementation by P.Sopicki

set pi [expr {acos(-1)} ]

# LumiCal eta range 3.0 - 4.0 (no beam crossing boost)
# Rear part
set PhiBins {}
  for {set i -24} {$i <= 24} {incr i} {
    add PhiBins [expr {$i * $pi/24.0} ]
  }

for {set i 0} {$i <= 64} {incr i} {
    set eta [expr {-4.0 + $i * 1.0/64.0} ]
    add EtaPhiBins $eta $PhiBins
  }
