#P.Sopicki: based on plots from D.Jeans:

set pi [expr {acos(-1)} ]

# LHCAL eta range 2.8 - 3.8
# Front part
set PhiBins {}
  for {set i -12} {$i <= 12} {incr i} {
    add PhiBins [expr {$i * $pi/12.0} ]
  }

for {set i 0} {$i <= 32} {incr i} {
    set eta [expr {2.8 + $i * 1.0/32.0} ]
    add EtaPhiBins $eta $PhiBins
  }
