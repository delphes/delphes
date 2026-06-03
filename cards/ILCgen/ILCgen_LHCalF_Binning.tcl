#P.Sopicki: based on plots from D.Jeans:


# LHCAL eta range 2.8 - 3.8
# Front part
set PhiBins 24

for {set i 0} {$i <= 32} {incr i} {
    set eta [expr {2.8 + $i * 1.0/32.0} ]
    add EtaPhiBins $eta $PhiBins
  }
