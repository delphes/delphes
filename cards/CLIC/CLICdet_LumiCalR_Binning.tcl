# Copy of ILCDelphes implementation by P.Sopicki


# LumiCal eta range 3.0 - 4.0 (no beam crossing boost)
# Rear part
set PhiBins 48

for {set i 0} {$i <= 64} {incr i} {
    set eta [expr {-4.0 + $i * 1.0/64.0} ]
    add EtaPhiBins $eta $PhiBins
  }
