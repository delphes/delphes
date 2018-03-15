module EnergyScale JetEnergyScale_VLCR05N2 {
    set InputArray  JetMomentumSmearing_VLCR05N2/JER_VLCjetsR05N2
    set OutputArray JES_VLCjetsR05N2

     # Scale Formula
    set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05)}
}
module EnergyScale JetEnergyScale_VLCR05N3 {                                       
 set InputArray  JetMomentumSmearing_VLCR05N3/JER_VLCjetsR05N3                                        
  set OutputArray JES_VLCjetsR05N3                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR05N4 {                                       
 set InputArray  JetMomentumSmearing_VLCR05N4/JER_VLCjetsR05N4                                        
  set OutputArray JES_VLCjetsR05N4                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR05N5 {                                       
 set InputArray  JetMomentumSmearing_VLCR05N5/JER_VLCjetsR05N5                                        
  set OutputArray JES_VLCjetsR05N5                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR05N6 {                                       
 set InputArray  JetMomentumSmearing_VLCR05N6/JER_VLCjetsR05N6                                        
  set OutputArray JES_VLCjetsR05N6                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR05_inclusive {                               
 set InputArray  JetMomentumSmearing_VLCR05_inclusive/JER_VLCjetsR05_inclusive                        
  set OutputArray JES_VLCjetsR05_inclusive                                                   
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR07N2 {                                       
 set InputArray  JetMomentumSmearing_VLCR07N2/JER_VLCjetsR07N2                                        
  set OutputArray JES_VLCjetsR07N2                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR07N3 {                                       
 set InputArray  JetMomentumSmearing_VLCR07N3/JER_VLCjetsR07N3                                        
  set OutputArray JES_VLCjetsR07N3                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR07N4 {                                       
 set InputArray  JetMomentumSmearing_VLCR07N4/JER_VLCjetsR07N4                                        
  set OutputArray JES_VLCjetsR07N4                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR07N5 {                                       
 set InputArray  JetMomentumSmearing_VLCR07N5/JER_VLCjetsR07N5                                        
  set OutputArray JES_VLCjetsR07N5                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR07N6 {                                       
 set InputArray  JetMomentumSmearing_VLCR07N6/JER_VLCjetsR07N6                                        
  set OutputArray JES_VLCjetsR07N6                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR07_inclusive {                               
 set InputArray  JetMomentumSmearing_VLCR07_inclusive/JER_VLCjetsR07_inclusive                        
  set OutputArray JES_VLCjetsR07_inclusive                                                   
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR10N2 {                                       
 set InputArray  JetMomentumSmearing_VLCR10N2/JER_VLCjetsR10N2                                        
  set OutputArray JES_VLCjetsR10N2                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR10N3 {                                       
 set InputArray  JetMomentumSmearing_VLCR10N3/JER_VLCjetsR10N3                                        
  set OutputArray JES_VLCjetsR10N3                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR10N4 {                                       
 set InputArray  JetMomentumSmearing_VLCR10N4/JER_VLCjetsR10N4                                        
  set OutputArray JES_VLCjetsR10N4                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR10N5 {                                       
 set InputArray  JetMomentumSmearing_VLCR10N5/JER_VLCjetsR10N5                                        
  set OutputArray JES_VLCjetsR10N5                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR10N6 {                                       
 set InputArray  JetMomentumSmearing_VLCR10N6/JER_VLCjetsR10N6                                        
  set OutputArray JES_VLCjetsR10N6                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR10_inclusive {                               
 set InputArray  JetMomentumSmearing_VLCR10_inclusive/JER_VLCjetsR10_inclusive                        
  set OutputArray JES_VLCjetsR10_inclusive                                                   
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR12N2 {                                       
 set InputArray  JetMomentumSmearing_VLCR12N2/JER_VLCjetsR12N2                                        
  set OutputArray JES_VLCjetsR12N2                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR12N3 {                                       
 set InputArray  JetMomentumSmearing_VLCR12N3/JER_VLCjetsR12N3                                        
  set OutputArray JES_VLCjetsR12N3                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR12N4 {                                       
 set InputArray  JetMomentumSmearing_VLCR12N4/JER_VLCjetsR12N4                                        
  set OutputArray JES_VLCjetsR12N4                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR12N5 {                                       
 set InputArray  JetMomentumSmearing_VLCR12N5/JER_VLCjetsR12N5                                        
  set OutputArray JES_VLCjetsR12N5                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR12N6 {                                       
 set InputArray  JetMomentumSmearing_VLCR12N6/JER_VLCjetsR12N6                                        
  set OutputArray JES_VLCjetsR12N6                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR12_inclusive {                               
 set InputArray  JetMomentumSmearing_VLCR12_inclusive/JER_VLCjetsR12_inclusive                        
  set OutputArray JES_VLCjetsR12_inclusive                                                   
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR15N2 {                                       
 set InputArray  JetMomentumSmearing_VLCR15N2/JER_VLCjetsR15N2                                        
  set OutputArray JES_VLCjetsR15N2                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR15N3 {                                       
 set InputArray  JetMomentumSmearing_VLCR15N3/JER_VLCjetsR15N3                                        
  set OutputArray JES_VLCjetsR15N3                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR15N4 {                                       
 set InputArray  JetMomentumSmearing_VLCR15N4/JER_VLCjetsR15N4                                        
  set OutputArray JES_VLCjetsR15N4                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR15N5 {                                       
 set InputArray  JetMomentumSmearing_VLCR15N5/JER_VLCjetsR15N5                                        
  set OutputArray JES_VLCjetsR15N5                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR15N6 {                                       
 set InputArray  JetMomentumSmearing_VLCR15N6/JER_VLCjetsR15N6                                        
  set OutputArray JES_VLCjetsR15N6                                                           
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                                                            
module EnergyScale JetEnergyScale_VLCR15_inclusive {                               
 set InputArray  JetMomentumSmearing_VLCR15_inclusive/JER_VLCjetsR15_inclusive                        
  set OutputArray JES_VLCjetsR15_inclusive                                                   
  set ScaleFormula { (abs(eta) < 0.76) * ( 1.0 ) + (abs(eta) >= 0.76 ) * (1.05) }                                                              
}                                                 
