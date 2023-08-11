######################################################################################################################
# 08/09/2023 Check-Ping Wong (cwong1@bnl.gov) and Jihee Kim (jkim11@bnl.gov)
# Added ePIC ECal and HCal resolutions
# Changed magnetic field stregth from 3 T to 1.7 T (just scaling down); set Bz 
# 08/10/2023 Jihee Kim
# Added ePIC Tracking resolution and DCA transverse pointing resolution
# 08/11/2023 Jihee Kim
# Added ePIC Tracking eifficiency and DCA longitudinal pointing resolution
######################################################################################################################
# ATHENA detector model. Based on parametrizations from G4 simulations made by the ATHENA collaboration
# email: miguel.arratia@ucr.edu, ssekula@mail.smu.edu
#######################################################################################################################

#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  BeamSpotFilter
  ParticlePropagator

  ChargedHadronTrackingEfficiency
  ElectronTrackingEfficiency
  MuonTrackingEfficiency

  ChargedHadronSmearing
  ElectronSmearing
  MuonSmearing

  TrackMerger

  ECal
  HCal

  Calorimeter
  EFlowMerger
  EFlowFilter

  PhotonEfficiency
  PhotonIsolation

  ElectronFilter
  ElectronEfficiency
  ElectronIsolation

  ChargedHadronFilter
  MissingET

  NeutrinoFilter
  GenJetFinder
  GenMissingET

  FastJetFinder

  JetEnergyScale

  JetFlavorAssociation
  GenJetFlavorAssociation

  UniqueObjectFinder

  ScalarHT

  TrackCountingBTagging

  pfRICH
  barrelDIRC_epid
  barrelDIRC_hpid  
  BTOF_e
  BTOF_h
  dualRICH_aerogel
  dualRICH_c2f6

  TreeWriter
}

#################################
# Propagate particles in cylinder
#################################

#####################################################
# GenBeamSpotFilter
# Saves a particle intended to represent the beamspot
#####################################################

module BeamSpotFilter BeamSpotFilter {
    set InputArray Delphes/stableParticles
    set OutputArray beamSpotParticle
}

module ParticlePropagator ParticlePropagator {
    set InputArray Delphes/stableParticles
    set OutputArray stableParticles
    set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
    set ChargedHadronOutputArray chargedHadrons
    set ElectronOutputArray electrons
    set MuonOutputArray muons
    # radius of the magnetic field coverage, in m
    # ePIC outer radius of magnet
    set Radius 1.77
    # half-length of the magnetic field coverage, in m
    set HalfLength 1.92
    # magnetic field, in T
    set Bz 1.7
}

####################################
# Common Tracking Efficiency Model
####################################
# Refer to ePIC tracking efficiency
# According to Stephen Maple
# In CraterLake version of ePIC tracking
# Individual tracking layers in simulation are 100% efficient
# When we use truth seeding 1 charged particle = 1 track (assuming it produces sufficient hits for track fitting)
# See slide 5 of Barak's presentation https://indico.bnl.gov/event/20126/contributions/78820/attachments/48724/83521/real_seeding_072023.pdf
set CommonTrackingEfficiency {
    (eta > -3.5 && eta <= -3.0 ) * (pt*cosh(eta) > 1.25 && pt*cosh(eta)<6.0)   * (1.0) +
    (eta > -3.5 && eta <= -3.0 ) * (pt*cosh(eta) > 6.0 )                       * (1.0) +
    (eta > -3.0 && eta <= -2.5 ) * (pt*cosh(eta) > 0.55 && pt*cosh(eta)<2.0)   * (1.0) +
    (eta > -3.0 && eta <= -2.5 ) * (pt*cosh(eta) > 2.0)                        * (1.0) +
    (eta > -2.5 && eta <= -2.0)  * (pt*cosh(eta)> 0.45 && pt*cosh(eta)<0.6)    * (1.0) +
    (eta > -2.5 && eta <= -2.0)  * (pt*cosh(eta)>0.6)                          * (1.0) +
    (eta > -2.0 && eta <= -1.5)  * (pt*cosh(eta)> 0.250 && pt*cosh(eta)<0.500)* (1.0) +
    (eta > -2.0 && eta <= -1.5)  * (pt*cosh(eta)>0.500)                        * (1.0) +    
    (eta > -1.5 && eta <= -1.0)  * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.300)* (1.0) +
    (eta > -1.5 && eta <= -1.0)  * (pt*cosh(eta) > 0.300)                      * (1.0) +
    (eta > -1.0 && eta <= -0.5)  * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.200)* (1.0) +
    (eta > -1.0 && eta <= -0.5)  * (pt*cosh(eta) > 0.200)                      * (1.0) +
    (eta > -0.5 && eta <= 0.0)   * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.200)*(1.0) +
    (eta > -0.5 && eta <= 0.0)   * (pt*cosh(eta) > 0.200 )                     * (1.0) +
    (eta > 0.0 && eta <= 0.5 )   * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.200)* (1.0) +
    (eta > 0.0 && eta <= 0.5 )   * (pt*cosh(eta) > 0.200)                      * (1.0) +
    (eta > 0.5 && eta <= 1.0 )   * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.200)* (1.0) +
    (eta > 0.5 && eta <= 1.0 )   * (pt*cosh(eta) > 0.200)                      * (1.0) +
    (eta > 1.0 && eta <= 1.5)    * (pt*cosh(eta) > 0.150 && pt*cosh(eta)<0.200)* (1.0) +
    (eta > 1.0 && eta <= 1.5)    * (pt*cosh(eta) > 0.200)                      * (1.0) +
    (eta > 1.5 && eta <= 2.0)    * (pt*cosh(eta) > 0.250 && pt*cosh(eta)<0.500)* (1.0) +
    (eta > 1.5 && eta <= 2.0)    * (pt*cosh(eta) > 0.500)                      * (1.0) +
    (eta > 2.0 && eta <= 2.5)    * (pt*cosh(eta) > 0.350 && pt*cosh(eta)<0.700)* (1.0) +
    (eta > 2.0 && eta <= 2.5)    * (pt*cosh(eta) > 0.700)                      * (1.0) +
    (eta > 2.5 && eta <= 3.0)    * (pt*cosh(eta) > 0.550 && pt*cosh(eta)<2.0)  * (1.0) +
    (eta > 2.5 && eta <= 3.0)    * (pt*cosh(eta) > 2.0)                        * (1.0) +
    (eta > 3.0 && eta <= 3.5)    * (pt*cosh(eta) > 0.850 && pt*cosh(eta)<4.0)  * (1.0) +
    (eta > 3.0 && eta <= 3.5)    * (pt*cosh(eta) > 4.0)                        * (1.0) +
    (abs(eta) > 3.5)  * (0.00)+
    0.0
}

# Refer to ePIC tracking resolution
# Provided by Stephen Maple
# Used CraterLake version of the ePIC tracking system with truth seeding
set CommonTrackingResolution {
    (eta<=-3.0 && eta>-3.5)  * (sqrt( (7.4e-2)^2 + (pt*cosh(eta)*4.92e-3)^2  ) )  +
    (eta<=-2.5 && eta>-3.0)  * (sqrt( (4.0e-2)^2 + (pt*cosh(eta)*1.55e-3)^2  ) )  +
    (eta<=-2.0 && eta>-2.5)  * (sqrt( (2.3e-2)^2 + (pt*cosh(eta)*0.0)^2  ) )  +
    (eta<=-1.5 && eta>-2.0)  * (sqrt( (1.4e-2)^2 + (pt*cosh(eta)*0.0)^2  ) )  +
    (eta<=-1.0 && eta>-1.5)  * (sqrt( (1.2e-2)^2 + (pt*cosh(eta)*0.0)^2  ) )  +
    (eta<=-0.5 && eta>-1.0)  * (sqrt( (0.5e-2)^2 + (pt*cosh(eta)*4.8e-4)^2  ) )  +
    (eta<= 0.0 && eta>-0.5)  * (sqrt( (0.4e-2)^2 + (pt*cosh(eta)*5.3e-4)^2  ) )  +

    (eta<=0.5 && eta>0)  * (sqrt( (0.4e-2)^2 + (pt*cosh(eta)*5.3e-4)^2  ) )  +
    (eta<=1.0 && eta>0.5) * (sqrt( (0.5e-2)^2 + (pt*cosh(eta)*5.1e-4)^2   ) )  +
    (eta<=1.5 && eta>1.0) * (sqrt( (1.2e-2)^2 + (pt*cosh(eta)*7.0e-5)^2   ) )  +
    (eta<=2.0 && eta>1.5) * (sqrt( (1.3e-2)^2 + (pt*cosh(eta)*0.0)^2   ) )  +
    (eta<=2.5 && eta>2.0) * (sqrt( (1.8e-2)^2 + (pt*cosh(eta)*0.0)^2   ) )  +
    (eta<=3.0 && eta>2.5) * (sqrt( (3.0e-2)^2 + (pt*cosh(eta)*3.7e-4)^2   ) )  +
    (eta<=3.5 && eta>3.0) * (sqrt( (5.4e-2)^2 + (pt*cosh(eta)*2.31e-3)^2  ) )  
}

####################################
# Charged hadron tracking efficiency
####################################

module Efficiency ChargedHadronTrackingEfficiency {
  set InputArray ParticlePropagator/chargedHadrons
  set OutputArray chargedHadrons
  set EfficiencyFormula $CommonTrackingEfficiency
}

##############################
# Electron tracking efficiency
##############################

module Efficiency ElectronTrackingEfficiency {
  set InputArray ParticlePropagator/electrons
  set OutputArray electrons
  set EfficiencyFormula $CommonTrackingEfficiency
}

##############################
# Muon tracking efficiency
##############################
module Efficiency MuonTrackingEfficiency {
  set InputArray ParticlePropagator/muons
  set OutputArray muons
  set EfficiencyFormula $CommonTrackingEfficiency

}

########################################
# Smearing for charged hadrons
########################################

module TrackSmearing ChargedHadronSmearing {
  set InputArray ChargedHadronTrackingEfficiency/chargedHadrons
  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
  set OutputArray chargedHadrons
#  set ApplyToPileUp true
  # magnetic field, in T
  set Bz 1.7
  set PResolutionFormula $CommonTrackingResolution
  set CtgThetaResolutionFormula { 0.0 }
  set PhiResolutionFormula { 0.0 }

    
  # Refer to ePIC tracking resolution
  # Provided by Stephen Maple
  # Used CraterLake version of the ePIC tracking system with truth seeding
  # DCA resolution in transverse (D0) [mm]
  set D0ResolutionFormula "
    (eta<=0.0 && eta>-0.5)    * (sqrt( (0.0066)^2 +   (0.0247/pt)^2   ) )  +
    (eta<=-0.5 && eta>-1.0)   * (sqrt( (0.0067)^2 +   (0.0279/pt)^2   ) )  +
    (eta<=-1.0 && eta>-1.5)   * (sqrt( (0.0093)^2 +   (0.0438/pt)^2   ) )  +
    (eta<=-1.5 && eta>-2.0)   * (sqrt( (0.0126)^2 +   (0.0583/pt)^2   ) )  +
    (eta<=-2.0 && eta>-2.5)   * (sqrt( (0.0156)^2 +   (0.0678/pt)^2   ) )  +
    (eta<=-2.5 && eta>-3.0)   * (sqrt( (0.0398)^2 +   (0.0711/pt)^2   ) )  +
    (eta<=-3.0 && eta>-3.5)   * (sqrt( (0.1032)^2 +   (0.0535/pt)^2   ) )  +

    (eta<=0.5 && eta>0)     * (sqrt( (0.0065)^2 +   (0.0241/pt)^2   ) )  +                                                                                                   
    (eta<=1.0 && eta>0.5)   * (sqrt( (0.0068)^2 +   (0.0283/pt)^2   ) )  +                                                                                                              
    (eta<=1.5 && eta>1.0)   * (sqrt( (0.0084)^2 +   (0.0433/pt)^2   ) )  +                                                                                                             
    (eta<=2.0 && eta>1.5)   * (sqrt( (0.0133)^2 +   (0.0579/pt)^2   ) )  +                                                                                                                 
    (eta<=2.5 && eta>2.0)   * (sqrt( (0.0141)^2 +   (0.0659/pt)^2   ) )  +                                                                                                              
    (eta<=3.0 && eta>2.5)   * (sqrt( (0.0308)^2 +   (0.0685/pt)^2   ) )  +                                                                                                               
    (eta<=3.5 && eta>3.0)   * (sqrt( (0.0770)^2 +   (0.0575/pt)^2   ) )  
  "
  
  # Refer to ePIC tracking resolution
  # Provided by Stephen Maple
  # Used CraterLake version of the ePIC tracking system with truth seeding
  # DCA resolution in longitudinal (DZ) [mm]
  set DZResolutionFormula "
    (eta <= -2.5 && eta > -3.0) * (sqrt( (0.1322)^2 + (0.4925/pt)^2)) + 
    (eta <= -2.0 && eta > -2.5) * (sqrt( (0.0828)^2 + (0.2337/pt)^2)) + 
    (eta <= -1.5 && eta > -2.0) * (sqrt( (0.0380)^2 + (0.1099/pt)^2)) + 
    (eta <= -1.0 && eta > -1.5) * (sqrt( (0.0160)^2 + (0.0666/pt)^2)) + 
    (eta <= -0.5 && eta > -1.0) * (sqrt( (0.0055)^2 + (0.0351/pt)^2)) + 
    (eta <=  0.0 && eta > -0.5) * (sqrt( (0.0047)^2 + (0.0260/pt)^2)) + 
    (eta <=  0.5 && eta >  0.0) * (sqrt( (0.0047)^2 + (0.0264/pt)^2)) + 
    (eta <=  1.0 && eta >  0.5) * (sqrt( (0.0057)^2 + (0.0351/pt)^2)) + 
    (eta <=  1.5 && eta >  1.0) * (sqrt( (0.0160)^2 + (0.0658/pt)^2)) + 
    (eta <=  2.0 && eta >  1.5) * (sqrt( (0.0382)^2 + (0.1122/pt)^2)) + 
    (eta <=  2.5 && eta >  2.0) * (sqrt( (0.0864)^2 + (0.2372/pt)^2)) + 
    (eta <=  3.0 && eta >  2.5) * (sqrt( (0.1433)^2 + (0.4948/pt)^2))
  "

}

###################################
# Smearing for muons
###################################

module TrackSmearing MuonSmearing {
  set InputArray MuonTrackingEfficiency/muons
  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
  set OutputArray muons
#  set ApplyToPileUp true
  # magnetic field in T
  set Bz 1.7
  set PResolutionFormula $CommonTrackingResolution
  set CtgThetaResolutionFormula { 0.0 }
  set PhiResolutionFormula { 0.0 }

  # Refer to ePIC tracking resolution
  # Provided by Stephen Maple
  # Used CraterLake version of the ePIC tracking system with truth seeding
  # DCA resolution in transverse (D0) [mm]
  set D0ResolutionFormula "
    (eta<=0.0 && eta>-0.5)    * (sqrt( (0.0066)^2 +   (0.0247/pt)^2   ) )  +
    (eta<=-0.5 && eta>-1.0)   * (sqrt( (0.0067)^2 +   (0.0279/pt)^2   ) )  +
    (eta<=-1.0 && eta>-1.5)   * (sqrt( (0.0093)^2 +   (0.0438/pt)^2   ) )  +
    (eta<=-1.5 && eta>-2.0)   * (sqrt( (0.0126)^2 +   (0.0583/pt)^2   ) )  +
    (eta<=-2.0 && eta>-2.5)   * (sqrt( (0.0156)^2 +   (0.0678/pt)^2   ) )  +
    (eta<=-2.5 && eta>-3.0)   * (sqrt( (0.0398)^2 +   (0.0711/pt)^2   ) )  +
    (eta<=-3.0 && eta>-3.5)   * (sqrt( (0.1032)^2 +   (0.0535/pt)^2   ) )  +

    (eta<=0.5 && eta>0)     * (sqrt( (0.0065)^2 +   (0.0241/pt)^2   ) )  +                                                                                                   
    (eta<=1.0 && eta>0.5)   * (sqrt( (0.0068)^2 +   (0.0283/pt)^2   ) )  +                                                                                                              
    (eta<=1.5 && eta>1.0)   * (sqrt( (0.0084)^2 +   (0.0433/pt)^2   ) )  +                                                                                                             
   
    (eta<=2.0 && eta>1.5)   * (sqrt( (0.0133)^2 +   (0.0579/pt)^2   ) )  +                                                                                                                 
    (eta<=2.5 && eta>2.0)   * (sqrt( (0.0141)^2 +   (0.0659/pt)^2   ) )  +                                                                                                              
    (eta<=3.0 && eta>2.5)   * (sqrt( (0.0308)^2 +   (0.0685/pt)^2   ) )  +                                                                                                               
    (eta<=3.5 && eta>3.0)   * (sqrt( (0.0770)^2 +   (0.0575/pt)^2   ) )  
  "
  # Refer to ePIC tracking resolution
  # Provided by Stephen Maple
  # Used CraterLake version of the ePIC tracking system with truth seeding
  # DCA resolution in longitudinal (DZ) [mm]
  set DZResolutionFormula "
    (eta <= -2.5 && eta > -3.0) * (sqrt( (0.1322)^2 + (0.4925/pt)^2)) + 
    (eta <= -2.0 && eta > -2.5) * (sqrt( (0.0828)^2 + (0.2337/pt)^2)) +    
    (eta <= -1.5 && eta > -2.0) * (sqrt( (0.0380)^2 + (0.1099/pt)^2)) +    
    (eta <= -1.0 && eta > -1.5) * (sqrt( (0.0160)^2 + (0.0666/pt)^2)) +    
    (eta <= -0.5 && eta > -1.0) * (sqrt( (0.0055)^2 + (0.0351/pt)^2)) +    
    (eta <=  0.0 && eta > -0.5) * (sqrt( (0.0047)^2 + (0.0260/pt)^2)) +    
    (eta <=  0.5 && eta >  0.0) * (sqrt( (0.0047)^2 + (0.0264/pt)^2)) +    
    (eta <=  1.0 && eta >  0.5) * (sqrt( (0.0057)^2 + (0.0351/pt)^2)) +    
    (eta <=  1.5 && eta >  1.0) * (sqrt( (0.0160)^2 + (0.0658/pt)^2)) +    
    (eta <=  2.0 && eta >  1.5) * (sqrt( (0.0382)^2 + (0.1122/pt)^2)) +    
    (eta <=  2.5 && eta >  2.0) * (sqrt( (0.0864)^2 + (0.2372/pt)^2)) +    
    (eta <=  3.0 && eta >  2.5) * (sqrt( (0.1433)^2 + (0.4948/pt)^2))    
  "

}

###################################
# Smearing for electrons
###################################


module TrackSmearing ElectronSmearing {
  set InputArray ElectronTrackingEfficiency/electrons
  set BeamSpotInputArray BeamSpotFilter/beamSpotParticle
  set OutputArray electrons
#  set ApplyToPileUp true
  # magnetic field, in T
  set Bz 1.7
  set PResolutionFormula $CommonTrackingResolution
  set CtgThetaResolutionFormula { 0.0 }
  set PhiResolutionFormula { 0.0 }

  # Refer to ePIC tracking resolution
  # Provided by Stephen Maple
  # Used CraterLake version of the ePIC tracking system with truth seeding
  # DCA resolution in transverse (D0) [mm]
  set D0ResolutionFormula "
    (eta<=0.0 && eta>-0.5)    * (sqrt( (0.0066)^2 +   (0.0247/pt)^2   ) )  +
    (eta<=-0.5 && eta>-1.0)   * (sqrt( (0.0067)^2 +   (0.0279/pt)^2   ) )  +
    (eta<=-1.0 && eta>-1.5)   * (sqrt( (0.0093)^2 +   (0.0438/pt)^2   ) )  +
    (eta<=-1.5 && eta>-2.0)   * (sqrt( (0.0126)^2 +   (0.0583/pt)^2   ) )  +
    (eta<=-2.0 && eta>-2.5)   * (sqrt( (0.0156)^2 +   (0.0678/pt)^2   ) )  +
    (eta<=-2.5 && eta>-3.0)   * (sqrt( (0.0398)^2 +   (0.0711/pt)^2   ) )  +
    (eta<=-3.0 && eta>-3.5)   * (sqrt( (0.1032)^2 +   (0.0535/pt)^2   ) )  +

    (eta<=0.5 && eta>0)     * (sqrt( (0.0065)^2 +   (0.0241/pt)^2   ) )  +                                                                                                   
    (eta<=1.0 && eta>0.5)   * (sqrt( (0.0068)^2 +   (0.0283/pt)^2   ) )  +                                                                                                              
    (eta<=1.5 && eta>1.0)   * (sqrt( (0.0084)^2 +   (0.0433/pt)^2   ) )  +                                                                                                             
   
    (eta<=2.0 && eta>1.5)   * (sqrt( (0.0133)^2 +   (0.0579/pt)^2   ) )  +                                                                                                                 
    (eta<=2.5 && eta>2.0)   * (sqrt( (0.0141)^2 +   (0.0659/pt)^2   ) )  +                                                                                                              
    (eta<=3.0 && eta>2.5)   * (sqrt( (0.0308)^2 +   (0.0685/pt)^2   ) )  +                                                                                                               
    (eta<=3.5 && eta>3.0)   * (sqrt( (0.0770)^2 +   (0.0575/pt)^2   ) )  
  "

  # Refer to ePIC tracking resolution
  # Provided by Stephen Maple
  # Used CraterLake version of the ePIC tracking system with truth seeding
  # DCA resolution in longitudinal (DZ) [mm]
  set DZResolutionFormula "
    (eta <= -2.5 && eta > -3.0) * (sqrt( (0.1322)^2 + (0.4925/pt)^2)) +    
    (eta <= -2.0 && eta > -2.5) * (sqrt( (0.0828)^2 + (0.2337/pt)^2)) +    
    (eta <= -1.5 && eta > -2.0) * (sqrt( (0.0380)^2 + (0.1099/pt)^2)) +    
    (eta <= -1.0 && eta > -1.5) * (sqrt( (0.0160)^2 + (0.0666/pt)^2)) +    
    (eta <= -0.5 && eta > -1.0) * (sqrt( (0.0055)^2 + (0.0351/pt)^2)) +    
    (eta <=  0.0 && eta > -0.5) * (sqrt( (0.0047)^2 + (0.0260/pt)^2)) +    
    (eta <=  0.5 && eta >  0.0) * (sqrt( (0.0047)^2 + (0.0264/pt)^2)) +    
    (eta <=  1.0 && eta >  0.5) * (sqrt( (0.0057)^2 + (0.0351/pt)^2)) +    
    (eta <=  1.5 && eta >  1.0) * (sqrt( (0.0160)^2 + (0.0658/pt)^2)) +    
    (eta <=  2.0 && eta >  1.5) * (sqrt( (0.0382)^2 + (0.1122/pt)^2)) +    
    (eta <=  2.5 && eta >  2.0) * (sqrt( (0.0864)^2 + (0.2372/pt)^2)) +    
    (eta <=  3.0 && eta >  2.5) * (sqrt( (0.1433)^2 + (0.4948/pt)^2))  
  "
}

##############
# Track merger
##############

module Merger TrackMerger {
# add InputArray InputArray
  add InputArray ChargedHadronSmearing/chargedHadrons
  add InputArray ElectronSmearing/electrons
  add InputArray MuonSmearing/muons
  set OutputArray tracks
}

#############
#   ECAL
#############

module SimpleCalorimeter ECal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray TrackMerger/tracks

  set TowerOutputArray ecalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowPhotons

  set IsEcal true
  set EnergyMin 0.050
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # ATHENA model
  # lists of the edges of each tower in eta and phi
  # each list starts with the lower edge of the first tower
  # the list ends with the higher edged of the last tower

    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    for {set i -10} {$i <=10} {incr i} {
	set eta [expr {$i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }

    
   ## assume 0.1 x 0.1 (real cell size will be smaller, so this is to represent some cluster)
    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }

    foreach eta {-4.0 -3.5 -3.41102102 -3.26996837 -3.14642305 -3.03653567 -2.93760447 -2.84766006 -2.76522251 -2.68915144 -2.61854952 -2.55269788 -2.49101173 -2.43300894 -2.3782873  -2.3265078  -2.27738197 -2.23066235 -2.18613503 -2.14361383 -2.10293569 -2.063957 -2.02655061 -1.99060337 -1.95601417 -1.92269228 -1.89055593 -1.8595312  -1.82955102 -1.80055436 -1.77248548 -1.74529337 -1.71893119 -1.69335587 -1.66852765 -1.64440978 -1.62096821 -1.59817135 -1.57598979 -1.55439612 -1.53336478 -1.51287184 -1.4928949  -1.47341295 -1.45440623 -1.43585618 -1.41774529 -1.40005705 -1.38277588 -1.36588703 -1.34937654 -1.33323117 -1.31743839 -1.30198626 -1.28686345 -1.27205918 -1.25756317 -1.24336562 -1.22945719 -1.21582897 -1.20247241 -1.18937936 -1.17654201 -1.16395288 -1.15160481 -1.13949092 -1.12760462 -1.11593955 -1.10448965 -1.09324904 -1.08221211	-1.07137342 -1.06072776 -1.0502701  -1.03999558} {
	add EtaPhiBins $eta $PhiBins
    }
    foreach eta {1.0 1.0502701  1.06072776 1.07137342 1.08221211 1.09324904 1.10448965 1.11593955 1.12760462 1.13949092 1.15160481 1.16395288 1.17654201 1.18937936 1.20247241 1.21582897 1.22945719 1.24336562 1.25756317 1.27205918 1.28686345 1.30198626 1.31743839 1.33323117 1.34937654 1.36588703 1.38277588 1.40005705 1.41774529 1.43585618 1.45440623 1.47341295 1.4928949  1.51287184 1.53336478 1.55439612 1.57598979 1.59817135 1.62096821 1.64440978 1.66852765 1.69335587 1.71893119 1.74529337 1.77248548 1.80055436 1.82955102 1.8595312 1.89055593 1.92269228 1.95601417 1.99060337 2.02655061 2.063957 2.10293569 2.14361383 2.18613503 2.23066235 2.27738197 2.3265078 2.3782873  2.43300894 2.49101173 2.55269788 2.61854952 2.68915144 2.76522251 2.84766006 2.93760447 3.03653567 3.14642305 3.26996837 3.41102102 3.5 4.0} {
   	add EtaPhiBins $eta $PhiBins
    }
  # ATHENA model 
  add EnergyFraction {0} {0.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {1.0}
  add EnergyFraction {22} {1.0}
  add EnergyFraction {111} {1.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.3}
  add EnergyFraction {3122} {0.3}

  # Refer to ePIC ECal resolution (presented at Hard Probes 2023)
  # https://wiki.bnl.gov/EPIC/index.php?title=File:ECalResolution_electron_Fitted.pdf
  set ResolutionFormula {          (eta <= -1.5 && eta > -3.5)                         * sqrt(energy^2*0.012^2 + energy*0.019^2 )+ \
				   (eta <=  1.5 && eta > -1.5)                        * sqrt(energy^2*0.001^2 + energy*0.049^2 )+ \
				   (eta <=  3.5 && eta >  1.5)                         * sqrt(energy^2*0.032^2 + energy*0.106^2 )}

}


#############
#   HCAL
#############

module SimpleCalorimeter HCal {
  set ParticleInputArray ParticlePropagator/stableParticles
  set TrackInputArray ECal/eflowTracks

  set TowerOutputArray hcalTowers
  set EFlowTrackOutputArray eflowTracks
  set EFlowTowerOutputArray eflowNeutralHadrons

  set IsEcal false

  set EnergyMin 0.300
  set EnergySignificanceMin 1.0

  set SmearTowerCenter true

  set pi [expr {acos(-1)}]

  # ATHENA model
    set PhiBins {}
    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    for {set i -10} {$i <=10} {incr i} {
	set eta [expr {$i * 0.1}]
	add EtaPhiBins $eta $PhiBins
    }

    for {set i -30} {$i <=30} {incr i} {
	add PhiBins [expr {$i * $pi/30.0}]
    }
    foreach eta {-4.0 -3.5 -2.95880652 -2.68264484 -2.46773612 -2.29224349 -2.14432155 -2.01681569 -1.90506801 -1.80587261 -1.71692581 -1.63651428 -1.56332731 -1.49633825 -1.43472677 -1.37782606 -1.325086   -1.27604684 -1.23031998 -1.18757364 -1.14752205 -1.10991713 -1.07454199 -1.04120583 -1.00} {
	add EtaPhiBins $eta $PhiBins
    }
    
    foreach eta {1.0 1.04 1.075 1.1099 1.14752205 1.18757364 1.23031998 1.27604684 1.325086 1.37782606 1.43472677 1.49633825 1.56332731 1.63651428 1.71692581 1.80587261 1.90506801 2.01681569 2.14432155 2.29224349 2.46773612 2.68264484 2.95880652 3.5 4.0} {
	add EtaPhiBins $eta $PhiBins
    }

  # ATHENA model
  add EnergyFraction {0} {1.0}
  # energy fractions for e, gamma and pi0
  add EnergyFraction {11} {0.0}
  add EnergyFraction {22} {0.0}
  add EnergyFraction {111} {0.0}
  # energy fractions for muon, neutrinos and neutralinos
  add EnergyFraction {12} {0.0}
  add EnergyFraction {13} {0.0}
  add EnergyFraction {14} {0.0}
  add EnergyFraction {16} {0.0}
  add EnergyFraction {1000022} {0.0}
  add EnergyFraction {1000023} {0.0}
  add EnergyFraction {1000025} {0.0}
  add EnergyFraction {1000035} {0.0}
  add EnergyFraction {1000045} {0.0}
  # energy fractions for K0short and Lambda
  add EnergyFraction {310} {0.7}
  add EnergyFraction {3122} {0.7}

  # Refer to ePIC HCal resolution (presented at Hard Probes 2023)
  # https://wiki.bnl.gov/EPIC/index.php?title=File:HCalResolution_SinglePion.pdf
  # set HCalResolutionFormula {resolution formula as a function of eta and energy}
  set ResolutionFormula {    (eta <= 1.0 && eta > -1.0)                       * sqrt(energy^2*0.145^2 + energy*0.75^2)+
                             (eta <= 4.0 && eta >  1.0)                       * sqrt(energy^2*0.055^2 + energy*0.443^2)
  }

}


#################
# Electron filter
#################

module PdgCodeFilter ElectronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray electrons
  set Invert true
  add PdgCode {11}
  add PdgCode {-11}
}

######################
# ChargedHadronFilter
######################

module PdgCodeFilter ChargedHadronFilter {
  set InputArray HCal/eflowTracks
  set OutputArray chargedHadrons

  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}

###################################################
# Tower Merger (in case not using e-flow algorithm)
###################################################

module Merger Calorimeter {
# add InputArray InputArray
  add InputArray ECal/ecalTowers
  add InputArray HCal/hcalTowers
  set OutputArray towers
}

####################
# Energy flow merger
####################

module Merger EFlowMerger {
# add InputArray InputArray
  add InputArray HCal/eflowTracks
  add InputArray ECal/eflowPhotons
  add InputArray HCal/eflowNeutralHadrons
  set OutputArray eflow
}

######################
# EFlowFilter
######################

module PdgCodeFilter EFlowFilter {
  set InputArray EFlowMerger/eflow
  set OutputArray eflow

  add PdgCode {11}
  add PdgCode {-11}
  add PdgCode {13}
  add PdgCode {-13}
}

###################
# Photon efficiency
###################

module Efficiency PhotonEfficiency {
  set InputArray ECal/eflowPhotons
  set OutputArray photons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for photons
    set EfficiencyFormula { 1}
}

##################
# Photon isolation
##################

module Isolation PhotonIsolation {
  set CandidateInputArray PhotonEfficiency/photons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray photons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
}


#####################
# Electron efficiency
#####################

module Efficiency ElectronEfficiency {
  set InputArray ElectronFilter/electrons
  set OutputArray electrons

  # set EfficiencyFormula {efficiency formula as a function of eta and pt}

  # efficiency formula for electrons
    set EfficiencyFormula {1}
}

####################
# Electron isolation
####################

module Isolation ElectronIsolation {
  set CandidateInputArray ElectronEfficiency/electrons
  set IsolationInputArray EFlowFilter/eflow

  set OutputArray electrons

  set DeltaRMax 0.5

  set PTMin 0.5

  set PTRatioMax 0.12
}

###################
# Missing ET merger
###################

module Merger MissingET {
# add InputArray InputArray
  add InputArray EFlowMerger/eflow
  set MomentumOutputArray momentum
}

##################
# Scalar HT merger
##################

module Merger ScalarHT {
# add InputArray InputArray
  add InputArray UniqueObjectFinder/jets
  add InputArray UniqueObjectFinder/electrons
  add InputArray UniqueObjectFinder/photons

  set EnergyOutputArray energy
}

#####################
# Neutrino Filter
#####################

module PdgCodeFilter NeutrinoFilter {

  set InputArray Delphes/stableParticles
  set OutputArray filteredParticles

  set PTMin 0.0

  add PdgCode {12}
  add PdgCode {14}
  add PdgCode {16}
  add PdgCode {-12}
  add PdgCode {-14}
  add PdgCode {-16}

}

#####################
# MC truth jet finder
#####################

module FastJetFinder GenJetFinder {
  set InputArray NeutrinoFilter/filteredParticles

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.0

  set JetPTMin 3.0
}

#########################
# Gen Missing ET merger
########################

module Merger GenMissingET {
# add InputArray InputArray
  add InputArray NeutrinoFilter/filteredParticles
  set MomentumOutputArray momentum
}

############
# Jet finder
############

module FastJetFinder FastJetFinder {
#  set InputArray Calorimeter/towers
  set InputArray EFlowMerger/eflow

  set OutputArray jets

  # algorithm: 1 CDFJetClu, 2 MidPoint, 3 SIScone, 4 kt, 5 Cambridge/Aachen, 6 antikt
  set JetAlgorithm 6
  set ParameterR 1.0

  set ComputeNsubjettiness 1
  set Beta 1.0
  set AxisMode 4

  set ComputeTrimming 1
  set RTrim 0.4
  set PtFracTrim 0.20
  #set PtFracTrim 0.05

  set ComputePruning 1
  set ZcutPrun 0.1
  set RcutPrun 0.5
  set RPrun 0.8

  set ComputeSoftDrop 1
  set BetaSoftDrop 0.0
  set SymmetryCutSoftDrop 0.1
  set R0SoftDrop 0.8

  set JetPTMin 3.0}

##################
# Jet Energy Scale
##################

module EnergyScale JetEnergyScale {
  set InputArray FastJetFinder/jets
  set OutputArray jets

  # scale formula for jets (do not apply it)
  set ScaleFormula {1.0}
}

########################
# Jet Flavor Association
########################

module JetFlavorAssociation JetFlavorAssociation {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray JetEnergyScale/jets

  set DeltaR 0.5
  set PartonPTMin 4.0
}

module JetFlavorAssociation GenJetFlavorAssociation {
  set PartonInputArray Delphes/partons
  set ParticleInputArray Delphes/allParticles
  set ParticleLHEFInputArray Delphes/allParticlesLHEF
  set JetInputArray GenJetFinder/jets

  set DeltaR 0.5
  set PartonPTMin 1.0
  set PartonEtaMax 4.0

}

#####################################################
# Find uniquely identified photons/electrons/tau/jets
#####################################################

module UniqueObjectFinder UniqueObjectFinder {
# earlier arrays take precedence over later ones
# add InputArray InputArray OutputArray
  add InputArray PhotonIsolation/photons photons
  add InputArray ElectronIsolation/electrons electrons
  add InputArray JetEnergyScale/jets jets
}

############################
# b-tagging (track counting)
############################

module TrackCountingBTagging TrackCountingBTagging {
    set JetInputArray JetEnergyScale/jets
    set TrackInputArray HCal/eflowTracks
    set BitNumber 0
    # maximum distance between jet and track
    set DeltaR 0.5
    # minimum pt of tracks
    set TrackPtMin 1.0
    # maximum transverse impact parameter (in mm)
    set TrackIPMax 3
    # minimum ip significance for the track to be counted
    set SigMin 2.0
    set Use3D true
    # alternate setting for 2D IP (default)
    #  set SigMin 1.3
    #  set Use3D false
    # minimum number of tracks (high efficiency n=2, high purity n=3)
    #set Ntracks 3
}

##################
# Particle ID Systems
##################
# ATHENA model
source pfRICH_0.25mrad.tcl
#source pfRICH_0.0mrad.tcl
source barrelDIRC_epid.tcl
source barrelDIRC_hpid.tcl
#source BTOF.epid.t020ps.tcl
#source BTOF.hpid.t020ps.tcl
source BTOF.epid.t030ps.tcl
source BTOF.hpid.t030ps.tcl
#source dualRICH_aerogel_0.00mrad.tcl
#source dualRICH_aerogel_0.50mrad.tcl
source dualRICH_aerogel_0.25mrad.tcl
source dualRICH_c2f6_0.25mrad.tcl
#source dualRICH_c2f6_0.00mrad.tcl
#source dualRICH_c2f6_0.50mrad.tcl

##################
# ROOT tree writer
##################

# tracks, towers and eflow objects are not stored by default in the output.
# if needed (for jet constituent or other studies), uncomment the relevant
# "add Branch ..." lines.

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle
  add Branch BeamSpotFilter/beamSpotParticle BeamSpot GenParticle

  add Branch TrackMerger/tracks Track Track
  add Branch Calorimeter/towers Tower Tower

  add Branch HCal/eflowTracks EFlowTrack Track
  add Branch ECal/eflowPhotons EFlowPhoton Tower
  add Branch HCal/eflowNeutralHadrons EFlowNeutralHadron Tower

  add Branch pfRICH/tracks pfRICHTrack Track
  
  add Branch barrelDIRC_epid/tracks barrelDIRC_epidTrack Track
  add Branch barrelDIRC_hpid/tracks barrelDIRC_hpidTrack Track

  add Branch BTOF_e/tracks BTOF_eTrack Track
  add Branch BTOF_h/tracks BTOF_hTrack Track
  add Branch dualRICH_aerogel/tracks dualRICHagTrack Track
  add Branch dualRICH_c2f6/tracks dualRICHcfTrack Track

  add Branch GenJetFinder/jets GenJet Jet
  add Branch GenMissingET/momentum GenMissingET MissingET

  add Branch UniqueObjectFinder/jets Jet Jet
  add Branch UniqueObjectFinder/electrons Electron Electron
  add Branch UniqueObjectFinder/photons Photon Photon

  add Branch MissingET/momentum MissingET MissingET
  add Branch ScalarHT/energy ScalarHT ScalarHT
}
