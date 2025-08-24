import ROOT
import libcppyy as _backend

ROOT.gSystem.Load("libDelphes")

#####################################################
### Definition of the string conversion methods   ###
#####################################################


def _Event__str__(self):
    theString = "Event number %i\n" % self.Number
    theString += "Read in %f and processed in %f.\n" % (self.ReadTime, self.ProcTime)
    return theString


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("Event").__str__ = _Event__str__
else:
    _backend.CreateScopeProxy("Event").__str__ = _Event__str__


def _LHCOEvent__str__(self):
    theString = "Trigger word: %b\n" % self.Trigger
    return theString


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("LHCOEvent").__str__ = _LHCOEvent__str__
else:
    _backend.CreateScopeProxy("LHCOEvent").__str__ = _LHCOEvent__str__


def _LHEFEvent__str__(self):
    hepmcstr = "Process ID: %i\n" % self.ProcessID
    hepmcstr += "ScalePDF = %f\n" % self.ScalePDF
    hepmcstr += "Alpha QED = %f\n" % self.AlphaQED
    hepmcstr += "Alpha DCQ = %f\n" % self.AlphaQCD
    hepmcstr += "Event weight = %f\n" % self.Weight
    return hepmcstr


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("LHEFEvent").__str__ = _LHEFEvent__str__
else:
    _backend.CreateScopeProxy("LHEFEvent").__str__ = _LHEFEvent__str__


def _HepMCEvent__str__(self):
    hepmcstr = "Process ID: %i\n" % self.ProcessID
    hepmcstr += "Energy scale: %f\n" % self.Scale
    hepmcstr += "Alpha QED = %f\n" % self.AlphaQED
    hepmcstr += "Alpha QCD = %f\n" % self.AlphaQCD
    hepmcstr += "Number of multi parton interactions: %f\n" % self.MPI
    hepmcstr += """Q-scale used in evaluation of PDFs (in GeV): %f\n""" % self.ScalePDF
    hepmcstr += """Fraction of beam momentum carried by %i parton ("beam side"): %f\n""" % (self.ID1, self.X1)
    hepmcstr += """Fraction of beam momentum carried by %i parton ("target side"): %f\n""" % (self.ID2, self.X2)
    hepmcstr += "PDF (id1, x1, Q) = %f\n" % self.PDF1
    hepmcstr += "PDF (id2, x2, Q) = %f\n" % self.PDF2
    return hepmcstr


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("HepMCEvent").__str__ = _HepMCEvent__str__
else:
    _backend.CreateScopeProxy("HepMCEvent").__str__ = _HepMCEvent__str__


def _GenParticle__str__(self):
    thestring = "Particle HEP ID number = %d\n" % self.PID
    thestring += "Particle status: %d\n" % self.Status
    if self.IsPU:
        thestring += "That particle is from PU.\n"
    thestring += "Particle 1st mother: %d\n" % self.M1
    thestring += "Particle 2nd mother: %d\n" % self.M2
    thestring += "Particle 1st daughter: %d\n" % self.D1
    thestring += "Particle last daughter: %d\n" % self.D2
    thestring += "(Px,Py,Pz,E) = (%f,%f,%f,%f)\n" % (self.Px, self.Py, self.Pz, self.E)
    thestring += "(Pt,Eta,Phi) = (%f,%f,%f)\n" % (self.PT, self.Eta, self.Phi)
    thestring += "Rapidity = %f\n" % self.Rapidity
    thestring += "Charge: %d\n" % self.Charge
    thestring += "Mass: %f\n" % self.Mass
    thestring += "Vertex position: (x,y,z,t) = (%f,%f,%f,%f)\n" % (self.X, self.Y, self.Z, self.T)
    return thestring


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("GenParticle").__str__ = _GenParticle__str__
else:
    _backend.CreateScopeProxy("GenParticle").__str__ = _GenParticle__str__


def _GenParticle_printDecay(self, db, particles, pre=""):
    thestring = "%s (%f,%f,%f,%f) (status %d)\n" % (db.GetParticle(self.PID).GetName(), self.Px, self.Py, self.Pz, self.E, self.Status)
    for daughter in range(self.D1, self.D2 + 1):
        if daughter == -1:
            break
        if daughter < self.D2:
            thestring += "%s+->%s" % (pre, particles.At(daughter).printDecay(db, particles, pre + "|   "))
        else:
            thestring += "%s+->%s" % (pre, particles.At(daughter).printDecay(db, particles, pre + "    "))
    return thestring


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("GenParticle").printDecay = _GenParticle_printDecay
else:
    _backend.CreateScopeProxy("GenParticle").printDecay = _GenParticle_printDecay


def _MissingET__str__(self):
    thestring = "Mising transverse energy: %f\n" % self.MET
    thestring += "Mising energy azimuthal angle: %f\n" % self.Phi
    return thestring


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("MissingET").__str__ = _MissingET__str__
else:
    _backend.CreateScopeProxy("MissingET").__str__ = _MissingET__str__


def _ScalarHT__str__(self):
    thestring = "Scalar sum of transverse momenta: %f\n" % self.HT
    return thestring


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("ScalarHT").__str__ = _ScalarHT__str__
else:
    _backend.CreateScopeProxy("ScalarHT").__str__ = _ScalarHT__str__


def _Photon__str__(self):
    thestring = "(Pt,Eta,Phi) =(%f,%f,%f)\n" % (self.PT, self.Eta, self.Phi)
    return thestring


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("Photon").__str__ = _Photon__str__
else:
    _backend.CreateScopeProxy("Photon").__str__ = _Photon__str__


def _Electron__str__(self):
    thestring = "(Pt,Eta,Phi) =(%f,%f,%f)\n" % (self.PT, self.Eta, self.Phi)
    thestring += "charge: %d\n" % self.Charge
    return thestring


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("Electron").__str__ = _Electron__str__
else:
    _backend.CreateScopeProxy("Electron").__str__ = _Electron__str__


def _Muon__str__(self):
    thestring = "(Pt,Eta,Phi) =(%f,%f,%f)\n" % (self.PT, self.Eta, self.Phi)
    thestring += "charge: %d\n" % self.Charge
    return thestring


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("Muon").__str__ = _Muon__str__
else:
    _backend.CreateScopeProxy("Muon").__str__ = _Muon__str__


def _Jet__str__(self):
    thestring = "(Pt,Eta,Phi) =(%f,%f,%f)\n" % (self.PT, self.Eta, self.Phi)
    thestring += "Charge: %d\n" % self.Charge
    thestring += "Mass: %f\n" % self.Mass
    if self.BTag:
        thestring += "This jet is b-tagged.\n"
    if self.TauTag:
        thestring += "This jet is tagged as a tau-jet.\n"
    thestring += "Jet radius in (eta,phi) = (%f,%f)\n" % (self.DeltaEta, self.DeltaPhi)
    thestring += "Number of constituents: %d\n" % self.Constituents.GetEntries()
    thestring += "H/E ratio in the calorimeters: %f\n" % self.EhadOverEem
    nch = 0
    ech = 0.0
    for constit in self.Constituents:
        if constit.Charge != 0:
            fourvector = ROOT.TLorentzVector()
            fourvector.SetPtEtaPhiM(constit.PT, constit.Eta, constit.Phi, 0.0)
            nch += 1
            ech += fourvector.E()
    thestring += "Charged Energy fraction: %f\n" % (ech / self.P4().E())
    thestring += "Charged Multiplicity: %d\n" % nch
    return thestring


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("Jet").__str__ = _Jet__str__
else:
    _backend.CreateScopeProxy("Jet").__str__ = _Jet__str__


def _Track__str__(self):
    thestring = "Particle HEP ID number = %d\n" % self.PID
    thestring += "(Pt,Eta,Phi) = (%f,%f,%f)\n" % (self.PT, self.Eta, self.Phi)
    thestring += "Charge: %d\n" % self.Charge
    thestring += "At calorimeter surface: (Pt,Eta,Phi) = (%f,%f,%f)\n" % (self.PT, self.EtaOuter, self.PhiOuter)
    thestring += "Vertex position: (x,y,z) = (%f,%f,%f)\n" % (self.X, self.Y, self.Z)
    thestring += "Position at calorimer surface: (x,y,z) = (%f,%f,%f)\n" % (self.XOuter, self.YOuter, self.ZOuter)
    return thestring


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("Track").__str__ = _Track__str__
else:
    _backend.CreateScopeProxy("Track").__str__ = _Track__str__


def _Tower__str__(self):
    thestring = "(Et,Eta,Phi) = (%f,%f,%f)\n" % (self.ET, self.Eta, self.Phi)
    thestring += "(E,Eem,Ehad) = (%f,%f,%f)\n" % (self.E, self.Eem, self.Ehad)
    thestring += "Edges: %f, %f, %f, %f\n" % (self.Edges[0], self.Edges[1], self.Edges[2], self.Edges[3])
    return thestring


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("Tower").__str__ = _Tower__str__
else:
    _backend.CreateScopeProxy("Tower").__str__ = _Tower__str__

#####################################################
### Definition of additional ROOT methods         ###
#####################################################


def _lorentzVector__str__(self):
    theString = "(pt, eta, phi) = (%f,%f,%f)\n" % (self.Pt(), self.Eta(), self.Phi())
    theString += "mass = %f, p = %f, mt = %f\n" % (self.M(), self.P(), self.Mt())
    return theString


if hasattr(_backend, "MakeRootClass"):
    _backend.MakeRootClass("TLorentzVector").__str__ = _lorentzVector__str__
else:
    _backend.CreateScopeProxy("TLorentzVector").__str__ = _lorentzVector__str__
