import sys

import urllib2
import cStringIO
import zipfile 

import ROOT

import ProMCHeader_pb2
import ProMC_pb2

################################################################################

class HttpFile(object):
  def __init__(self, url):
    self.url = url
    self.offset = 0
    self.length = -1

  def size(self):
    if self.length < 0:
      f = urllib2.urlopen(self.url)
      self.length = int(f.headers["Content-length"])
    return self.length

  def read(self, count = -1):
    req = urllib2.Request(self.url)
    if count < 0:
      end = self.size() - 1
    else:
      end = self.offset + count - 1
    req.headers['Range'] = "bytes=%s-%s" % (self.offset, end)
    f = urllib2.urlopen(req)
    data = f.read()
    # FIXME: should check that we got the range expected, etc.
    chunk = len(data)
    if count >= 0:
      assert chunk == count
    self.offset += chunk
    return data

  def seek(self, offset, whence = 0):
    if whence == 0:
      self.offset = offset
    elif whence == 1:
      self.offset += offset
    elif whence == 2:
      self.offset = self.size() + offset
    else:
      raise Exception("Invalid whence")

  def tell(self):
    return self.offset

################################################################################

def ConvertInput(name, momentumUnit, lengthUnit, branch, factory,
  allParticleOutputArray, stableParticleOutputArray, partonOutputArray):

  pdg = ROOT.TDatabasePDG.Instance()

  data = ProMC_pb2.ProMCEvent()
  data.ParseFromString(zip.read(name))
  
  event = data.event
  particles = data.particles

  element = branch.NewEntry()


# element.Number = event.Number
# element.ProcessID = event.Process_ID
# element.MPI = event.MPI
# element.Weight = event.Weight
# element.Scale = event.Scale
# element.AlphaQED = event.Alpha_QED
# element.AlphaQCD = event.Alpha_QCD
# element.ID1 = event.ID1
# element.ID2 = event.ID2
# element.X1 = event.X1
# element.X2 = event.X2
# element.ScalePDF = event.Scale_PDF
# element.PDF1 = event.PDF1
# element.PDF2 = event.PDF2

  for i in range(len(particles.pdg_id)):
    pid = particles.pdg_id[i]
    status = particles.status[i]
    px = particles.Px[i]
    py = particles.Py[i]
    pz = particles.Pz[i]
    mass = particles.mass[i]

    x = particles.X[i]
    y = particles.Y[i]
    z = particles.Z[i]
    t = particles.T[i]

    candidate = factory.NewCandidate()

    candidate.PID = pid
    pdgCode = ROOT.TMath.Abs(candidate.PID)

    candidate.Status = status

#    candidate.M1 = particles.mother1[i]
#    candidate.M2 = particles.mother2[i]

#    candidate.D1 = particles.daughter1[i]
#    candidate.D2 = particles.daughter2[i]

    candidate.Mass = mass

    candidate.Momentum.SetXYZM(px, py, pz, mass)

    candidate.Position.SetXYZT(x, y, z, t)

    pdgParticle = pdg.GetParticle(pid)
    if pdgParticle == None:
      candidate.Charge = -999
      continue
    else:
      candidate.Charge = int(pdgParticle.Charge()/3.0)

    allParticleOutputArray.Add(candidate)

    if status == 1:
      stableParticleOutputArray.Add(candidate)
    elif pdgCode <= 5 or pdgCode == 21 or pdgCode == 15:
      partonOutputArray.Add(candidate)

################################################################################

if len(sys.argv) < 2:
  print " Usage: DelphesProMC.py config_file output_file input_file(s)"
  sys.exit(1)

ROOT.gSystem.Load('libDelphes')

outputFile = ROOT.TFile(sys.argv[2], "CREATE")

if not outputFile.IsOpen():
  print "** ERROR: can't open ", sys.argv[2]
  sys.exit(1)

treeWriter = ROOT.ExRootTreeWriter(outputFile, "Delphes")

branchEvent = treeWriter.NewBranch("Event", ROOT.HepMCEvent.Class())

confReader = ROOT.ExRootConfReader()
confReader.ReadFile(sys.argv[1])

modularDelphes = ROOT.Delphes("Delphes")
modularDelphes.SetConfReader(confReader)
modularDelphes.SetTreeWriter(treeWriter)

factory = modularDelphes.GetFactory()
allParticleOutputArray = modularDelphes.ExportArray("allParticles")
stableParticleOutputArray = modularDelphes.ExportArray("stableParticles")
partonOutputArray = modularDelphes.ExportArray("partons")

modularDelphes.InitTask()

for fileName in sys.argv[3:]:
  print "** Reading " + fileName

  file = HttpFile(fileName)
  zip = zipfile.ZipFile(file)

  # retrive information from the header file
  header = ProMCHeader_pb2.ProMCHeader()
  header.ParseFromString(zip.read("header"))
  momentumUnit = float(header.MomentumUnit)
  lengthUnit = float(header.LengthUnit)

  modularDelphes.Clear()
  treeWriter.Clear()
  for name in zip.namelist():
    if name.isdigit():
      ConvertInput(name, momentumUnit, lengthUnit, branchEvent, factory,
        allParticleOutputArray, stableParticleOutputArray, partonOutputArray)

      modularDelphes.ProcessTask()

      treeWriter.Fill()

      modularDelphes.Clear()
      treeWriter.Clear()

  zip.close()

modularDelphes.FinishTask()
treeWriter.Write()

del treeWriter
