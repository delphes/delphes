#!/usr/bin/env python3

import sys

from urllib.request import Request, urlopen
import zipfile

import ROOT

import ProMCHeader_pb2
import ProMC_pb2

################################################################################
# HttpFile class is written by Eli Carter (retracile)
# http://stackoverflow.com/a/7852229


class HttpFile(object):
    def __init__(self, url):
        self.url = url
        self.offset = 0
        self.length = -1

    def seekable(self):
        return True

    def size(self):
        if self.length < 0:
            f = urlopen(self.url)
            self.length = int(f.headers["Content-length"])
        return self.length

    def read(self, count=-1):
        req = Request(self.url)
        if count < 0:
            end = self.size() - 1
        else:
            end = self.offset + count - 1
        req.headers["Range"] = "bytes=%s-%s" % (self.offset, end)
        f = urlopen(req)
        data = f.read()
        # FIXME: should check that we got the range expected, etc.
        chunk = len(data)
        if count >= 0:
            assert chunk == count
        self.offset += chunk
        return data

    def seek(self, offset, whence=0):
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


def ConvertInput(input, momentumUnit, lengthUnit, branch, factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray):

    pdg = ROOT.TDatabasePDG.Instance()

    data = ProMC_pb2.ProMCEvent()
    data.ParseFromString(input)

    event = data.event
    particles = data.particles

    element = branch.NewEntry()

    element.Number = event.Number if event.Number else 0
    element.ProcessID = event.Process_ID if event.Process_ID else 0
    element.MPI = event.MPI if event.MPI else 0
    element.Weight = event.Weight if event.Weight else 0
    element.Scale = event.Scale if event.Scale else 0
    element.AlphaQED = event.Alpha_QED if event.Alpha_QED else 0
    element.AlphaQCD = event.Alpha_QCD if event.Alpha_QCD else 0
    element.ID1 = event.ID1 if event.ID1 else 0
    element.ID2 = event.ID2 if event.ID2 else 0
    element.X1 = event.X1 if event.X1 else 0
    element.X2 = event.X2 if event.X2 else 0
    element.ScalePDF = event.Scale_PDF if event.Scale_PDF else 0
    element.PDF1 = event.PDF1 if event.PDF1 else 0
    element.PDF2 = event.PDF2 if event.PDF2 else 0

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

        candidate.Status = status

        candidate.M1 = particles.mother1[i] if particles.mother1[i] < sys.maxsize else -1
        candidate.M2 = particles.mother2[i] if particles.mother2[i] < sys.maxsize else -1

        candidate.D1 = particles.daughter1[i] if particles.daughter1[i] < sys.maxsize else -1
        candidate.D2 = particles.daughter2[i] if particles.daughter2[i] < sys.maxsize else -1

        candidate.Mass = mass

        candidate.Momentum.SetXYZM(px, py, pz, mass)

        candidate.Position.SetXYZT(x, y, z, t)

        pdgParticle = pdg.GetParticle(pid)

        candidate.Charge = int(pdgParticle.Charge() / 3.0) if pdgParticle else -999

        allParticleOutputArray.Add(candidate)

        if not pdgParticle:
            continue

        pdgCode = abs(pid)

        if status == 1:
            stableParticleOutputArray.Add(candidate)
        elif pdgCode <= 5 or pdgCode == 21 or pdgCode == 15:
            partonOutputArray.Add(candidate)


################################################################################

if len(sys.argv) < 4:
    print(" Usage: DelphesProMC.py config_file output_file input_file(s)")
    sys.exit(1)

ROOT.gSystem.Load("libDelphes")

outputFile = ROOT.TFile(sys.argv[2], "CREATE")

if not outputFile.IsOpen():
    print("** ERROR: can't open", sys.argv[2])
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
    print("** Reading", fileName)

    if fileName.startswith("http://") or fileName.startswith("https://"):
        file = HttpFile(fileName)
    else:
        file = fileName

    with zipfile.ZipFile(file) as zf:

        numberOfEvents = len(zf.namelist())
        if numberOfEvents <= 0:
            continue

        progressBar = ROOT.ExRootProgressBar(numberOfEvents - 1)

        # retrive information from the header file
        header = ProMCHeader_pb2.ProMCHeader()
        header.ParseFromString(zf.read("header"))
        momentumUnit = float(header.MomentumUnit)
        lengthUnit = float(header.LengthUnit)

        modularDelphes.Clear()
        treeWriter.Clear()
        eventCounter = 0
        for name in zf.namelist():
            eventCounter += 1
            if not name.isdigit():
                continue

            ConvertInput(zf.read(name), momentumUnit, lengthUnit, branchEvent, factory, allParticleOutputArray, stableParticleOutputArray, partonOutputArray)

            modularDelphes.ProcessTask()

            treeWriter.Fill()

            modularDelphes.Clear()
            treeWriter.Clear()
            progressBar.Update(eventCounter)

        progressBar.Update(eventCounter, eventCounter, ROOT.kTRUE)
        progressBar.Finish()

modularDelphes.FinishTask()
treeWriter.Write()

del treeWriter
