
#
# Makefile for ExRootAnalysis
#
# Author: P. Demin - UCL, Louvain-la-Neuve
#
# multi-platform configuration is taken from ROOT (root/test/Makefile.arch)
#

include doc/Makefile.arch

ifeq ($(ARCH),macosx64)
UNDEFOPT = dynamic_lookup
endif

SrcSuf = cc

CXXFLAGS += $(ROOTCFLAGS) -Wno-write-strings -D_FILE_OFFSET_BITS=64 -DDROP_CGAL -I. -Iexternal -Iexternal/tcl
DELPHES_LIBS = $(shell $(RC) --libs) -lEG $(SYSLIBS)
DISPLAY_LIBS = $(shell $(RC) --evelibs) $(SYSLIBS)

ifneq ($(PYTHIA8),)
HAS_PYTHIA8 = true
CXXFLAGS += -DHAS_PYTHIA8 -I$(PYTHIA8)/include
DELPHES_LIBS += -L$(PYTHIA8)/lib -lpythia8 -llhapdfdummy
else
ifneq ($(PYTHIA8DATA),)
HAS_PYTHIA8 = true
CXXFLAGS += -DHAS_PYTHIA8 -I$(PYTHIA8DATA)/../include
DELPHES_LIBS += -L$(PYTHIA8DATA)/../lib -lpythia8 -llhapdfdummy
else
ifneq ($(PYTHIA8DATA),)
HAS_PYTHIA8 = true
CXXFLAGS += -DHAS_PYTHIA8 -I$(PYTHIA8DATA)/../include
DELPHES_LIBS += -L$(PYTHIA8DATA)/../lib/archive -lpythia8 -llhapdfdummy
endif
endif


###

DELPHES = libDelphes.$(DllSuf)
DELPHESLIB = libDelphes.lib

DISPLAY = libDelphesDisplay.$(DllSuf)
DISPLAYLIB = libDelphesDisplay.lib

VERSION = $(shell cat VERSION)
DISTDIR = Delphes-$(VERSION)
DISTTAR = $(DISTDIR).tar.gz

all:


DelphesLHEF$(ExeSuf): \
	tmp/readers/DelphesLHEF.$(ObjSuf)

tmp/readers/DelphesLHEF.$(ObjSuf): \
	readers/DelphesLHEF.cpp \
	modules/Delphes.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesLHEFReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
DelphesHepMC$(ExeSuf): \
	tmp/readers/DelphesHepMC.$(ObjSuf)

tmp/readers/DelphesHepMC.$(ObjSuf): \
	readers/DelphesHepMC.cpp \
	modules/Delphes.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesHepMCReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
DelphesSTDHEP$(ExeSuf): \
	tmp/readers/DelphesSTDHEP.$(ObjSuf)

tmp/readers/DelphesSTDHEP.$(ObjSuf): \
	readers/DelphesSTDHEP.cpp \
	modules/Delphes.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesSTDHEPReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
lhco2root$(ExeSuf): \
	tmp/converters/lhco2root.$(ObjSuf)

tmp/converters/lhco2root.$(ObjSuf): \
	converters/lhco2root.cpp \
	modules/Delphes.h \
	classes/DelphesStream.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
root2pileup$(ExeSuf): \
	tmp/converters/root2pileup.$(ObjSuf)

tmp/converters/root2pileup.$(ObjSuf): \
	converters/root2pileup.cpp \
	classes/DelphesClasses.h \
	classes/DelphesPileUpWriter.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootProgressBar.h
root2lhco$(ExeSuf): \
	tmp/converters/root2lhco.$(ObjSuf)

tmp/converters/root2lhco.$(ObjSuf): \
	converters/root2lhco.cpp \
	classes/DelphesClasses.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootProgressBar.h
stdhep2pileup$(ExeSuf): \
	tmp/converters/stdhep2pileup.$(ObjSuf)

tmp/converters/stdhep2pileup.$(ObjSuf): \
	converters/stdhep2pileup.cpp \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesSTDHEPReader.h \
	classes/DelphesPileUpWriter.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
hepmc2pileup$(ExeSuf): \
	tmp/converters/hepmc2pileup.$(ObjSuf)

tmp/converters/hepmc2pileup.$(ObjSuf): \
	converters/hepmc2pileup.cpp \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesHepMCReader.h \
	classes/DelphesPileUpWriter.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
pileup2root$(ExeSuf): \
	tmp/converters/pileup2root.$(ObjSuf)

tmp/converters/pileup2root.$(ObjSuf): \
	converters/pileup2root.cpp \
	classes/DelphesStream.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesPileUpReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootProgressBar.h
Example1$(ExeSuf): \
	tmp/examples/Example1.$(ObjSuf)

tmp/examples/Example1.$(ObjSuf): \
	examples/Example1.cpp \
	classes/DelphesClasses.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootUtilities.h
EXECUTABLE =  \
	DelphesLHEF$(ExeSuf) \
	DelphesHepMC$(ExeSuf) \
	DelphesSTDHEP$(ExeSuf) \
	lhco2root$(ExeSuf) \
	root2pileup$(ExeSuf) \
	root2lhco$(ExeSuf) \
	stdhep2pileup$(ExeSuf) \
	hepmc2pileup$(ExeSuf) \
	pileup2root$(ExeSuf) \
	Example1$(ExeSuf)

EXECUTABLE_OBJ =  \
	tmp/readers/DelphesLHEF.$(ObjSuf) \
	tmp/readers/DelphesHepMC.$(ObjSuf) \
	tmp/readers/DelphesSTDHEP.$(ObjSuf) \
	tmp/converters/lhco2root.$(ObjSuf) \
	tmp/converters/root2pileup.$(ObjSuf) \
	tmp/converters/root2lhco.$(ObjSuf) \
	tmp/converters/stdhep2pileup.$(ObjSuf) \
	tmp/converters/hepmc2pileup.$(ObjSuf) \
	tmp/converters/pileup2root.$(ObjSuf) \
	tmp/examples/Example1.$(ObjSuf)

tmp/classes/ClassesDict.$(SrcSuf): \
	classes/ClassesLinkDef.h \
	classes/DelphesModule.h \
	classes/DelphesFactory.h \
	classes/SortableObject.h \
	classes/DelphesClasses.h
tmp/modules/ModulesDict.$(SrcSuf): \
	modules/ModulesLinkDef.h \
	modules/Delphes.h \
	modules/FastJetFinder.h \
	modules/ParticlePropagator.h \
	modules/Efficiency.h \
	modules/EnergySmearing.h \
	modules/MomentumSmearing.h \
	modules/Calorimeter.h \
	modules/Isolation.h \
	modules/EnergyScale.h \
	modules/UniqueObjectFinder.h \
	modules/BTagging.h \
	modules/TauTagging.h \
	modules/TreeWriter.h \
	modules/Merger.h \
	modules/LeptonDressing.h \
	modules/PileUpMerger.h \
	modules/JetPileUpSubtractor.h \
	modules/TrackPileUpSubtractor.h \
	modules/ConstituentFilter.h \
	modules/StatusPidFilter.h \
	modules/Cloner.h \
	modules/Weighter.h \
	modules/ExampleModule.h \
	modules/PileUpMergerPythia8.h
tmp/external/ExRootAnalysis/ExRootAnalysisDict.$(SrcSuf): \
	external/ExRootAnalysis/ExRootAnalysisLinkDef.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootUtilities.h \
	external/ExRootAnalysis/ExRootClassifier.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootProgressBar.h \
	external/ExRootAnalysis/ExRootConfReader.h \
	external/ExRootAnalysis/ExRootTask.h
DELPHES_DICT =  \
	tmp/classes/ClassesDict.$(SrcSuf) \
	tmp/modules/ModulesDict.$(SrcSuf) \
	tmp/external/ExRootAnalysis/ExRootAnalysisDict.$(SrcSuf)

DELPHES_DICT_OBJ =  \
	tmp/classes/ClassesDict.$(ObjSuf) \
	tmp/modules/ModulesDict.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootAnalysisDict.$(ObjSuf)

tmp/display/DisplayDict.$(SrcSuf): \
	display/DisplayLinkDef.h \
	display/DelphesDisplay.h \
	display/DelphesCaloData.h
DISPLAY_DICT =  \
	tmp/display/DisplayDict.$(SrcSuf)

DISPLAY_DICT_OBJ =  \
	tmp/display/DisplayDict.$(ObjSuf)

tmp/classes/DelphesHepMCReader.$(ObjSuf): \
	classes/DelphesHepMCReader.$(SrcSuf) \
	classes/DelphesHepMCReader.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesStream.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/classes/DelphesLHEFReader.$(ObjSuf): \
	classes/DelphesLHEFReader.$(SrcSuf) \
	classes/DelphesLHEFReader.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesStream.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/classes/DelphesFactory.$(ObjSuf): \
	classes/DelphesFactory.$(SrcSuf) \
	classes/DelphesFactory.h \
	classes/DelphesClasses.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/classes/DelphesFormula.$(ObjSuf): \
	classes/DelphesFormula.$(SrcSuf) \
	classes/DelphesFormula.h
tmp/classes/DelphesSTDHEPReader.$(ObjSuf): \
	classes/DelphesSTDHEPReader.$(SrcSuf) \
	classes/DelphesSTDHEPReader.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/classes/DelphesStream.$(ObjSuf): \
	classes/DelphesStream.$(SrcSuf) \
	classes/DelphesStream.h
tmp/classes/DelphesPileUpWriter.$(ObjSuf): \
	classes/DelphesPileUpWriter.$(SrcSuf) \
	classes/DelphesPileUpWriter.h
tmp/classes/DelphesModule.$(ObjSuf): \
	classes/DelphesModule.$(SrcSuf) \
	classes/DelphesModule.h \
	classes/DelphesFactory.h \
	external/ExRootAnalysis/ExRootTreeReader.h \
	external/ExRootAnalysis/ExRootTreeBranch.h \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootResult.h
tmp/classes/DelphesClasses.$(ObjSuf): \
	classes/DelphesClasses.$(SrcSuf) \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/SortableObject.h
tmp/classes/DelphesPileUpReader.$(ObjSuf): \
	classes/DelphesPileUpReader.$(SrcSuf) \
	classes/DelphesPileUpReader.h
tmp/modules/LeptonDressing.$(ObjSuf): \
	modules/LeptonDressing.$(SrcSuf) \
	modules/LeptonDressing.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Efficiency.$(ObjSuf): \
	modules/Efficiency.$(SrcSuf) \
	modules/Efficiency.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/UniqueObjectFinder.$(ObjSuf): \
	modules/UniqueObjectFinder.$(SrcSuf) \
	modules/UniqueObjectFinder.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/ExampleModule.$(ObjSuf): \
	modules/ExampleModule.$(SrcSuf) \
	modules/ExampleModule.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/ParticlePropagator.$(ObjSuf): \
	modules/ParticlePropagator.$(SrcSuf) \
	modules/ParticlePropagator.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/EnergySmearing.$(ObjSuf): \
	modules/EnergySmearing.$(SrcSuf) \
	modules/EnergySmearing.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/FastJetFinder.$(ObjSuf): \
	modules/FastJetFinder.$(SrcSuf) \
	modules/FastJetFinder.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h \
	external/fastjet/PseudoJet.hh \
	external/fastjet/JetDefinition.hh \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/Selector.hh \
	external/fastjet/ClusterSequenceArea.hh \
	external/fastjet/tools/JetMedianBackgroundEstimator.hh \
	external/fastjet/plugins/SISCone/fastjet/SISConePlugin.hh \
	external/fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh \
	external/fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh
tmp/modules/PileUpMergerPythia8.$(ObjSuf): \
	modules/PileUpMergerPythia8.$(SrcSuf) \
	modules/PileUpMergerPythia8.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	classes/DelphesPileUpReader.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/StatusPidFilter.$(ObjSuf): \
	modules/StatusPidFilter.$(SrcSuf) \
	modules/StatusPidFilter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/ConstituentFilter.$(ObjSuf): \
	modules/ConstituentFilter.$(SrcSuf) \
	modules/ConstituentFilter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/EnergyScale.$(ObjSuf): \
	modules/EnergyScale.$(SrcSuf) \
	modules/EnergyScale.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/TauTagging.$(ObjSuf): \
	modules/TauTagging.$(SrcSuf) \
	modules/TauTagging.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Merger.$(ObjSuf): \
	modules/Merger.$(SrcSuf) \
	modules/Merger.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/BTagging.$(ObjSuf): \
	modules/BTagging.$(SrcSuf) \
	modules/BTagging.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/TreeWriter.$(ObjSuf): \
	modules/TreeWriter.$(SrcSuf) \
	modules/TreeWriter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/modules/Weighter.$(ObjSuf): \
	modules/Weighter.$(SrcSuf) \
	modules/Weighter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Delphes.$(ObjSuf): \
	modules/Delphes.$(SrcSuf) \
	modules/Delphes.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h \
	external/ExRootAnalysis/ExRootConfReader.h \
	external/ExRootAnalysis/ExRootTreeWriter.h
tmp/modules/Calorimeter.$(ObjSuf): \
	modules/Calorimeter.$(SrcSuf) \
	modules/Calorimeter.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Isolation.$(ObjSuf): \
	modules/Isolation.$(SrcSuf) \
	modules/Isolation.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/PileUpMerger.$(ObjSuf): \
	modules/PileUpMerger.$(SrcSuf) \
	modules/PileUpMerger.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	classes/DelphesPileUpReader.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/TrackPileUpSubtractor.$(ObjSuf): \
	modules/TrackPileUpSubtractor.$(SrcSuf) \
	modules/TrackPileUpSubtractor.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/Cloner.$(ObjSuf): \
	modules/Cloner.$(SrcSuf) \
	modules/Cloner.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/JetPileUpSubtractor.$(ObjSuf): \
	modules/JetPileUpSubtractor.$(SrcSuf) \
	modules/JetPileUpSubtractor.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/modules/MomentumSmearing.$(ObjSuf): \
	modules/MomentumSmearing.$(SrcSuf) \
	modules/MomentumSmearing.h \
	classes/DelphesClasses.h \
	classes/DelphesFactory.h \
	classes/DelphesFormula.h \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/external/ExRootAnalysis/ExRootFilter.$(ObjSuf): \
	external/ExRootAnalysis/ExRootFilter.$(SrcSuf) \
	external/ExRootAnalysis/ExRootFilter.h \
	external/ExRootAnalysis/ExRootClassifier.h
tmp/external/ExRootAnalysis/ExRootConfReader.$(ObjSuf): \
	external/ExRootAnalysis/ExRootConfReader.$(SrcSuf) \
	external/ExRootAnalysis/ExRootConfReader.h \
	external/tcl/tcl.h
tmp/external/ExRootAnalysis/ExRootTreeWriter.$(ObjSuf): \
	external/ExRootAnalysis/ExRootTreeWriter.$(SrcSuf) \
	external/ExRootAnalysis/ExRootTreeWriter.h \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/external/ExRootAnalysis/ExRootResult.$(ObjSuf): \
	external/ExRootAnalysis/ExRootResult.$(SrcSuf) \
	external/ExRootAnalysis/ExRootResult.h \
	external/ExRootAnalysis/ExRootUtilities.h
tmp/external/ExRootAnalysis/ExRootTreeBranch.$(ObjSuf): \
	external/ExRootAnalysis/ExRootTreeBranch.$(SrcSuf) \
	external/ExRootAnalysis/ExRootTreeBranch.h
tmp/external/ExRootAnalysis/ExRootTreeReader.$(ObjSuf): \
	external/ExRootAnalysis/ExRootTreeReader.$(SrcSuf) \
	external/ExRootAnalysis/ExRootTreeReader.h
tmp/external/ExRootAnalysis/ExRootTask.$(ObjSuf): \
	external/ExRootAnalysis/ExRootTask.$(SrcSuf) \
	external/ExRootAnalysis/ExRootTask.h \
	external/ExRootAnalysis/ExRootConfReader.h
tmp/external/ExRootAnalysis/ExRootProgressBar.$(ObjSuf): \
	external/ExRootAnalysis/ExRootProgressBar.$(SrcSuf) \
	external/ExRootAnalysis/ExRootProgressBar.h
tmp/external/ExRootAnalysis/ExRootUtilities.$(ObjSuf): \
	external/ExRootAnalysis/ExRootUtilities.$(SrcSuf) \
	external/ExRootAnalysis/ExRootUtilities.h
tmp/external/fastjet/Dnn2piCylinder.$(ObjSuf): \
	external/fastjet/Dnn2piCylinder.$(SrcSuf) \
	external/fastjet/internal/Dnn2piCylinder.hh
tmp/external/fastjet/GhostedAreaSpec.$(ObjSuf): \
	external/fastjet/GhostedAreaSpec.$(SrcSuf) \
	external/fastjet/GhostedAreaSpec.hh \
	external/fastjet/Error.hh
tmp/external/fastjet/ClusterSequenceActiveArea.$(ObjSuf): \
	external/fastjet/ClusterSequenceActiveArea.$(SrcSuf) \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/ClusterSequenceActiveArea.hh \
	external/fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh
tmp/external/fastjet/ClusterSequence_Delaunay.$(ObjSuf): \
	external/fastjet/ClusterSequence_Delaunay.$(SrcSuf) \
	external/fastjet/Error.hh \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/internal/Dnn4piCylinder.hh \
	external/fastjet/internal/Dnn3piCylinder.hh \
	external/fastjet/internal/Dnn2piCylinder.hh
tmp/external/fastjet/ClusterSequenceArea.$(ObjSuf): \
	external/fastjet/ClusterSequenceArea.$(SrcSuf) \
	external/fastjet/ClusterSequenceArea.hh
tmp/external/fastjet/Voronoi.$(ObjSuf): \
	external/fastjet/Voronoi.$(SrcSuf) \
	external/fastjet/internal/Voronoi.hh
tmp/external/fastjet/Selector.$(ObjSuf): \
	external/fastjet/Selector.$(SrcSuf) \
	external/fastjet/Selector.hh \
	external/fastjet/GhostedAreaSpec.hh
tmp/external/fastjet/ClusterSequence_TiledN2.$(ObjSuf): \
	external/fastjet/ClusterSequence_TiledN2.$(SrcSuf) \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/internal/MinHeap.hh
tmp/external/fastjet/ClusterSequenceActiveAreaExplicitGhosts.$(ObjSuf): \
	external/fastjet/ClusterSequenceActiveAreaExplicitGhosts.$(SrcSuf) \
	external/fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh
tmp/external/fastjet/LimitedWarning.$(ObjSuf): \
	external/fastjet/LimitedWarning.$(SrcSuf) \
	external/fastjet/LimitedWarning.hh
tmp/external/fastjet/CompositeJetStructure.$(ObjSuf): \
	external/fastjet/CompositeJetStructure.$(SrcSuf)
tmp/external/fastjet/DnnPlane.$(ObjSuf): \
	external/fastjet/DnnPlane.$(SrcSuf) \
	external/fastjet/internal/DnnPlane.hh
tmp/external/fastjet/ClusterSequence_DumbN3.$(ObjSuf): \
	external/fastjet/ClusterSequence_DumbN3.$(SrcSuf) \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequence.hh
tmp/external/fastjet/ClusterSequencePassiveArea.$(ObjSuf): \
	external/fastjet/ClusterSequencePassiveArea.$(SrcSuf) \
	external/fastjet/ClusterSequencePassiveArea.hh \
	external/fastjet/ClusterSequenceVoronoiArea.hh
tmp/external/fastjet/BasicRandom.$(ObjSuf): \
	external/fastjet/BasicRandom.$(SrcSuf) \
	external/fastjet/internal/BasicRandom.hh
tmp/external/fastjet/ClusterSequenceAreaBase.$(ObjSuf): \
	external/fastjet/ClusterSequenceAreaBase.$(SrcSuf) \
	external/fastjet/ClusterSequenceAreaBase.hh
tmp/external/fastjet/PseudoJetStructureBase.$(ObjSuf): \
	external/fastjet/PseudoJetStructureBase.$(SrcSuf) \
	external/fastjet/PseudoJetStructureBase.hh \
	external/fastjet/Error.hh \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/ClusterSequenceAreaBase.hh
tmp/external/fastjet/ClusterSequence.$(ObjSuf): \
	external/fastjet/ClusterSequence.$(SrcSuf) \
	external/fastjet/Error.hh \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/ClusterSequenceStructure.hh \
	external/fastjet/version.hh
tmp/external/fastjet/JetDefinition.$(ObjSuf): \
	external/fastjet/JetDefinition.$(SrcSuf) \
	external/fastjet/JetDefinition.hh \
	external/fastjet/Error.hh \
	external/fastjet/CompositeJetStructure.hh
tmp/external/fastjet/Error.$(ObjSuf): \
	external/fastjet/Error.$(SrcSuf) \
	external/fastjet/Error.hh \
	external/fastjet/config.h
tmp/external/fastjet/RangeDefinition.$(ObjSuf): \
	external/fastjet/RangeDefinition.$(SrcSuf) \
	external/fastjet/RangeDefinition.hh
tmp/external/fastjet/ClusterSequence_N2.$(ObjSuf): \
	external/fastjet/ClusterSequence_N2.$(SrcSuf)
tmp/external/fastjet/ClusterSequenceStructure.$(ObjSuf): \
	external/fastjet/ClusterSequenceStructure.$(SrcSuf) \
	external/fastjet/ClusterSequenceStructure.hh \
	external/fastjet/Error.hh \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/ClusterSequenceAreaBase.hh
tmp/external/fastjet/Dnn4piCylinder.$(ObjSuf): \
	external/fastjet/Dnn4piCylinder.$(SrcSuf) \
	external/fastjet/internal/Dnn4piCylinder.hh
tmp/external/fastjet/ClusterSequence1GhostPassiveArea.$(ObjSuf): \
	external/fastjet/ClusterSequence1GhostPassiveArea.$(SrcSuf) \
	external/fastjet/ClusterSequence1GhostPassiveArea.hh
tmp/external/fastjet/MinHeap.$(ObjSuf): \
	external/fastjet/MinHeap.$(SrcSuf) \
	external/fastjet/internal/MinHeap.hh
tmp/external/fastjet/PseudoJet.$(ObjSuf): \
	external/fastjet/PseudoJet.$(SrcSuf) \
	external/fastjet/Error.hh \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/ClusterSequenceAreaBase.hh \
	external/fastjet/CompositeJetStructure.hh
tmp/external/fastjet/Dnn3piCylinder.$(ObjSuf): \
	external/fastjet/Dnn3piCylinder.$(SrcSuf) \
	external/fastjet/internal/Dnn3piCylinder.hh
tmp/external/fastjet/ClusterSequence_CP2DChan.$(ObjSuf): \
	external/fastjet/ClusterSequence_CP2DChan.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/internal/ClosestPair2D.hh
tmp/external/fastjet/ClusterSequenceVoronoiArea.$(ObjSuf): \
	external/fastjet/ClusterSequenceVoronoiArea.$(SrcSuf) \
	external/fastjet/ClusterSequenceVoronoiArea.hh \
	external/fastjet/internal/Voronoi.hh
tmp/external/fastjet/ClosestPair2D.$(ObjSuf): \
	external/fastjet/ClosestPair2D.$(SrcSuf) \
	external/fastjet/internal/ClosestPair2D.hh
tmp/external/fastjet/FunctionOfPseudoJet.$(ObjSuf): \
	external/fastjet/FunctionOfPseudoJet.$(SrcSuf)
tmp/external/fastjet/AreaDefinition.$(ObjSuf): \
	external/fastjet/AreaDefinition.$(SrcSuf) \
	external/fastjet/AreaDefinition.hh
tmp/external/fastjet/tools/CASubJetTagger.$(ObjSuf): \
	external/fastjet/tools/CASubJetTagger.$(SrcSuf)
tmp/external/fastjet/tools/TopTaggerBase.$(ObjSuf): \
	external/fastjet/tools/TopTaggerBase.$(SrcSuf)
tmp/external/fastjet/tools/BackgroundEstimatorBase.$(ObjSuf): \
	external/fastjet/tools/BackgroundEstimatorBase.$(SrcSuf) \
	external/fastjet/tools/BackgroundEstimatorBase.hh
tmp/external/fastjet/tools/Subtractor.$(ObjSuf): \
	external/fastjet/tools/Subtractor.$(SrcSuf) \
	external/fastjet/tools/Subtractor.hh
tmp/external/fastjet/tools/RestFrameNSubjettinessTagger.$(ObjSuf): \
	external/fastjet/tools/RestFrameNSubjettinessTagger.$(SrcSuf)
tmp/external/fastjet/tools/JetMedianBackgroundEstimator.$(ObjSuf): \
	external/fastjet/tools/JetMedianBackgroundEstimator.$(SrcSuf) \
	external/fastjet/tools/JetMedianBackgroundEstimator.hh
tmp/external/fastjet/tools/MassDropTagger.$(ObjSuf): \
	external/fastjet/tools/MassDropTagger.$(SrcSuf)
tmp/external/fastjet/tools/JHTopTagger.$(ObjSuf): \
	external/fastjet/tools/JHTopTagger.$(SrcSuf)
tmp/external/fastjet/tools/GridMedianBackgroundEstimator.$(ObjSuf): \
	external/fastjet/tools/GridMedianBackgroundEstimator.$(SrcSuf) \
	external/fastjet/tools/GridMedianBackgroundEstimator.hh
tmp/external/fastjet/tools/Filter.$(ObjSuf): \
	external/fastjet/tools/Filter.$(SrcSuf) \
	external/fastjet/tools/Filter.hh
tmp/external/fastjet/tools/Pruner.$(ObjSuf): \
	external/fastjet/tools/Pruner.$(SrcSuf) \
	external/fastjet/tools/Pruner.hh \
	external/fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh \
	external/fastjet/Selector.hh
tmp/external/fastjet/plugins/ATLASCone/ATLASConePlugin.$(ObjSuf): \
	external/fastjet/plugins/ATLASCone/ATLASConePlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh
tmp/external/fastjet/plugins/ATLASCone/JetConeFinderTool.$(ObjSuf): \
	external/fastjet/plugins/ATLASCone/JetConeFinderTool.$(SrcSuf)
tmp/external/fastjet/plugins/ATLASCone/Jet.$(ObjSuf): \
	external/fastjet/plugins/ATLASCone/Jet.$(SrcSuf)
tmp/external/fastjet/plugins/ATLASCone/JetSplitMergeTool.$(ObjSuf): \
	external/fastjet/plugins/ATLASCone/JetSplitMergeTool.$(SrcSuf)
tmp/external/fastjet/plugins/NestedDefs/NestedDefsPlugin.$(ObjSuf): \
	external/fastjet/plugins/NestedDefs/NestedDefsPlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh
tmp/external/fastjet/plugins/D0RunIICone/D0RunIIConePlugin.$(ObjSuf): \
	external/fastjet/plugins/D0RunIICone/D0RunIIConePlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/Error.hh
tmp/external/fastjet/plugins/TrackJet/TrackJetPlugin.$(ObjSuf): \
	external/fastjet/plugins/TrackJet/TrackJetPlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh
tmp/external/fastjet/plugins/CDFCones/CDFMidPointPlugin.$(ObjSuf): \
	external/fastjet/plugins/CDFCones/CDFMidPointPlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/Error.hh
tmp/external/fastjet/plugins/CDFCones/MidPointAlgorithm.$(ObjSuf): \
	external/fastjet/plugins/CDFCones/MidPointAlgorithm.$(SrcSuf)
tmp/external/fastjet/plugins/CDFCones/CDFJetCluPlugin.$(ObjSuf): \
	external/fastjet/plugins/CDFCones/CDFJetCluPlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh
tmp/external/fastjet/plugins/CDFCones/JetCluAlgorithm.$(ObjSuf): \
	external/fastjet/plugins/CDFCones/JetCluAlgorithm.$(SrcSuf)
tmp/external/fastjet/plugins/D0RunICone/D0RunIBaseConePlugin.$(ObjSuf): \
	external/fastjet/plugins/D0RunICone/D0RunIBaseConePlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/Error.hh
tmp/external/fastjet/plugins/SISCone/geom_2d.$(ObjSuf): \
	external/fastjet/plugins/SISCone/geom_2d.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/SISConePlugin.$(ObjSuf): \
	external/fastjet/plugins/SISCone/SISConePlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh
tmp/external/fastjet/plugins/SISCone/siscone.$(ObjSuf): \
	external/fastjet/plugins/SISCone/siscone.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/hash.$(ObjSuf): \
	external/fastjet/plugins/SISCone/hash.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/quadtree.$(ObjSuf): \
	external/fastjet/plugins/SISCone/quadtree.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/area.$(ObjSuf): \
	external/fastjet/plugins/SISCone/area.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/momentum.$(ObjSuf): \
	external/fastjet/plugins/SISCone/momentum.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/reference.$(ObjSuf): \
	external/fastjet/plugins/SISCone/reference.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/ranlux.$(ObjSuf): \
	external/fastjet/plugins/SISCone/ranlux.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/protocones.$(ObjSuf): \
	external/fastjet/plugins/SISCone/protocones.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/split_merge.$(ObjSuf): \
	external/fastjet/plugins/SISCone/split_merge.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/vicinity.$(ObjSuf): \
	external/fastjet/plugins/SISCone/vicinity.$(SrcSuf)
tmp/external/fastjet/plugins/SISCone/siscone_error.$(ObjSuf): \
	external/fastjet/plugins/SISCone/siscone_error.$(SrcSuf)
tmp/external/fastjet/plugins/CMSIterativeCone/CMSIterativeConePlugin.$(ObjSuf): \
	external/fastjet/plugins/CMSIterativeCone/CMSIterativeConePlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh
tmp/external/fastjet/plugins/GridJet/GridJetPlugin.$(ObjSuf): \
	external/fastjet/plugins/GridJet/GridJetPlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh
tmp/external/fastjet/plugins/Jade/JadePlugin.$(ObjSuf): \
	external/fastjet/plugins/Jade/JadePlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/NNH.hh
tmp/external/fastjet/plugins/EECambridge/EECambridgePlugin.$(ObjSuf): \
	external/fastjet/plugins/EECambridge/EECambridgePlugin.$(SrcSuf) \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/NNH.hh
DELPHES_OBJ =  \
	tmp/classes/DelphesHepMCReader.$(ObjSuf) \
	tmp/classes/DelphesLHEFReader.$(ObjSuf) \
	tmp/classes/DelphesFactory.$(ObjSuf) \
	tmp/classes/DelphesFormula.$(ObjSuf) \
	tmp/classes/DelphesSTDHEPReader.$(ObjSuf) \
	tmp/classes/DelphesStream.$(ObjSuf) \
	tmp/classes/DelphesPileUpWriter.$(ObjSuf) \
	tmp/classes/DelphesModule.$(ObjSuf) \
	tmp/classes/DelphesClasses.$(ObjSuf) \
	tmp/classes/DelphesPileUpReader.$(ObjSuf) \
	tmp/modules/LeptonDressing.$(ObjSuf) \
	tmp/modules/Efficiency.$(ObjSuf) \
	tmp/modules/UniqueObjectFinder.$(ObjSuf) \
	tmp/modules/ExampleModule.$(ObjSuf) \
	tmp/modules/ParticlePropagator.$(ObjSuf) \
	tmp/modules/EnergySmearing.$(ObjSuf) \
	tmp/modules/FastJetFinder.$(ObjSuf) \
	tmp/modules/StatusPidFilter.$(ObjSuf) \
	tmp/modules/ConstituentFilter.$(ObjSuf) \
	tmp/modules/EnergyScale.$(ObjSuf) \
	tmp/modules/TauTagging.$(ObjSuf) \
	tmp/modules/Merger.$(ObjSuf) \
	tmp/modules/BTagging.$(ObjSuf) \
	tmp/modules/TreeWriter.$(ObjSuf) \
	tmp/modules/Weighter.$(ObjSuf) \
	tmp/modules/Delphes.$(ObjSuf) \
	tmp/modules/Calorimeter.$(ObjSuf) \
	tmp/modules/Isolation.$(ObjSuf) \
	tmp/modules/PileUpMerger.$(ObjSuf) \
	tmp/modules/TrackPileUpSubtractor.$(ObjSuf) \
	tmp/modules/Cloner.$(ObjSuf) \
	tmp/modules/JetPileUpSubtractor.$(ObjSuf) \
	tmp/modules/MomentumSmearing.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootFilter.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootConfReader.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootTreeWriter.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootResult.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootTreeBranch.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootTreeReader.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootTask.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootProgressBar.$(ObjSuf) \
	tmp/external/ExRootAnalysis/ExRootUtilities.$(ObjSuf) \
	tmp/external/fastjet/Dnn2piCylinder.$(ObjSuf) \
	tmp/external/fastjet/GhostedAreaSpec.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequenceActiveArea.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequence_Delaunay.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequenceArea.$(ObjSuf) \
	tmp/external/fastjet/Voronoi.$(ObjSuf) \
	tmp/external/fastjet/Selector.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequence_TiledN2.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequenceActiveAreaExplicitGhosts.$(ObjSuf) \
	tmp/external/fastjet/LimitedWarning.$(ObjSuf) \
	tmp/external/fastjet/CompositeJetStructure.$(ObjSuf) \
	tmp/external/fastjet/DnnPlane.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequence_DumbN3.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequencePassiveArea.$(ObjSuf) \
	tmp/external/fastjet/BasicRandom.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequenceAreaBase.$(ObjSuf) \
	tmp/external/fastjet/PseudoJetStructureBase.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequence.$(ObjSuf) \
	tmp/external/fastjet/JetDefinition.$(ObjSuf) \
	tmp/external/fastjet/Error.$(ObjSuf) \
	tmp/external/fastjet/RangeDefinition.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequence_N2.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequenceStructure.$(ObjSuf) \
	tmp/external/fastjet/Dnn4piCylinder.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequence1GhostPassiveArea.$(ObjSuf) \
	tmp/external/fastjet/MinHeap.$(ObjSuf) \
	tmp/external/fastjet/PseudoJet.$(ObjSuf) \
	tmp/external/fastjet/Dnn3piCylinder.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequence_CP2DChan.$(ObjSuf) \
	tmp/external/fastjet/ClusterSequenceVoronoiArea.$(ObjSuf) \
	tmp/external/fastjet/ClosestPair2D.$(ObjSuf) \
	tmp/external/fastjet/FunctionOfPseudoJet.$(ObjSuf) \
	tmp/external/fastjet/AreaDefinition.$(ObjSuf) \
	tmp/external/fastjet/tools/CASubJetTagger.$(ObjSuf) \
	tmp/external/fastjet/tools/TopTaggerBase.$(ObjSuf) \
	tmp/external/fastjet/tools/BackgroundEstimatorBase.$(ObjSuf) \
	tmp/external/fastjet/tools/Subtractor.$(ObjSuf) \
	tmp/external/fastjet/tools/RestFrameNSubjettinessTagger.$(ObjSuf) \
	tmp/external/fastjet/tools/JetMedianBackgroundEstimator.$(ObjSuf) \
	tmp/external/fastjet/tools/MassDropTagger.$(ObjSuf) \
	tmp/external/fastjet/tools/JHTopTagger.$(ObjSuf) \
	tmp/external/fastjet/tools/GridMedianBackgroundEstimator.$(ObjSuf) \
	tmp/external/fastjet/tools/Filter.$(ObjSuf) \
	tmp/external/fastjet/tools/Pruner.$(ObjSuf) \
	tmp/external/fastjet/plugins/ATLASCone/ATLASConePlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/ATLASCone/JetConeFinderTool.$(ObjSuf) \
	tmp/external/fastjet/plugins/ATLASCone/Jet.$(ObjSuf) \
	tmp/external/fastjet/plugins/ATLASCone/JetSplitMergeTool.$(ObjSuf) \
	tmp/external/fastjet/plugins/NestedDefs/NestedDefsPlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/D0RunIICone/D0RunIIConePlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/TrackJet/TrackJetPlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/CDFCones/CDFMidPointPlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/CDFCones/MidPointAlgorithm.$(ObjSuf) \
	tmp/external/fastjet/plugins/CDFCones/CDFJetCluPlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/CDFCones/JetCluAlgorithm.$(ObjSuf) \
	tmp/external/fastjet/plugins/D0RunICone/D0RunIBaseConePlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/geom_2d.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/SISConePlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/siscone.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/hash.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/quadtree.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/area.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/momentum.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/reference.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/ranlux.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/protocones.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/split_merge.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/vicinity.$(ObjSuf) \
	tmp/external/fastjet/plugins/SISCone/siscone_error.$(ObjSuf) \
	tmp/external/fastjet/plugins/CMSIterativeCone/CMSIterativeConePlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/GridJet/GridJetPlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/Jade/JadePlugin.$(ObjSuf) \
	tmp/external/fastjet/plugins/EECambridge/EECambridgePlugin.$(ObjSuf)

ifeq ($(HAS_PYTHIA8),true)
DELPHES_OBJ +=  \
	tmp/modules/PileUpMergerPythia8.$(ObjSuf)
endif
tmp/display/DelphesDisplay.$(ObjSuf): \
	display/DelphesDisplay.$(SrcSuf) \
	display/DelphesDisplay.h
tmp/display/DelphesCaloData.$(ObjSuf): \
	display/DelphesCaloData.$(SrcSuf) \
	display/DelphesCaloData.h
DISPLAY_OBJ =  \
	tmp/display/DelphesDisplay.$(ObjSuf) \
	tmp/display/DelphesCaloData.$(ObjSuf)

ifeq ($(HAS_PYTHIA8),true)
DISPLAY_OBJ +=  \
	
endif
tmp/external/tcl/tclObj.$(ObjSuf): \
	external/tcl/tclObj.c
tmp/external/tcl/tclUtil.$(ObjSuf): \
	external/tcl/tclUtil.c
tmp/external/tcl/tclAsync.$(ObjSuf): \
	external/tcl/tclAsync.c
tmp/external/tcl/tclPosixStr.$(ObjSuf): \
	external/tcl/tclPosixStr.c
tmp/external/tcl/tclCompile.$(ObjSuf): \
	external/tcl/tclCompile.c
tmp/external/tcl/tclHistory.$(ObjSuf): \
	external/tcl/tclHistory.c
tmp/external/tcl/tclResolve.$(ObjSuf): \
	external/tcl/tclResolve.c
tmp/external/tcl/tclParse.$(ObjSuf): \
	external/tcl/tclParse.c
tmp/external/tcl/tclVar.$(ObjSuf): \
	external/tcl/tclVar.c
tmp/external/tcl/tclIndexObj.$(ObjSuf): \
	external/tcl/tclIndexObj.c
tmp/external/tcl/tclCkalloc.$(ObjSuf): \
	external/tcl/tclCkalloc.c
tmp/external/tcl/tclListObj.$(ObjSuf): \
	external/tcl/tclListObj.c
tmp/external/tcl/tclHash.$(ObjSuf): \
	external/tcl/tclHash.c
tmp/external/tcl/tclCmdIL.$(ObjSuf): \
	external/tcl/tclCmdIL.c
tmp/external/tcl/tclStringObj.$(ObjSuf): \
	external/tcl/tclStringObj.c
tmp/external/tcl/tclAlloc.$(ObjSuf): \
	external/tcl/tclAlloc.c
tmp/external/tcl/tclCompExpr.$(ObjSuf): \
	external/tcl/tclCompExpr.c
tmp/external/tcl/tclLink.$(ObjSuf): \
	external/tcl/tclLink.c
tmp/external/tcl/tclCmdAH.$(ObjSuf): \
	external/tcl/tclCmdAH.c
tmp/external/tcl/panic.$(ObjSuf): \
	external/tcl/panic.c
tmp/external/tcl/tclBasic.$(ObjSuf): \
	external/tcl/tclBasic.c
tmp/external/tcl/tclPreserve.$(ObjSuf): \
	external/tcl/tclPreserve.c
tmp/external/tcl/tclGet.$(ObjSuf): \
	external/tcl/tclGet.c
tmp/external/tcl/tclNamesp.$(ObjSuf): \
	external/tcl/tclNamesp.c
tmp/external/tcl/tclProc.$(ObjSuf): \
	external/tcl/tclProc.c
tmp/external/tcl/tclExecute.$(ObjSuf): \
	external/tcl/tclExecute.c
tmp/external/tcl/tclCmdMZ.$(ObjSuf): \
	external/tcl/tclCmdMZ.c
TCL_OBJ =  \
	tmp/external/tcl/tclObj.$(ObjSuf) \
	tmp/external/tcl/tclUtil.$(ObjSuf) \
	tmp/external/tcl/tclAsync.$(ObjSuf) \
	tmp/external/tcl/tclPosixStr.$(ObjSuf) \
	tmp/external/tcl/tclCompile.$(ObjSuf) \
	tmp/external/tcl/tclHistory.$(ObjSuf) \
	tmp/external/tcl/tclResolve.$(ObjSuf) \
	tmp/external/tcl/tclParse.$(ObjSuf) \
	tmp/external/tcl/tclVar.$(ObjSuf) \
	tmp/external/tcl/tclIndexObj.$(ObjSuf) \
	tmp/external/tcl/tclCkalloc.$(ObjSuf) \
	tmp/external/tcl/tclListObj.$(ObjSuf) \
	tmp/external/tcl/tclHash.$(ObjSuf) \
	tmp/external/tcl/tclCmdIL.$(ObjSuf) \
	tmp/external/tcl/tclStringObj.$(ObjSuf) \
	tmp/external/tcl/tclAlloc.$(ObjSuf) \
	tmp/external/tcl/tclCompExpr.$(ObjSuf) \
	tmp/external/tcl/tclLink.$(ObjSuf) \
	tmp/external/tcl/tclCmdAH.$(ObjSuf) \
	tmp/external/tcl/panic.$(ObjSuf) \
	tmp/external/tcl/tclBasic.$(ObjSuf) \
	tmp/external/tcl/tclPreserve.$(ObjSuf) \
	tmp/external/tcl/tclGet.$(ObjSuf) \
	tmp/external/tcl/tclNamesp.$(ObjSuf) \
	tmp/external/tcl/tclProc.$(ObjSuf) \
	tmp/external/tcl/tclExecute.$(ObjSuf) \
	tmp/external/tcl/tclCmdMZ.$(ObjSuf)

external/fastjet/internal/ClosestPair2D.hh: \
	external/fastjet/internal/ClosestPair2DBase.hh \
	external/fastjet/internal/SearchTree.hh \
	external/fastjet/internal/MinHeap.hh
	@touch $@

external/fastjet/ClusterSequence.hh: \
	external/fastjet/internal/DynamicNearestNeighbours.hh \
	external/fastjet/PseudoJet.hh \
	external/fastjet/Error.hh \
	external/fastjet/JetDefinition.hh \
	external/fastjet/SharedPtr.hh \
	external/fastjet/LimitedWarning.hh \
	external/fastjet/FunctionOfPseudoJet.hh \
	external/fastjet/ClusterSequenceStructure.hh
	@touch $@

external/fastjet/internal/MinHeap.hh: \
	external/fastjet/internal/base.hh
	@touch $@

modules/EnergySmearing.h: \
	classes/DelphesModule.h
	@touch $@

modules/LeptonDressing.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/internal/Voronoi.hh: \
	external/fastjet/LimitedWarning.hh
	@touch $@

external/fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh: \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequenceAreaBase.hh \
	external/fastjet/GhostedAreaSpec.hh \
	external/fastjet/LimitedWarning.hh
	@touch $@

external/fastjet/JetDefinition.hh: \
	external/fastjet/internal/numconsts.hh \
	external/fastjet/PseudoJet.hh
	@touch $@

modules/ConstituentFilter.h: \
	classes/DelphesModule.h
	@touch $@

modules/Calorimeter.h: \
	classes/DelphesModule.h
	@touch $@

classes/DelphesModule.h: \
	external/ExRootAnalysis/ExRootTask.h
	@touch $@

modules/Isolation.h: \
	classes/DelphesModule.h
	@touch $@

modules/EnergyScale.h: \
	classes/DelphesModule.h
	@touch $@

modules/Merger.h: \
	classes/DelphesModule.h
	@touch $@

modules/ExampleModule.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/internal/Dnn2piCylinder.hh: \
	external/fastjet/internal/DynamicNearestNeighbours.hh \
	external/fastjet/internal/DnnPlane.hh \
	external/fastjet/internal/numconsts.hh
	@touch $@

external/fastjet/Selector.hh: \
	external/fastjet/PseudoJet.hh \
	external/fastjet/RangeDefinition.hh
	@touch $@

modules/JetPileUpSubtractor.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/Error.hh: \
	external/fastjet/internal/base.hh
	@touch $@

modules/Efficiency.h: \
	classes/DelphesModule.h
	@touch $@

modules/TrackPileUpSubtractor.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/tools/GridMedianBackgroundEstimator.hh: \
	external/fastjet/tools/BackgroundEstimatorBase.hh
	@touch $@

external/fastjet/internal/DnnPlane.hh: \
	external/fastjet/internal/Triangulation.hh \
	external/fastjet/internal/DynamicNearestNeighbours.hh
	@touch $@

external/fastjet/ClusterSequenceArea.hh: \
	external/fastjet/ClusterSequenceAreaBase.hh \
	external/fastjet/ClusterSequenceActiveArea.hh \
	external/fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh \
	external/fastjet/ClusterSequencePassiveArea.hh \
	external/fastjet/ClusterSequenceVoronoiArea.hh \
	external/fastjet/AreaDefinition.hh
	@touch $@

external/fastjet/ClusterSequence1GhostPassiveArea.hh: \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequenceAreaBase.hh \
	external/fastjet/ClusterSequenceActiveArea.hh
	@touch $@

modules/PileUpMerger.h: \
	classes/DelphesModule.h
	@touch $@

modules/Cloner.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/PseudoJet.hh: \
	external/fastjet/internal/numconsts.hh \
	external/fastjet/internal/IsBase.hh \
	external/fastjet/SharedPtr.hh \
	external/fastjet/Error.hh \
	external/fastjet/PseudoJetStructureBase.hh
	@touch $@

external/fastjet/tools/Pruner.hh: \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/WrappedStructure.hh \
	external/fastjet/tools/Transformer.hh
	@touch $@

external/fastjet/version.hh: \
	external/fastjet/config.h
	@touch $@

modules/MomentumSmearing.h: \
	classes/DelphesModule.h
	@touch $@

modules/TauTagging.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/GhostedAreaSpec.hh: \
	external/fastjet/PseudoJet.hh \
	external/fastjet/internal/BasicRandom.hh \
	external/fastjet/Selector.hh \
	external/fastjet/LimitedWarning.hh
	@touch $@

external/fastjet/internal/Dnn4piCylinder.hh: \
	external/fastjet/internal/DynamicNearestNeighbours.hh \
	external/fastjet/internal/DnnPlane.hh \
	external/fastjet/internal/numconsts.hh
	@touch $@

modules/Delphes.h: \
	classes/DelphesModule.h
	@touch $@

modules/UniqueObjectFinder.h: \
	classes/DelphesModule.h
	@touch $@

modules/PileUpMergerPythia8.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/ClusterSequenceActiveArea.hh: \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequenceAreaBase.hh \
	external/fastjet/ClusterSequenceActiveAreaExplicitGhosts.hh
	@touch $@

modules/ParticlePropagator.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/plugins/CDFCones/fastjet/CDFMidPointPlugin.hh: \
	external/fastjet/JetDefinition.hh
	@touch $@

external/fastjet/RangeDefinition.hh: \
	external/fastjet/PseudoJet.hh \
	external/fastjet/Error.hh \
	external/fastjet/LimitedWarning.hh
	@touch $@

external/fastjet/PseudoJetStructureBase.hh: \
	external/fastjet/internal/base.hh
	@touch $@

external/fastjet/ClusterSequenceAreaBase.hh: \
	external/fastjet/ClusterSequence.hh \
	external/fastjet/LimitedWarning.hh \
	external/fastjet/Selector.hh
	@touch $@

external/fastjet/ClusterSequenceVoronoiArea.hh: \
	external/fastjet/PseudoJet.hh \
	external/fastjet/AreaDefinition.hh \
	external/fastjet/ClusterSequenceAreaBase.hh
	@touch $@

modules/BTagging.h: \
	classes/DelphesModule.h
	@touch $@

modules/Weighter.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/internal/BasicRandom.hh: \
	external/fastjet/internal/base.hh
	@touch $@

external/fastjet/plugins/CDFCones/fastjet/CDFJetCluPlugin.hh: \
	external/fastjet/JetDefinition.hh \
	external/fastjet/PseudoJet.hh
	@touch $@

external/ExRootAnalysis/ExRootTask.h: \
	external/ExRootAnalysis/ExRootConfReader.h
	@touch $@

external/fastjet/tools/Subtractor.hh: \
	external/fastjet/tools/Transformer.hh \
	external/fastjet/tools/BackgroundEstimatorBase.hh
	@touch $@

external/fastjet/AreaDefinition.hh: \
	external/fastjet/GhostedAreaSpec.hh
	@touch $@

external/fastjet/internal/Dnn3piCylinder.hh: \
	external/fastjet/internal/DynamicNearestNeighbours.hh \
	external/fastjet/internal/DnnPlane.hh \
	external/fastjet/internal/numconsts.hh
	@touch $@

modules/TreeWriter.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/ClusterSequenceStructure.hh: \
	external/fastjet/internal/base.hh \
	external/fastjet/SharedPtr.hh \
	external/fastjet/PseudoJetStructureBase.hh
	@touch $@

modules/StatusPidFilter.h: \
	classes/DelphesModule.h
	@touch $@

external/fastjet/config.h: \
	external/fastjet/config_win.h
	@touch $@

external/fastjet/LimitedWarning.hh: \
	external/fastjet/internal/base.hh
	@touch $@

classes/DelphesClasses.h: \
	classes/SortableObject.h
	@touch $@

external/fastjet/ClusterSequencePassiveArea.hh: \
	external/fastjet/PseudoJet.hh \
	external/fastjet/ClusterSequence1GhostPassiveArea.hh
	@touch $@

modules/FastJetFinder.h: \
	classes/DelphesModule.h
	@touch $@



###

all: $(DELPHES) $(EXECUTABLE)

display: $(DISPLAY)

$(DELPHES): $(DELPHES_DICT_OBJ) $(DELPHES_OBJ) $(TCL_OBJ)
	@mkdir -p $(@D)
	@echo ">> Building $@"
ifeq ($(ARCH),aix5)
	@$(MAKESHARED) $(OutPutOpt) $@ $(DELPHES_LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	@$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(DELPHES_LIBS)
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
	@ln -sf $@ $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
	@bindexplib $* $^ > $*.def
	@lib -nologo -MACHINE:IX86 $^ -def:$*.def $(OutPutOpt)$(DELPHESLIB)
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(DELPHES_LIBS) $(OutPutOpt)$@
	@$(MT_DLL)
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(DELPHES_LIBS)
	@$(MT_DLL)
endif
endif
endif

$(DISPLAY): $(DELPHES_DICT_OBJ) $(DISPLAY_DICT_OBJ) $(DELPHES_OBJ) $(DISPLAY_OBJ) $(TCL_OBJ)
	@mkdir -p $(@D)
	@echo ">> Building $@"
ifeq ($(ARCH),aix5)
	@$(MAKESHARED) $(OutPutOpt) $@ $(DISPLAY_LIBS) -p 0 $^
else
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	@$(LD) $(SOFLAGS)$@ $(LDFLAGS) $^ $(OutPutOpt) $@ $(DISPLAY_LIBS)
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
	@ln -sf $@ $(subst .$(DllSuf),.so,$@)
endif
endif
else
ifeq ($(PLATFORM),win32)
	@bindexplib $* $^ > $*.def
	@lib -nologo -MACHINE:IX86 $^ -def:$*.def $(OutPutOpt)$(DISPLAYLIB)
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(DISPLAY_LIBS) $(OutPutOpt)$@
	@$(MT_DLL)
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(DISPLAY_LIBS)
	@$(MT_DLL)
endif
endif
endif

clean:
	@rm -f $(DELPHES_DICT_OBJ) $(DISPLAY_DICT_OBJ) $(DELPHES_OBJ) $(DISPLAY_OBJ) $(TCL_OBJ) core
	@rm -rf tmp

distclean: clean
	@rm -f $(DELPHES) $(DELPHESLIB) $(DISPLAY) $(DISPLAYLIB) $(EXECUTABLE)

dist:
	@echo ">> Building $(DISTTAR)"
	@mkdir -p $(DISTDIR)
	@cp -a CREDITS README VERSION Makefile configure classes converters display doc examples external modules python readers $(DISTDIR)
	@find $(DISTDIR) -depth -name .\* -exec rm -rf {} \;
	@tar -czf $(DISTTAR) $(DISTDIR)
	@rm -rf $(DISTDIR)

###

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)

%Dict.$(SrcSuf):
	@mkdir -p $(@D)
	@echo ">> Generating $@"
	@rootcint -f $@ -c $(CXXFLAGS) $<
	@echo "#define private public" > $@.arch
	@echo "#define protected public" >> $@.arch
	@mv $@ $@.base
	@cat $@.arch $< $@.base > $@
	@rm $@.arch $@.base

$(DELPHES_OBJ): tmp/%.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(DISPLAY_OBJ): tmp/%.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(DELPHES_DICT_OBJ): %.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(DISPLAY_DICT_OBJ): %.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(TCL_OBJ): tmp/%.$(ObjSuf): %.c
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@gcc $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(EXECUTABLE_OBJ): tmp/%.$(ObjSuf): %.cpp
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(EXECUTABLE): %$(ExeSuf): $(DELPHES_DICT_OBJ) $(DELPHES_OBJ) $(TCL_OBJ)
	@echo ">> Building $@"
	@$(LD) $(LDFLAGS) $^ $(DELPHES_LIBS) $(OutPutOpt)$@

###


