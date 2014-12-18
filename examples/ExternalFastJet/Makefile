# 'make' builds everything
# 'make clean' deletes everything except source files and Makefile
#
# You need to set ROOT_DIR, DELPHES_DIR and FASTJET_DIR for your project.

# ROOT directory
ROOT_DIR = /work/root/root_5.34.18

# Delphes directory
DELPHES_DIR = ../..

# FastJet directory
FASTJET_DIR = /work/FastJet/fastjet_3.1.0

CXXFLAGS += $(shell $(FASTJET_DIR)/bin/fastjet-config --cxxflags)
CXXFLAGS += -I$(shell $(ROOT_DIR)/bin/root-config --incdir)
CXXFLAGS += -I$(DELPHES_DIR) -I$(DELPHES_DIR)/external

LDFLAGS += $(shell $(FASTJET_DIR)/bin/fastjet-config --libs)
LDFLAGS += $(shell $(ROOT_DIR)/bin/root-config --libs) -lEG
LDFLAGS += -Wl,-rpath=$(DELPHES_DIR) -L$(DELPHES_DIR) -lDelphesNoFastJet

all: ExternalFastJetBasic ExternalFastJetHepMC

clean:
	rm -f ExternalFastJetBasic ExternalFastJetHepMC

%: %.cpp
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@
