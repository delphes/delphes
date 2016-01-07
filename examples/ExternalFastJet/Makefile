# 'make' builds everything
# 'make clean' deletes everything except source files and Makefile
#
# You need to set ROOT_DIR, DELPHES_DIR and FASTJET_DIR for your project.

# FastJet directory
FASTJET_DIR = /opt/fastjet

# ROOT directory
ROOT_DIR = /opt/root

# Delphes directory
DELPHES_DIR = ../..

#

CXXFLAGS += $(shell $(FASTJET_DIR)/bin/fastjet-config --cxxflags)
CXXFLAGS += $(shell $(ROOT_DIR)/bin/root-config --cflags)
CXXFLAGS += -I$(DELPHES_DIR) -I$(DELPHES_DIR)/external

LDFLAGS += $(shell $(FASTJET_DIR)/bin/fastjet-config --libs)
LDFLAGS += $(shell $(ROOT_DIR)/bin/root-config --libs) -lEG
LDFLAGS += -Wl,-rpath,$(DELPHES_DIR) -L$(DELPHES_DIR) -lDelphesNoFastJet

all: ExternalFastJetBasic ExternalFastJetHepMC

clean:
	rm -f ExternalFastJetBasic ExternalFastJetHepMC

%: %.cpp
	$(CXX) $(CXXFLAGS) $^ $(LDFLAGS) -o $@
