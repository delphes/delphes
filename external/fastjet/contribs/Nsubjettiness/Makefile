# If you are using this Makefile standalone and fastjet-config is not
# in your path, edit this line to specify the full path
FASTJETCONFIG=fastjet-config
PREFIX=`$(FASTJETCONFIG) --prefix`
CXX=g++
CXXFLAGS= -O3 -Wall -Woverloaded-virtual -g -Wunused-parameter
install_script = $(SHELL) ../utils/install-sh
check_script = ../utils/check.sh

# global contrib-wide Makefile include may override some of the above
# variables (leading "-" means don't give an error if you can't find
# the file)
-include ../.Makefile.inc

#------------------------------------------------------------------------
# things that are specific to this contrib
NAME=Nsubjettiness
SRCS=Nsubjettiness.cc Njettiness.cc NjettinessPlugin.cc MeasureFunction.cc AxesFinder.cc WinnerTakeAllRecombiner.cc NjettinessDefinition.cc
EXAMPLES=example_basic_usage example_advanced_usage example_v1p0p3
INSTALLED_HEADERS=Nsubjettiness.hh Njettiness.hh NjettinessPlugin.hh MeasureFunction.hh AxesFinder.hh WinnerTakeAllRecombiner.hh NjettinessDefinition.hh
#------------------------------------------------------------------------

CXXFLAGS+= $(shell $(FASTJETCONFIG) --cxxflags)
LDFLAGS += -lm $(shell $(FASTJETCONFIG) --libs)

OBJS  = $(SRCS:.cc=.o)
EXAMPLES_SRCS  = $(EXAMPLES:=.cc)

install_HEADER  = $(install_script) -c -m 644
install_LIB     = $(install_script) -c -m 644
install_DIR     = $(install_script) -d
install_DATA    = $(install_script) -c -m 644
install_PROGRAM = $(install_script) -c -s
install_SCRIPT  = $(install_script) -c

.PHONY: clean distclean examples check install

# compilation of the code (default target)
all: lib$(NAME).a

lib$(NAME).a: $(OBJS)
	ar cru lib$(NAME).a $(OBJS)
	ranlib lib$(NAME).a

# building the examples
examples: $(EXAMPLES)

# the following construct alloews to build each of the examples listed
# in $EXAMPLES automatically
$(EXAMPLES): % : %.o all
	$(CXX) -o $@ $< -L. -l$(NAME) $(LDFLAGS)

# check that everything went fine
check: examples
	@for prog in $(EXAMPLES); do\
	  $(check_script) $${prog} ../data/single-event.dat || exit 1; \
	done
	@echo "All tests successful"

# cleaning the directory
clean:
	rm -f *~ *.o *.a

distclean: clean
	rm -f lib$(NAME).a $(EXAMPLES)

# install things in PREFIX/...
install: all
	$(install_DIR) $(PREFIX)/include/fastjet/contrib
	for header in $(INSTALLED_HEADERS); do\
	  $(install_HEADER) $$header $(PREFIX)/include/fastjet/contrib/;\
	done
	$(install_DIR) $(PREFIX)/lib
	$(install_LIB) lib$(NAME).a $(PREFIX)/lib

depend:
	makedepend -Y --   -- $(SRCS) $(EXAMPLES_SRCS)
# DO NOT DELETE

Nsubjettiness.o: Nsubjettiness.hh Njettiness.hh MeasureFunction.hh
Nsubjettiness.o: AxesFinder.hh WinnerTakeAllRecombiner.hh
Nsubjettiness.o: NjettinessDefinition.hh
Njettiness.o: Njettiness.hh MeasureFunction.hh AxesFinder.hh
Njettiness.o: WinnerTakeAllRecombiner.hh NjettinessDefinition.hh
NjettinessPlugin.o: NjettinessPlugin.hh Njettiness.hh MeasureFunction.hh
NjettinessPlugin.o: AxesFinder.hh WinnerTakeAllRecombiner.hh
NjettinessPlugin.o: NjettinessDefinition.hh
MeasureFunction.o: MeasureFunction.hh
AxesFinder.o: AxesFinder.hh WinnerTakeAllRecombiner.hh MeasureFunction.hh
WinnerTakeAllRecombiner.o: WinnerTakeAllRecombiner.hh
NjettinessDefinition.o: NjettinessDefinition.hh MeasureFunction.hh
NjettinessDefinition.o: AxesFinder.hh WinnerTakeAllRecombiner.hh
example_basic_usage.o: Nsubjettiness.hh Njettiness.hh MeasureFunction.hh
example_basic_usage.o: AxesFinder.hh WinnerTakeAllRecombiner.hh
example_basic_usage.o: NjettinessDefinition.hh NjettinessPlugin.hh
example_advanced_usage.o: Nsubjettiness.hh Njettiness.hh MeasureFunction.hh
example_advanced_usage.o: AxesFinder.hh WinnerTakeAllRecombiner.hh
example_advanced_usage.o: NjettinessDefinition.hh NjettinessPlugin.hh
example_v1p0p3.o: Nsubjettiness.hh Njettiness.hh MeasureFunction.hh
example_v1p0p3.o: AxesFinder.hh WinnerTakeAllRecombiner.hh
example_v1p0p3.o: NjettinessDefinition.hh NjettinessPlugin.hh
