#!/usr/bin/env tclsh

set prefix "tmp/"
set suffix " \\\n\t"

set srcSuf {.$(SrcSuf)}
set objSuf {.$(ObjSuf)}
set pcmSuf {$(PcmSuf)}
set exeSuf {$(ExeSuf)}

proc dependencies {fileName firstLine {force 1} {command {}}} {
  global suffix headerFiles sourceFiles

  if {[info exists sourceFiles($fileName)]} return

  set sourceFiles($fileName) 1

  set list {}
  set fid [open $fileName]
  while {! [eof $fid]} {
    set line [gets $fid]
    if [regexp -- {^\s*#include\s*"((\w+/)+\w+\.(h|hh))"} $line] {
      set elements [split $line {"}]
      set file [lindex $elements 1]
      if [file exists $file] {
        lappend list $file
        set headerFiles($file) 1
      } elseif [file exists external/$file] {
        lappend list external/$file
        set headerFiles(external/$file) 1
      }
    }
  }

  if {[llength $list] > 0} {
    puts -nonewline $firstLine
    foreach file $list {puts -nonewline $suffix$file}
    if {$command != {}} {
      puts {}
      puts $command
    }
    puts {}
  } elseif {$force} {
    puts -nonewline $firstLine
    if {$command != {}} {
      puts {}
      puts $command
    }
    puts {}
  }

  close $fid
}

proc dictDeps {dictPrefix args} {

  global prefix suffix srcSuf objSuf pcmSuf

  set dict [lsort [eval glob -nocomplain $args]]

  set dictSrcFiles {}
  set dictObjFiles {}

  foreach fileName $dict {
    regsub {LinkDef\.h} $fileName {Dict} dictName
    set dictName $prefix$dictName

    lappend dictSrcFiles $dictName$srcSuf
    lappend dictObjFiles $dictName$objSuf
    lappend dictPcmFiles [file tail $dictName$pcmSuf]

    dependencies $fileName "$dictName$srcSuf:$suffix$fileName"

    puts -nonewline [file tail $dictName$pcmSuf]:$suffix
    puts -nonewline $dictName$pcmSuf$suffix
    puts -nonewline $dictName$srcSuf
    puts {}
  }

  puts -nonewline "${dictPrefix}_OBJ += $suffix"
  puts [join $dictObjFiles $suffix]
  puts {}

  puts -nonewline "${dictPrefix}_PCM += $suffix"
  puts [join $dictPcmFiles $suffix]
  puts {}
}

proc sourceDeps {srcPrefix args} {

  global prefix suffix srcSuf objSuf

  set source [lsort [eval glob -nocomplain $args]]

  set srcObjFiles {}
  set srcObjFilesFastJet {}
  set srcObjFilesPythia8 {}

  foreach fileName $source {
    regsub {\.cc} $fileName {} srcName
    set srcObjName $prefix$srcName

    if {$fileName == "modules/PileUpMergerPythia8.cc"} {
      lappend srcObjFilesPythia8 $srcObjName$objSuf
    } elseif {([string match {modules/FastJet*.cc} $fileName] || [string match {modules/RunPUPPI.cc} $fileName]) && $srcPrefix != {FASTJET}} {
      continue
    } else {
      lappend srcObjFiles $srcObjName$objSuf
    }

    dependencies $fileName "$srcObjName$objSuf:$suffix$srcName$srcSuf"
  }

  puts -nonewline "${srcPrefix}_OBJ += $suffix"
  puts [join $srcObjFiles $suffix]
  puts {}

  puts {ifeq ($(HAS_PYTHIA8),true)}
  puts -nonewline "${srcPrefix}_OBJ += $suffix"
  puts [join $srcObjFilesPythia8 $suffix]
  puts {endif}
  puts {}
}

proc tclDeps {} {

  global prefix suffix srcSuf objSuf

  set source [lsort [glob -nocomplain {external/tcl/*.c}]]

  set srcObjFiles {}

  foreach fileName $source {
    if {$fileName == "tcl/tclc.c" || $fileName == "tcl/tcl.c"} continue

    regsub {\.c} $fileName {} srcName
    set srcObjName $prefix$srcName

    lappend srcObjFiles $srcObjName$objSuf

    dependencies $fileName "$srcObjName$objSuf:$suffix$fileName"
  }

  puts -nonewline "TCL_OBJ += $suffix"
  puts [join $srcObjFiles $suffix]
  puts {}
}

proc executableDeps {args} {

  global prefix suffix objSuf exeSuf

  set executable [lsort [eval glob -nocomplain $args]]

  set exeFiles {}

  foreach fileName $executable {
    regsub {\.cpp} $fileName {} exeObjName
    set exeObjName $prefix$exeObjName
    set exeName [file tail $exeObjName]

    lappend exeFiles $exeName$exeSuf
    lappend exeObjFiles $exeObjName$objSuf

    puts "$exeName$exeSuf:$suffix$exeObjName$objSuf"
    puts {}

    dependencies $fileName "$exeObjName$objSuf:$suffix$fileName"
  }

  if [info exists exeFiles] {
    puts -nonewline "EXECUTABLE += $suffix"
    puts [join $exeFiles $suffix]
    puts {}
  }
  if [info exists exeObjFiles] {
    puts -nonewline "EXECUTABLE_OBJ += $suffix"
    puts [join $exeObjFiles $suffix]
    puts {}
  }
}

proc headerDeps {} {
  global suffix headerFiles

  foreach fileName [array names headerFiles] {
    dependencies $fileName "$fileName:" 0 "\t@touch \$@"
  }
}

puts {
#
# Makefile for ExRootAnalysis
#
# Author: P. Demin - UCL, Louvain-la-Neuve
#
# multi-platform configuration is taken from ROOT (root/test/Makefile.arch)
#

include doc/Makefile.arch

ROOT_MAJOR := $(shell $(RC) --version | cut -d'.' -f1)

SrcSuf = cc
PcmSuf = _rdict.pcm

CXXFLAGS += $(ROOTCFLAGS) -Wno-write-strings -D_FILE_OFFSET_BITS=64 -DDROP_CGAL -I. -Iexternal -Iexternal/tcl
DELPHES_LIBS = $(shell $(RC) --libs) -lEG $(SYSLIBS)
DISPLAY_LIBS = $(shell $(RC) --evelibs) -lGuiHtml  $(SYSLIBS)

ifneq ($(CMSSW_FWLITE_INCLUDE_PATH),)
HAS_CMSSW = true
CXXFLAGS += -std=c++0x -I$(subst :, -I,$(CMSSW_FWLITE_INCLUDE_PATH))
OPT_LIBS += -L$(subst include,lib,$(subst :, -L,$(CMSSW_FWLITE_INCLUDE_PATH)))
ifneq ($(CMSSW_RELEASE_BASE),)
CXXFLAGS += -I$(CMSSW_RELEASE_BASE)/src
endif
ifneq ($(LD_LIBRARY_PATH),)
OPT_LIBS += -L$(subst include,lib,$(subst :, -L,$(LD_LIBRARY_PATH)))
endif
OPT_LIBS += -lGenVector -lFWCoreFWLite -lDataFormatsFWLite -lDataFormatsCommon -lDataFormatsPatCandidates -lDataFormatsLuminosity -lSimDataFormatsGeneratorProducts -lCommonToolsUtils -lDataFormatsCommon
endif

ifneq ($(PROMC),)
HAS_PROMC = true
CXXFLAGS += -I$(PROMC)/include -I$(PROMC)/src
OPT_LIBS += -L$(PROMC)/lib -lpromc -lprotoc -lprotobuf -lprotobuf-lite -lcbook -lz
endif

ifneq ($(PYTHIA8),)
#HAS_PYTHIA8 = true
CXXFLAGS += -I$(PYTHIA8)/include
CXXFLAGS += -I$(PYTHIA8)/include/Pythia8
OPT_LIBS += -L$(PYTHIA8)/lib -lpythia8 -ldl
endif

DELPHES_LIBS += $(OPT_LIBS)
DISPLAY_LIBS += $(OPT_LIBS)

###

NOFASTJET = libDelphesNoFastJet.$(DllSuf)
NOFASTJETLIB = libDelphesNoFastJet.lib

DELPHES = libDelphes.$(DllSuf)
DELPHESLIB = libDelphes.lib

DISPLAY = libDelphesDisplay.$(DllSuf)
DISPLAYLIB = libDelphesDisplay.lib

VERSION = $(shell cat VERSION)
DISTDIR = Delphes-$(VERSION)
DISTTAR = $(DISTDIR).tar.gz

all:

}

executableDeps {converters/*.cpp} {examples/*.cpp}

executableDeps {readers/DelphesHepMC.cpp} {readers/DelphesLHEF.cpp} {readers/DelphesSTDHEP.cpp}

puts {ifeq ($(HAS_CMSSW),true)}
executableDeps {readers/DelphesCMSFWLite.cpp}
puts {endif}
puts {}

puts {ifeq ($(HAS_PROMC),true)}
executableDeps {readers/DelphesProMC.cpp}
puts {endif}
puts {}

puts {ifeq ($(HAS_PYTHIA8),true)}
executableDeps {readers/DelphesPythia8.cpp}
dictDeps {DELPHES_DICT} {modules/Pythia8LinkDef.h}
puts {endif}
puts {}

dictDeps {DELPHES_DICT} {classes/ClassesLinkDef.h} {modules/ModulesLinkDef.h} {external/ExRootAnalysis/ExRootAnalysisLinkDef.h}

dictDeps {FASTJET_DICT} {modules/FastJetLinkDef.h}

dictDeps {DISPLAY_DICT} {display/DisplayLinkDef.h}

sourceDeps {DELPHES} {classes/*.cc} {modules/*.cc} {external/ExRootAnalysis/*.cc} {external/Hector/*.cc}

sourceDeps {FASTJET} {modules/FastJet*.cc} {modules/RunPUPPI.cc} {external/PUPPI/*.cc} {external/fastjet/*.cc} {external/fastjet/tools/*.cc} {external/fastjet/plugins/*/*.cc} {external/fastjet/contribs/*/*.cc} 

sourceDeps {DISPLAY} {display/*.cc}

tclDeps

headerDeps

puts {

###

ifeq ($(ROOT_MAJOR),6)
all: $(NOFASTJET) $(DELPHES) $(DELPHES_DICT_PCM) $(FASTJET_DICT_PCM) $(EXECUTABLE)
display: $(DISPLAY) $(DISPLAY_DICT_PCM)
else
all: $(NOFASTJET) $(DELPHES) $(EXECUTABLE)
display: $(DISPLAY)
endif

$(NOFASTJET): $(DELPHES_DICT_OBJ) $(DELPHES_OBJ) $(TCL_OBJ)
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
	@lib -nologo -MACHINE:IX86 $^ -def:$*.def $(OutPutOpt)$(NOFASTJETLIB)
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $*.exp $(DELPHES_LIBS) $(OutPutOpt)$@
	@$(MT_DLL)
else
	@$(LD) $(SOFLAGS) $(LDFLAGS) $^ $(OutPutOpt) $@ $(DELPHES_LIBS)
	@$(MT_DLL)
endif
endif
endif

$(DELPHES): $(DELPHES_DICT_OBJ) $(FASTJET_DICT_OBJ) $(DELPHES_OBJ) $(FASTJET_OBJ) $(TCL_OBJ)
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

$(DISPLAY): $(DELPHES_DICT_OBJ) $(FASTJET_DICT_OBJ) $(DISPLAY_DICT_OBJ) $(DELPHES_OBJ) $(FASTJET_OBJ) $(DISPLAY_OBJ) $(TCL_OBJ)
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
	@rm -f $(DELPHES_DICT_OBJ) $(DISPLAY_DICT_OBJ) $(DELPHES_OBJ) $(FASTJET_OBJ) $(DISPLAY_OBJ) $(TCL_OBJ) core
	@rm -rf tmp

distclean: clean
	@rm -f $(NOFASTJET) $(NOFASTJETLIB) $(DELPHES) $(DELPHESLIB) $(DELPHES_DICT_PCM) $(FASTJET_DICT_PCM) $(DISPLAY) $(DISPLAYLIB) $(DISPLAY_DICT_PCM) $(EXECUTABLE)

dist:
	@echo ">> Building $(DISTTAR)"
	@mkdir -p $(DISTDIR)
	@cp -a CHANGELOG CMakeLists.txt COPYING CREDITS DelphesEnv.sh README README_4LHCb VERSION Makefile MinBias.pileup configure cards classes converters display doc examples external modules python readers $(DISTDIR)
	@find $(DISTDIR) -depth -name .\* -exec rm -rf {} \;
	@tar -czf $(DISTTAR) $(DISTDIR)
	@rm -rf $(DISTDIR)

###

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf) $(PcmSuf)

%Dict.$(SrcSuf):
	@mkdir -p $(@D)
	@echo ">> Generating $@"
	@rootcint -f $@ -c -Iexternal $<
	@echo "#define private public" > $@.arch
	@echo "#define protected public" >> $@.arch
	@mv $@ $@.base
	@cat $@.arch $< $@.base > $@
	@rm $@.arch $@.base

%Dict$(PcmSuf):
	@echo ">> Copying $@"
	@cp $< $@

$(DELPHES_OBJ): tmp/%.$(ObjSuf): %.$(SrcSuf)
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(FASTJET_OBJ): tmp/%.$(ObjSuf): %.$(SrcSuf)
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

$(FASTJET_DICT_OBJ): %.$(ObjSuf): %.$(SrcSuf)
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
	@$(CC) $(patsubst -std=%,,$(CXXFLAGS)) -c $< $(OutPutOpt)$@

$(EXECUTABLE_OBJ): tmp/%.$(ObjSuf): %.cpp
	@mkdir -p $(@D)
	@echo ">> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< $(OutPutOpt)$@

$(EXECUTABLE): %$(ExeSuf): $(DELPHES_DICT_OBJ) $(FASTJET_DICT_OBJ) $(DELPHES_OBJ) $(FASTJET_OBJ) $(TCL_OBJ)
	@echo ">> Building $@"
	@$(LD) $(LDFLAGS) $^ $(DELPHES_LIBS) $(OutPutOpt)$@

###

}
