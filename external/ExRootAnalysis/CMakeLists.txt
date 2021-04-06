include_directories(
  ${CMAKE_SOURCE_DIR}/external
  ${ROOT_INCLUDE_DIRS}
)

file(GLOB sources *.cc)
file(GLOB headers *.h)

list(REMOVE_ITEM headers ${CMAKE_CURRENT_SOURCE_DIR}/ExRootAnalysisLinkDef.h)

# the macro invocation for ROOT6 ensures that the headers are relocatable
if (NOT ${ROOT_VERSION} VERSION_LESS "6.0.0")
  DELPHES_GENERATE_DICTIONARY(ExRootAnalysisDict
    ExRootAnalysis/ExRootClassifier.h
    ExRootAnalysis/ExRootConfReader.h
    ExRootAnalysis/ExRootFilter.h
    ExRootAnalysis/ExRootProgressBar.h
    ExRootAnalysis/ExRootResult.h
    ExRootAnalysis/ExRootTask.h
    ExRootAnalysis/ExRootTreeBranch.h
    ExRootAnalysis/ExRootTreeReader.h
    ExRootAnalysis/ExRootTreeWriter.h
    ExRootAnalysis/ExRootUtilities.h
    LINKDEF ExRootAnalysisLinkDef.h)
else()
  # for ROOT5 the above fails, keep the following as workaround
  DELPHES_GENERATE_DICTIONARY(ExRootAnalysisDict ${headers} LINKDEF ExRootAnalysisLinkDef.h)
endif()

add_library(ExRootAnalysis OBJECT ${sources} ExRootAnalysisDict.cxx)

# install headers
install(FILES ${headers} DESTINATION include/ExRootAnalysis)

# install pcms if they are created
if (NOT ${ROOT_VERSION} VERSION_LESS "6.0.0")
  install(FILES
      ${PROJECT_BINARY_DIR}/external/ExRootAnalysis/libExRootAnalysisDict_rdict.pcm
      DESTINATION lib)
endif()
