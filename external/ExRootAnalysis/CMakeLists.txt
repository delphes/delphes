include_directories(
  ${CMAKE_SOURCE_DIR}/external
  ${ROOT_INCLUDE_DIRS}
)

file(GLOB sources *.cc)
file(GLOB headers *.h)
list(REMOVE_ITEM headers ${CMAKE_CURRENT_SOURCE_DIR}/ExRootAnalysisLinkDef.h)

DELPHES_GENERATE_DICTIONARY(ExRootAnalysisDict ${headers} LINKDEF ExRootAnalysisLinkDef.h)

add_library(ExRootAnalysis OBJECT ${sources} ExRootAnalysisDict.cxx)

# install headers needed by public Delphes headers to include/
install(FILES ExRootTask.h ExRootConfReader.h ExRootTreeWriter.h ExRootTreeBranch.h
        DESTINATION include/ExRootAnalysis)

# install all LinkDef files into the same folder to ease user environment
install(FILES ExRootAnalysisLinkDef.h DESTINATION include)

# install pcms if they are created
if (${ROOT_VERSION} GREATER 6)
  install(FILES
      ${PROJECT_BINARY_DIR}/external/ExRootAnalysis/libExRootAnalysisDict_rdict.pcm
      DESTINATION lib)
endif()
