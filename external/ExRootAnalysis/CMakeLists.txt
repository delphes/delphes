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
install(FILES ExRootTask.h ExRootConfReader.h DESTINATION include/ExRootAnalysis)
