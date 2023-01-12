file(GLOB confDir *.tcl)
# copy *.tcl files into cards
install(FILES ${confDir} DESTINATION cards)

# Install FCC cards into subfolder, including scenarios subdirectory
install(DIRECTORY FCC DESTINATION cards
  FILES_MATCHING PATTERN "*.tcl")

# Install CMS_PhaseII cards into subfolder
install(DIRECTORY CMS_PhaseII DESTINATION cards
  FILES_MATCHING PATTERN "*.tcl")
