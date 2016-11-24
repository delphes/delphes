file(GLOB confDir *.tcl)
file(GLOB confFccDir FCC/*.tcl)
file(GLOB confCmsPhase2Dir CMS_PhaseII/*.tcl)

# copy *.tcl files into cards
install(FILES ${confDir} DESTINATION cards)
install(FILES ${confFccDir} DESTINATION cards/FCC)
install(FILES ${confCmsPhase2Dir} DESTINATION cards/CMS_PhaseII)
