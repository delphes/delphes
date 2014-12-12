#######################################
# Order of execution of various modules
#######################################

set ExecutionPath {
  TreeWriter
}

##################
# ROOT tree writer
##################

module TreeWriter TreeWriter {
# add Branch InputArray BranchName BranchClass
  add Branch Delphes/allParticles Particle GenParticle
}

