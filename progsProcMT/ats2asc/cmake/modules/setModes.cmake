#########################################
# Estabelece outros modos de compilação #
#########################################

include(setWarningFlags)
include(modeWarnDbg)
include(modeWarnRel)

# Update the documentation string of CMAKE_BUILD_TYPE for GUIs
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}"
  CACHE STRING
  "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel WarnDbg WarnRel."
  FORCE)
