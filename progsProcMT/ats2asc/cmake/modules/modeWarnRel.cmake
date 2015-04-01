################################################################################
# Modo WarnRel (Warning + Release)					       #
# 									       #
# Copia as flags do modo 'Release' e acrescenta todas as flags encontradas por #
# setWarningFlags e CMAKE_CXX_FLAGS_RELEASE em CMAKE_CXX_FLAGS_WARNREL	       #
################################################################################

set(CMAKE_CXX_FLAGS_WARNREL "${cxxWarningFlags} ${CMAKE_CXX_FLAGS_RELEASE}"
  CACHE STRING
  "Flags used by the C++ compiler during warnrel builds."
  FORCE)

set(CMAKE_C_FLAGS_WARNREL "${CMAKE_C_FLAGS_RELEASE}"
  CACHE STRING
  "Flags used by the C compiler during warnrel builds."
  FORCE)

set(CMAKE_EXE_LINKER_FLAGS_WARNREL
  "${CMAKE_EXE_LINKER_FLAGS_RELEASE}"
  CACHE STRING
  "Flags used for linking binaries during warnrel builds."
  FORCE)

set(CMAKE_SHARED_LINKER_FLAGS_WARNREL
  "${CMAKE_SHARED_LINKER_FLAGS_RELEASE}"
  CACHE STRING
  "Flags used by the shared libraries linker during warnrel builds."
  FORCE)

mark_as_advanced(
    CMAKE_CXX_FLAGS_WARNREL
    CMAKE_C_FLAGS_WARNREL
    CMAKE_EXE_LINKER_FLAGS_WARNREL
    CMAKE_SHARED_LINKER_FLAGS_WARNREL)
