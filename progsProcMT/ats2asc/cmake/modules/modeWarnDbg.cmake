################################################################################
# Modo Warndbg (Warning + Debug)					       #
# 									       #
# Copia as flags do modo 'Debug' e acrescenta todas as flags encontradas por #
# setWarningFlags e CMAKE_CXX_FLAGS_DEBUG em CMAKE_CXX_FLAGS_WARNDBG	       #
################################################################################

set(CMAKE_CXX_FLAGS_WARNDBG "${cxxWarningFlags} ${CMAKE_CXX_FLAGS_DEBUG}"
  CACHE STRING
  "Flags used by the C++ compiler during warndbg builds."
  FORCE)

set(CMAKE_C_FLAGS_WARNDBG "${CMAKE_C_FLAGS_DEBUG}"
  CACHE STRING
  "Flags used by the C compiler during warndbg builds."
  FORCE)

set(CMAKE_EXE_LINKER_FLAGS_WARNDBG
  "${CMAKE_EXE_LINKER_FLAGS_DEBUG}"
  CACHE STRING
  "Flags used for linking binaries during warndbg builds."
  FORCE)

set(CMAKE_SHARED_LINKER_FLAGS_WARNDBG
  "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}"
  CACHE STRING
  "Flags used by the shared libraries linker during warndbg builds."
  FORCE)

mark_as_advanced(
    CMAKE_CXX_FLAGS_WARNDBG
    CMAKE_C_FLAGS_WARNDBG
    CMAKE_EXE_LINKER_FLAGS_WARNDBG
    CMAKE_SHARED_LINKER_FLAGS_WARNDBG)
