# Define a variável GCC_VERSION com a versão do compilador gcc
function(gccVersion)

  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} -dumpversion
    OUTPUT_VARIABLE GCC_DUMPVERSION_)
  string(STRIP "${GCC_DUMPVERSION_}" GCC_DUMPVERSION_)
  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} --version
    COMMAND head -n1
    OUTPUT_VARIABLE GCC_PACKAGER_)
  string(REGEX REPLACE ".+ \\((.+)\\) .+" "\\1"
    GCC_PACKAGER_ "${GCC_PACKAGER_}")
  string(REGEX REPLACE "([^ ]+) .*" "\\1"
    GCC_PACKAGER_ID_ "${GCC_PACKAGER_}")
  if ("${GCC_PACKAGER_ID_}" STREQUAL "Ubuntu/Linaro")
    string(REGEX REPLACE "[^ ]+ (.*)" "\\1"
      GCC_VERSION_ "${GCC_PACKAGER_}")
  elseif ("${GCC_PACKAGER_}" STREQUAL "GCC")
    set(GCC_VERSION_ "${GCC_DUMPVERSION_}")
  endif()
  set(GCC_VERSION ${GCC_VERSION_} PARENT_SCOPE)
endfunction()
