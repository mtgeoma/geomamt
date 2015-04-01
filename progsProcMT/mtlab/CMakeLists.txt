CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# para esse if que muda o tipo de montagem padrão veja:
# http://www.cmake.org/Wiki/CMake_FAQ
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#Path to customs cmake modules
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMakeModules)

ADD_SUBDIRECTORY(src)
