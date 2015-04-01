CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#Path to customs cmake modules
SET(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# Montagem padrão: Release
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

#executáveis serão montados onde o make for executado
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# para evitar que o cmake use /usr/lib/libboost_filesystem-mt.a
SET(Boost_USE_MULTITHREADED FALSE)

ADD_SUBDIRECTORY(src)
