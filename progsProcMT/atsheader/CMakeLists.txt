CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# para esse if que muda o tipo de montagem padrão veja:
# http://www.cmake.org/Wiki/CMake_FAQ
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

PROJECT(atsheader)

SET(CMAKE_CXX_FLAGS "-pedantic -Wall -Wextra -Wconversion")
SET(CMAKE_CXX_FLAGS_DEBUG  "-g -DBZ_DEBUG")
SET(CMAKE_CXX_FLAGS_RELEASE  "-O3 -DNDEBUG")

SET(Boost_USE_STATIC_LIBS ON)
FIND_PACKAGE(Boost COMPONENTS date_time REQUIRED)
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIRS})
ADD_DEFINITIONS(-DBOOST_FILESYSTEM_NO_DEPRECATED)

SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

SET(APP_SRCS ${PROJECT_NAME}.cpp read_header.cpp write_header.cpp)

ADD_EXECUTABLE(${PROJECT_NAME} ${APP_SRCS})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES})
