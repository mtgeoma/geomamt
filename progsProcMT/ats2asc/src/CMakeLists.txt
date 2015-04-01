#Werror em uma linha em separado para facilitar a sua colocação em um
#comentário desativando a flag
SET(CMAKE_CXX_FLAGS "-Werror")

INCLUDE(setModes)

ADD_SUBDIRECTORY(app)
