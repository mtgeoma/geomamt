PROJECT(mtlab)

#FINDGMM
#http://websvn.kde.org/trunk/koffice/cmake/modules/FindGMM.cmake?view=markup
FIND_PACKAGE(GMM REQUIRED)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/lib)
ADD_LIBRARY(${PROJECT_NAME} import/import.cpp import/import_egbert.cpp
  import/import_jones invariants/erro_manual/base_xi-eta/bahr_invariants.cpp
  invariants/erro_manual/base_xi-eta/invariants_parameters.cpp
  invariants/erro_manual/base_xi-eta/wal_invariants.cpp
  invariants/dimensionalidade/inicializa_thresholds.cpp
  invariants/dimensionalidade/bahr_Q_dimensionalidade.cpp
  invariants/dimensionalidade/wal_dimensionalidade.cpp
  transfer_functions/transfer_functions.cpp)
