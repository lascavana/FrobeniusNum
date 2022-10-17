cmake_minimum_required(VERSION 3.3)
project(ls_reduction)

option(CXX "enable C++ compilation" ON)
### c++11 standards
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -std=c++11 -I/usr/local/include")

if(CXX)
  enable_language(CXX)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Debug")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

### NTL ###
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(GMP REQUIRED)
find_package(NTL REQUIRED)

### GUROBI ##
find_package(GUROBI REQUIRED)
include_directories(${GUROBI_INCLUDE_DIRS})


add_executable(ls_reduction
   src/main.cpp)

if(CXX)
  set(CMAKE_CXX_STANDARD 11)
  target_link_libraries(${CMAKE_PROJECT_NAME} optimized ${GUROBI_CXX_LIBRARY}
                        debug ${GUROBI_CXX_DEBUG_LIBRARY} ${NTL_LIBRARY} ${GMP_LIBRARIES} ${M_LIB})
endif()

target_link_libraries(${CMAKE_PROJECT_NAME} ${NTL_LIBRARY} ${GMP_LIBRARIES} ${M_LIB} ${GUROBI_LIBRARY})

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
  include(FeatureSummary)
  feature_summary(WHAT ALL)
endif()
