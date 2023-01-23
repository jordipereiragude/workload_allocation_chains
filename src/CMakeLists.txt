cmake_minimum_required(VERSION 2.8)

cmake_policy(SET CMP0048 NEW)
project(wap VERSION 1.0)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(CMakeBuildTypes)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build configuration: " ${CMAKE_BUILD_TYPE}) 

set(CMAKE_CXX_FLAGS "-DIL_STD -g -ggdb -Wall -Wno-unknown-pragmas -fexceptions -std=c++17")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -ggdb -fsanitize=address -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -ggdb -fno-omit-frame-pointer")

find_package(cplex)
set(Boost_NO_SYSTEM_PATHS ON)
find_package(Boost 1.75.0 COMPONENTS program_options timer system chrono)

include_directories(SYSTEM ${CPLEX_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
include_directories(AFTER .)
include_directories(../3rd/fmt/include)
link_directories(${Boost_LIBRARY_DIRS})

ADD_CUSTOM_TARGET(
  version
  ${CMAKE_COMMAND} -D SRC=${CMAKE_SOURCE_DIR}/version.hpp.in -D DST=${CMAKE_SOURCE_DIR}/version.hpp -P ${CMAKE_SOURCE_DIR}/cmake/version.cmake
)

add_executable(ipsolve ipsolve.cpp ip.cpp instance.cpp util.cpp random.cpp options.cpp statistics.cpp solution.cpp hs.cpp reducedIP.cpp beamSearch.cpp lb.cpp dp.cpp sadp.cpp ap.cpp solver.cpp cbfs.cpp memory.cpp cbfsdfs.cpp)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  target_link_libraries(ipsolve ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} dl)
  target_compile_definitions(ipsolve PUBLIC "IS_MACOS")
else()
  target_link_libraries(ipsolve ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} rt stdc++fs dl)
endif()

add_executable(info info.cpp instance.cpp util.cpp random.cpp options.cpp statistics.cpp solution.cpp)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  target_link_libraries(info ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} dl)
  target_compile_definitions(info PUBLIC "IS_MACOS")
else()
  target_link_libraries(info ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} rt stdc++fs dl)
endif()

add_executable(generate generate.cpp util.cpp random.cpp)
if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  target_link_libraries(generate ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} dl)
  target_compile_definitions(generate PUBLIC "IS_MACOS")
else()
  target_link_libraries(generate ${CPLEX_LIBRARIES} ${Boost_LIBRARIES} rt stdc++fs dl)
endif()

add_dependencies(ipsolve version)
