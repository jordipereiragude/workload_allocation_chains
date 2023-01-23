# Add new build types

## profiling run
message("* Adding build type \"profile\"")
set(CMAKE_CXX_FLAGS_PROFILE
    "-no-pie ${GCC_DEBUG_FLAGS} -pg"
    CACHE STRING "Flags used by the C++ compiler during profile builds."
    FORCE )
set(CMAKE_C_FLAGS_PROFILE
    "${GCC_DEBUG_FLAGS} --pg"
    CACHE STRING "Flags used by the C compiler during profile builds."
    FORCE )
set(CMAKE_EXE_LINKER_FLAGS_PROFILE
    ""
    CACHE STRING "Flags used for linking binaries during profile builds."
    FORCE )
set(CMAKE_SHARED_LINKER_FLAGS_PROFILE
    ""
    CACHE STRING "Flags used by the shared libraries linker during profile builds."
    FORCE )
mark_as_advanced(
    CMAKE_CXX_FLAGS_PROFILE
    CMAKE_C_FLAGS_PROFILE
    CMAKE_EXE_LINKER_FLAGS_PROFILE
    CMAKE_SHARED_LINKER_FLAGS_PROFILE )

## profiling w/ full optimization
message("* Adding build type \"optprofile\"")
set(CMAKE_CXX_FLAGS_OPTPROFILE
    "-no-pie -O3 -g -pg -DNDEBUG"
    CACHE STRING "Flags used by the C++ compiler during optimized profile builds."
    FORCE )
set(CMAKE_C_FLAGS_OPTPROFILE
    "-O3 -g -pg -DNDEBUG"
    CACHE STRING "Flags used by the C compiler during optimized profile builds."
    FORCE )
set(CMAKE_EXE_LINKER_FLAGS_OPTPROFILE
    ""
    CACHE STRING "Flags used for linking binaries during optimized profile builds."
    FORCE )
set(CMAKE_SHARED_LINKER_FLAGS_OPTPROFILE
    ""
    CACHE STRING "Flags used by the shared libraries linker during optimized profile builds."
    FORCE )
mark_as_advanced(
    CMAKE_CXX_FLAGS_OPTPROFILE
    CMAKE_C_FLAGS_OPTPROFILE
    CMAKE_EXE_LINKER_FLAGS_OPTPROFILE
    CMAKE_SHARED_LINKER_FLAGS_OPTPROFILE )

## profiling w/ simple optimization
message("* Adding build type \"soptprofile\"")
set(CMAKE_CXX_FLAGS_SOPTPROFILE
    "-no-pie -O -g -pg -DNDEBUG"
    CACHE STRING "Flags used by the C++ compiler during optimized profile builds."
    FORCE )
set(CMAKE_C_FLAGS_SOPTPROFILE
    "-O -g -pg -DNDEBUG"
    CACHE STRING "Flags used by the C compiler during optimized profile builds."
    FORCE )
set(CMAKE_EXE_LINKER_FLAGS_SOPTPROFILE
    ""
    CACHE STRING "Flags used for linking binaries during optimized profile builds."
    FORCE )
set(CMAKE_SHARED_LINKER_FLAGS_SOPTPROFILE
    ""
    CACHE STRING "Flags used by the shared libraries linker during optimized profile builds."
    FORCE )
mark_as_advanced(
    CMAKE_CXX_FLAGS_SOPTPROFILE
    CMAKE_C_FLAGS_SOPTPROFILE
    CMAKE_EXE_LINKER_FLAGS_SOPTPROFILE
    CMAKE_SHARED_LINKER_FLAGS_SOPTPROFILE )
