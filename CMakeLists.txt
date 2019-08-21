project                     (diy2-vtk7)
cmake_minimum_required      (VERSION 2.8)

option                      (diy_thread     "Enable diy threading"                                  OFF)
option                      (WITH_VTK                                                               OFF)
option                      (WITH_TIMEINFO                                                          ON)
option                      (WITH_MPE       "Enable MPE logging"                                    OFF)
option                      (wrapped_mpi    "MPI compiler wrapper requires no further MPI libs"     OFF)
option                      (log            "Build DIY with logging"                                OFF)
option                      (profile        "Build DIY with profiling"                              OFF)
option                      (caliper        "Build DIY with caliper"                                OFF)



# Default to Release
if                          (NOT CMAKE_BUILD_TYPE)
    set                     (CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
    set_property            (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif                       (NOT CMAKE_BUILD_TYPE)


if                          (bgq)
    add_definitions         (-DBGQ)
endif                       (bgq)

# Logging
if                          (log)
    add_definitions         (-DDIY_USE_SPDLOG)
    find_path               (SPDLOG_INCLUDE_DIR     spdlog/spdlog.h)
    include_directories     (${SPDLOG_INCLUDE_DIR})
endif                       ()

# C++11
set                         (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if                          (profile)
    add_definitions         (-DDIY_PROFILE)
endif                       ()

if                          (caliper)
    add_definitions         (-DDIY_USE_CALIPER)

    find_package            (caliper)
    include_directories     (${caliper_INCLUDE_DIR})
    set                     (libraries  ${libraries} caliper caliper-mpi)
endif()

# OSX flags
if                          (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  add_definitions           (-DMAC_OSX)
  set                       (CMAKE_MACOSX_RPATH on)

    # --- following RPATH settings are for Sierra w/ Clang, hopefully they don't hurt other versions
    # ref: https://cmake.org/Wiki/CMake_RPATH_handling
    # use, i.e. don't skip, the full RPATH for the build tree
    set                     (CMAKE_SKIP_BUILD_RPATH            false)
    # when building, don't use the install RPATH already (but later on when installing)
    set                     (CMAKE_BUILD_WITH_INSTALL_RPATH    false)
    # set RPATH to install path
    set                     (CMAKE_INSTALL_RPATH               "${CMAKE_INSTALL_PREFIX}/lib")
    # add the automatically determined parts of the RPATH
    # which point to directories outside the build tree to the install RPATH
    set                     (CMAKE_INSTALL_RPATH_USE_LINK_PATH true)
    # the RPATH to be used when installing, but only if it's not a system directory
    list                    (FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
                                                               "${CMAKE_INSTALL_PREFIX}/lib"
                                                               isSystemDir)
    if                      ("${isSystemDir}" STREQUAL         "-1")
      set                   (CMAKE_INSTALL_RPATH               "${CMAKE_INSTALL_PREFIX}/lib")
    endif                   ()
endif                       (${CMAKE_SYSTEM_NAME} MATCHES      "Darwin")

# MPI
 if                         (NOT wrapped_mpi)
     find_package           (MPI REQUIRED)
     set                    (libraries ${libraries} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
 endif                      ()

# Threads
find_package                (Threads)
if                          (NOT diy_thread)
    message                 ("Diy threading is disabled; setting diy threads will have no effect")
    add_definitions         (-DDIY_NO_THREADS)
endif                       (NOT diy_thread)

# DIY
find_path                   (DIY_INCLUDE_DIRS diy/types.h)

if (WITH_VTK)
# VTK
find_package                (VTK REQUIRED)
include                     (${VTK_USE_FILE})
endif ()

# Include dirs
set                         (CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem")
include_directories         (${DIY_INCLUDE_DIRS}
                             ${VTK_INCLUDE_DIRS}
                             SYSTEM ${MPI_INCLUDE_PATH})

# Libraries
set                         (libraries
                             ${libraries}
                             ${CMAKE_THREAD_LIBS_INIT}
                             ${VTK_LIBRARIES})


# Find PnetCDF
find_path                   (PNETCDF_DIR lib/libpnetcdf.a)
set                         (PNETCDF_INCLUDE_DIR ${PNETCDF_DIR}/include)
set                         (PNETCDF_LIBRARY ${PNETCDF_DIR}/lib/libpnetcdf.a)
include_directories         (${PNETCDF_INCLUDE_DIR})

# MPE
if (WITH_MPE)
    set                     (MPE_LIB ${MPE_DIR}/lib/libmpe.a)
    set                     (LMPE_LIB ${MPE_DIR}/lib/liblmpe.a)
    set                     (libraries
                            ${MPE_LIB}
                            ${LMPE_LIB}
                            ${libraries})
    message                 ("MPE logging is enabled.")
    add_definitions         (-DMPE)
    include_directories     (${MPE_DIR}/include)
endif()

# find_package (NetCDF REQUIRED)
# include_directories (${NETCDF_INCLUDE_DIR})

# Subdirectories
add_subdirectory           (examples)


