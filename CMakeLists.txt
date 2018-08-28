project                     (diy2-vtk7)
cmake_minimum_required      (VERSION 2.8)

option                      (diy_thread     "Enable diy threading"             OFF)
option (WITH_VTK OFF)
option (WITH_TIMEINFO ON)

if                          (bgq)
    add_definitions         (-DBGQ)
endif                       (bgq)

# C++11
set                         (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# OSX flags
if                          (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    add_definitions	    (-DMAC_OSX)
    set			    (CMAKE_MACOSX_RPATH	    on)
endif                       (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

 # MPI
find_package                (MPI REQUIRED)
if                          (NOT bgq)
    set                     (libraries ${libraries} ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
endif                       (NOT bgq)

# Threads
find_package		    (Threads)
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
find_path(PNETCDF_DIR lib/libpnetcdf.a)
set (PNETCDF_INCLUDE_DIR ${PNETCDF_DIR}/include)
set (PNETCDF_LIBRARY ${PNETCDF_DIR}/lib/libpnetcdf.a)
include_directories (${PNETCDF_INCLUDE_DIR})


# find_package (NetCDF REQUIRED)
# include_directories (${NETCDF_INCLUDE_DIR})

# Subdirectories
add_subdirectory           (examples)


