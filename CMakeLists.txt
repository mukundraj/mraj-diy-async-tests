project                     (diy2-vtk7)
cmake_minimum_required      (VERSION 2.8)

option                      (diy_thread     "Enable diy threading"             OFF)

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

## VTK
#find_package                (VTK REQUIRED)
#include                     (${VTK_USE_FILE})

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

# Subdirectories
add_subdirectory           (examples)

