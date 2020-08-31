include_directories(/home/mraj/libraries/Zoltan_v3.83/install/include)

add_executable              (ptrace-iexchange ptrace.cpp advect.cpp misc.cpp )
target_compile_definitions  (ptrace-iexchange PUBLIC IEXCHANGE=1)
add_executable              (ptrace-exchange ptrace.cpp advect.cpp misc.cpp)
target_compile_definitions  (ptrace-exchange PUBLIC IEXCHANGE=0)


add_executable(btrace btrace.cpp misc.cpp btrace/partitioning.cpp btrace/io.cpp btrace/advection.cpp btrace/bmisc.cpp btrace/bblock.cpp)


target_link_libraries       (ptrace-iexchange ${libraries} ${PNETCDF_LIBRARY} )
target_link_libraries       (ptrace-exchange  ${libraries} ${PNETCDF_LIBRARY})
target_link_libraries       (btrace  ${libraries} ${PNETCDF_LIBRARY} /home/mraj/libraries/Zoltan_v3.83/install/lib/libzoltan.a)



if (WITH_VTK)
target_compile_definitions(ptrace-iexchange PRIVATE WITH_VTK)
target_compile_definitions(ptrace-exchange  PRIVATE WITH_VTK)
endif ()

if (WITH_TIMEINFO)
target_compile_definitions(ptrace-iexchange PRIVATE WITH_TIMEINFO)
target_compile_definitions(ptrace-exchange  PRIVATE WITH_TIMEINFO)
endif()



install(TARGETS ptrace-iexchange
        DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/particle-tracing
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_WRITE WORLD_EXECUTE)

install(TARGETS ptrace-exchange
        DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/particle-tracing
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_WRITE WORLD_EXECUTE)

        install(TARGETS btrace
        DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/particle-tracing
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_WRITE WORLD_EXECUTE)

install(FILES PLUME_TEST TORNADO_TEST NEK_TEST1 plot_counters.py compare_segments.py bebop0.sh bebop1.sh scripts/theta1.sh scripts/bebop_isabel.sh scripts/baseline.sh
        DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/particle-tracing
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_WRITE WORLD_EXECUTE)
