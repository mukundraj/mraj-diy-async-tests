add_executable              (ptrace-iexchange ptrace.cpp advect.cpp)
target_compile_definitions  (ptrace-iexchange PUBLIC IEXCHANGE=1)
add_executable              (ptrace-exchange ptrace.cpp advect.cpp)
target_compile_definitions  (ptrace-exchange PUBLIC IEXCHANGE=0)


target_link_libraries       (ptrace-iexchange ${libraries} ${PNETCDF_LIBRARY})
target_link_libraries       (ptrace-exchange  ${libraries} ${PNETCDF_LIBRARY})



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

install(FILES PLUME_TEST TORNADO_TEST NEK_TEST1 plot_counters.py compare_segments.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/particle-tracing
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_WRITE WORLD_EXECUTE)
