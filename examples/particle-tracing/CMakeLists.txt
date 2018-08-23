add_executable              (ptrace ptrace.cpp advect.cpp)
target_link_libraries       (ptrace ${libraries})

install(TARGETS ptrace
        DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/particle-tracing
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_WRITE WORLD_EXECUTE)

install(FILES PLUME_TEST TORNADO_TEST
        DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/particle-tracing
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_WRITE WORLD_EXECUTE)
