add_executable              (rgrid rgrid.cpp)
add_executable              (diy-rgrid diy-rgrid.cpp)
target_link_libraries       (rgrid ${libraries})
target_link_libraries       (diy-rgrid ${libraries})

install(TARGETS rgrid diy-rgrid
        DESTINATION ${CMAKE_INSTALL_PREFIX}/examples/vtk-data-models
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_WRITE GROUP_EXECUTE
        WORLD_READ WORLD_WRITE WORLD_EXECUTE)
