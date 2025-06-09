# Create library targets
if(BUILD_CPP_API)
    add_library(beamprofile_cpp SHARED
        beamprofile_impl.cpp
    )
    target_include_directories(beamprofile_cpp
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(beamprofile_cpp PUBLIC cxx_std_17)
    
    # Set output name and properties
    set_target_properties(beamprofile_cpp PROPERTIES
        OUTPUT_NAME beamprofile_cpp
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
endif()

if(BUILD_C_API)
    add_library(beamprofile_c SHARED
        beamprofile_c.cpp
    )
    target_include_directories(beamprofile_c
        PUBLIC
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}
    )
    target_compile_features(beamprofile_c PUBLIC cxx_std_17)
    if(BUILD_CPP_API)
        target_link_libraries(beamprofile_c PRIVATE beamprofile_cpp)
    endif()
    
    # Set output name and properties
    set_target_properties(beamprofile_c PROPERTIES
        OUTPUT_NAME beamprofile_c
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        WINDOWS_EXPORT_ALL_SYMBOLS ON
    )
endif()

# Install targets
if(BUILD_CPP_API)
    install(TARGETS beamprofile_cpp
        EXPORT beamprofile-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
    )
endif()

if(BUILD_C_API)
    install(TARGETS beamprofile_c
        EXPORT beamprofile-targets
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
        INCLUDES DESTINATION include
    )
endif()

# Install headers
install(DIRECTORY include/beamprofile
    DESTINATION include
    FILES_MATCHING PATTERN "*.h*"
)

# Export targets
install(EXPORT beamprofile-targets
    FILE beamprofile-targets.cmake
    NAMESPACE beamprofile::
    DESTINATION lib/cmake/beamprofile
) 