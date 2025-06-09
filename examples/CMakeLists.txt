# Find libtiff
find_package(TIFF REQUIRED)

# Add analyze_image example
add_executable(analyze_image analyze_image.cpp)
target_link_libraries(analyze_image
    PRIVATE
        beamprofile_c
        TIFF::TIFF
) 