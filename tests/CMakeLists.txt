# Find required packages
find_package(GTest REQUIRED)
find_package(TIFF REQUIRED)
find_package(nlohmann_json REQUIRED)

# Add test executable
add_executable(beamprofile_test beamprofile_test.cpp)
target_link_libraries(beamprofile_test
    PRIVATE
        beamprofile_c
        GTest::GTest
        GTest::Main
        TIFF::TIFF
        nlohmann_json::nlohmann_json
)

# Add test
add_test(NAME beamprofile_test COMMAND beamprofile_test) 