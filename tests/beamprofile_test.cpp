#include <beamprofile/beamprofile.h>
#include <gtest/gtest.h>
#include <tiffio.h>
#include <fstream>
#include <nlohmann/json.hpp>
#include <filesystem>

using json = nlohmann::json;

class BeamProfileTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Load baseline metrics
        std::ifstream baseline_file("baseline_metrics.json");
        if (!baseline_file.is_open()) {
            FAIL() << "Could not open baseline_metrics.json";
        }
        baseline_metrics = json::parse(baseline_file);
    }

    json baseline_metrics;
};

TEST_F(BeamProfileTest, CompareWithBaseline) {
    // Get all TIFF files in the data directory
    std::vector<std::string> tiff_files;
    std::string data_dir_mosa = "../data/mosa";
    std::string data_dir_oap = "../data/oap";

    for (const auto& entry : std::filesystem::directory_iterator(data_dir_oap)) {
        if (entry.path().extension() == ".tif") {
            tiff_files.push_back(entry.path().string());
        }
    }

    for (const auto& entry : std::filesystem::directory_iterator(data_dir_mosa)) {
        if (entry.path().extension() == ".tif") {
            tiff_files.push_back(entry.path().string());
        }
    }

    for (const auto& tiff_file : tiff_files) {
        // Open TIFF file
        TIFF* tif = TIFFOpen(tiff_file.c_str(), "r");
        ASSERT_NE(tif, nullptr) << "Could not open TIFF file: " << tiff_file;

        // Get image dimensions
        uint32_t width, height;
        TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
        TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);

        // Get bits per sample
        uint16_t bits_per_sample;
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);

        // Get samples per pixel
        uint16_t samples_per_pixel;
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);

        // Determine image format
        BeamProfileImageFormat format;
        if (samples_per_pixel == 1) {
            switch (bits_per_sample) {
                case 8: format = BEAMPROFILE_MONO8; break;
                case 12: format = BEAMPROFILE_MONO12; break;
                case 16: format = BEAMPROFILE_MONO16; break;
                default:
                    FAIL() << "Unsupported bits per sample: " << bits_per_sample;
            }
        } else if (samples_per_pixel == 3 && bits_per_sample == 8) {
            format = BEAMPROFILE_RGB8;
        } else {
            FAIL() << "Unsupported image format";
        }

        // Allocate buffer for image data
        size_t bytes_per_sample = (bits_per_sample + 7) / 8;
        size_t bytes_per_pixel = bytes_per_sample * samples_per_pixel;
        size_t buffer_size = width * height * bytes_per_pixel;
        std::vector<uint8_t> buffer(buffer_size);

        // Read image data
        ASSERT_NE(TIFFReadEncodedStrip(tif, 0, buffer.data(), buffer_size), -1)
            << "Could not read image data";

        // Create beam profiler configuration
        BeamProfileConfig* config = beamprofile_create_config(width, height, format);
        ASSERT_NE(config, nullptr) << "Could not create beam profiler configuration";

        // Analyze beam profile
        BeamProfileMetrics* metrics = beamprofile_analyze(buffer.data(), config);
        ASSERT_NE(metrics, nullptr) << "Could not analyze beam profile";

        // Get baseline metrics for this file
        std::string filename = std::filesystem::path(tiff_file).filename().string();
        ASSERT_TRUE(baseline_metrics.contains(filename))
            << "No baseline metrics found for " << filename;
        const auto& baseline = baseline_metrics[filename];

        // Compare metrics with baseline
        const double tolerance = 1e-3;  // 0.1% tolerance

        EXPECT_NEAR(metrics->centroid_x, baseline["centroid_x"], tolerance)
            << "Centroid X mismatch for " << filename;
        EXPECT_NEAR(metrics->centroid_y, baseline["centroid_y"], tolerance)
            << "Centroid Y mismatch for " << filename;
        EXPECT_NEAR(metrics->sigma_xx, baseline["sigma_xx"], tolerance)
            << "Sigma XX mismatch for " << filename;
        EXPECT_NEAR(metrics->sigma_yy, baseline["sigma_yy"], tolerance)
            << "Sigma YY mismatch for " << filename;
        EXPECT_NEAR(metrics->sigma_xy, baseline["sigma_xy"], tolerance)
            << "Sigma XY mismatch for " << filename;
        EXPECT_NEAR(metrics->d4sigma_x, baseline["d4sigma_x"], tolerance)
            << "D4σ X mismatch for " << filename;
        EXPECT_NEAR(metrics->d4sigma_y, baseline["d4sigma_y"], tolerance)
            << "D4σ Y mismatch for " << filename;
        EXPECT_NEAR(metrics->rotation, baseline["rotation"], tolerance)
            << "Rotation mismatch for " << filename;
        EXPECT_NEAR(metrics->ellipticity, baseline["ellipticity"], tolerance)
            << "Ellipticity mismatch for " << filename;
        EXPECT_NEAR(metrics->skewness_x, baseline["skewness_x"], tolerance)
            << "Skewness X mismatch for " << filename;
        EXPECT_NEAR(metrics->skewness_y, baseline["skewness_y"], tolerance)
            << "Skewness Y mismatch for " << filename;
        EXPECT_NEAR(metrics->kurtosis_x, baseline["kurtosis_x"], tolerance)
            << "Kurtosis X mismatch for " << filename;
        EXPECT_NEAR(metrics->kurtosis_y, baseline["kurtosis_y"], tolerance)
            << "Kurtosis Y mismatch for " << filename;
        EXPECT_NEAR(metrics->edge_steepness, baseline["edge_steepness"], tolerance)
            << "Edge steepness mismatch for " << filename;
        EXPECT_NEAR(metrics->flatness, baseline["flatness"], tolerance)
            << "Flatness mismatch for " << filename;
        EXPECT_NEAR(metrics->rms_uniformity, baseline["rms_uniformity"], tolerance)
            << "RMS uniformity mismatch for " << filename;
        EXPECT_NEAR(metrics->d86_radius, baseline["d86_radius"], tolerance)
            << "D86 radius mismatch for " << filename;

        // Clean up
        beamprofile_destroy_metrics(metrics);
        beamprofile_destroy_config(config);
        TIFFClose(tif);
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 