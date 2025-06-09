#include "beamprofile/beamprofile.h"
#include <tiffio.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Helper function to print metrics
void print_metrics(const BeamProfileMetrics* metrics) {
    if (!metrics) {
        std::cerr << "Error: Invalid metrics" << std::endl;
        return;
    }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nBeam Profile Analysis Results:\n";
    std::cout << "=============================\n\n";

    // Second moments (ISO 11146)
    std::cout << "Second Moments (ISO 11146):\n";
    std::cout << "  Centroid: (" << metrics->centroid_x << ", " << metrics->centroid_y << ")\n";
    std::cout << "  D4σ Width: " << metrics->d4sigma_x << " x " << metrics->d4sigma_y << "\n";
    std::cout << "  Rotation: " << metrics->rotation * 180.0 / M_PI << "°\n";
    std::cout << "  Ellipticity: " << metrics->ellipticity << "\n\n";

    // Higher moments
    std::cout << "Higher Moments:\n";
    std::cout << "  Skewness: (" << metrics->skewness_x << ", " << metrics->skewness_y << ")\n";
    std::cout << "  Kurtosis: (" << metrics->kurtosis_x << ", " << metrics->kurtosis_y << ")\n\n";

    // Gaussian fit
    std::cout << "Gaussian Fit:\n";
    std::cout << "  Amplitude: " << metrics->gaussian_amplitude << "\n";
    std::cout << "  Center: (" << metrics->gaussian_center_x << ", " << metrics->gaussian_center_y << ")\n";
    std::cout << "  Sigma: (" << metrics->gaussian_sigma_x << ", " << metrics->gaussian_sigma_y << ")\n";
    std::cout << "  Rotation: " << metrics->gaussian_theta * 180.0 / M_PI << "°\n";
    std::cout << "  Fit Error: " << metrics->gaussian_fit_error << "\n\n";

    // ISO 13694 metrics
    std::cout << "ISO 13694 Metrics:\n";
    std::cout << "  Edge Steepness: " << metrics->edge_steepness << "\n";
    std::cout << "  Flatness: " << metrics->flatness << "\n";
    std::cout << "  RMS Uniformity: " << metrics->rms_uniformity << "\n";
    std::cout << "  D86 Radius: " << metrics->d86_radius << "\n\n";

    // Power in bucket
    if (metrics->power_in_bucket && metrics->power_in_bucket_size > 0) {
        std::cout << "Power in Bucket:\n";
        for (int i = 0; i < metrics->power_in_bucket_size; ++i) {
            std::cout << "  Radius " << i << ": " << metrics->power_in_bucket[i] << "\n";
        }
    }
}

// Helper function to print jitter results
void print_jitter(const BeamProfileJitter* jitter) {
    if (!jitter) {
        std::cerr << "Error: Invalid jitter results" << std::endl;
        return;
    }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nPointing Jitter Analysis:\n";
    std::cout << "=======================\n\n";
    std::cout << "RMS Jitter: " << jitter->rms_jitter << "\n\n";

    if (jitter->psd && jitter->psd_size > 0) {
        std::cout << "Power Spectral Density:\n";
        for (int i = 0; i < jitter->psd_size; ++i) {
            std::cout << "  Frequency " << i << ": " << jitter->psd[i] << "\n";
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <image.tif> [dark_frame.tif] [gain_map.tif]" << std::endl;
        return 1;
    }

    // Open and read TIFF image
    TIFF* tif = TIFFOpen(argv[1], "r");
    if (!tif) {
        std::cerr << "Error: Could not open image file" << std::endl;
        return 1;
    }

    uint32_t width, height;
    uint16_t bits_per_sample, samples_per_pixel;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);

    // Determine image format
    BeamProfileImageFormat format;
    if (samples_per_pixel == 1) {
        switch (bits_per_sample) {
            case 8: format = BEAMPROFILE_MONO8; break;
            case 12: format = BEAMPROFILE_MONO12; break;
            case 16: format = BEAMPROFILE_MONO16; break;
            default:
                std::cerr << "Error: Unsupported bit depth" << std::endl;
                TIFFClose(tif);
                return 1;
        }
    } else if (samples_per_pixel == 3 && bits_per_sample == 8) {
        format = BEAMPROFILE_RGB8;
    } else {
        std::cerr << "Error: Unsupported image format" << std::endl;
        TIFFClose(tif);
        return 1;
    }

    // Allocate image buffer
    size_t image_size = width * height * samples_per_pixel * (bits_per_sample / 8);
    std::vector<uint8_t> image_data(image_size);
    TIFFReadRawStrip(tif, 0, image_data.data(), image_size);
    TIFFClose(tif);

    // Create configuration
    auto config = beamprofile_create_config(width, height, format);
    if (!config) {
        std::cerr << "Error: Could not create configuration" << std::endl;
        return 1;
    }

    // Load dark frame if provided
    if (argc > 2) {
        TIFF* dark_tif = TIFFOpen(argv[2], "r");
        if (dark_tif) {
            std::vector<uint8_t> dark_frame(image_size);
            TIFFReadRawStrip(dark_tif, 0, dark_frame.data(), image_size);
            TIFFClose(dark_tif);
            beamprofile_config_set_dark_frame(config, dark_frame.data());
        }
    }

    // Load gain map if provided
    if (argc > 3) {
        TIFF* gain_tif = TIFFOpen(argv[3], "r");
        if (gain_tif) {
            std::vector<float> gain_map(width * height);
            TIFFReadRawStrip(gain_tif, 0, gain_map.data(), width * height * sizeof(float));
            TIFFClose(gain_tif);
            beamprofile_config_set_gain_map(config, gain_map.data());
        }
    }

    // Analyze beam profile
    auto metrics = beamprofile_analyze(image_data.data(), config);
    if (!metrics) {
        std::cerr << "Error: Analysis failed" << std::endl;
        beamprofile_destroy_config(config);
        return 1;
    }

    // Print results
    print_metrics(metrics);

    // Example of M² calculation
    std::vector<double> z_positions = {0.0, 50.0, 100.0, 150.0, 200.0};
    std::vector<double> beam_widths = {
        metrics->d4sigma_x,
        metrics->d4sigma_x * 1.2,
        metrics->d4sigma_x * 1.5,
        metrics->d4sigma_x * 1.8,
        metrics->d4sigma_x * 2.0
    };

    double m2 = beamprofile_calculate_m2(z_positions.data(), beam_widths.data(), z_positions.size());
    if (m2 >= 0.0) {
        std::cout << "\nM² Calculation:\n";
        std::cout << "==============\n\n";
        std::cout << "M² = " << m2 << "\n";
    }

    // Example of jitter calculation
    std::vector<double> x_positions = {metrics->centroid_x, metrics->centroid_x + 0.1, metrics->centroid_x - 0.1};
    std::vector<double> y_positions = {metrics->centroid_y, metrics->centroid_y + 0.1, metrics->centroid_y - 0.1};

    auto jitter = beamprofile_calculate_jitter(x_positions.data(), y_positions.data(), x_positions.size());
    if (jitter) {
        print_jitter(jitter);
        beamprofile_destroy_jitter(jitter);
    }

    // Cleanup
    beamprofile_destroy_metrics(metrics);
    beamprofile_destroy_config(config);

    return 0;
} 