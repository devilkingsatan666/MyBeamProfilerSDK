#include "beamprofile/beamprofile.h"
#include "beamprofile/beamprofile.hpp"
#include <memory>
#include <vector>

extern "C" {

BeamProfileConfig* beamprofile_create_config(int width, int height, BeamProfileImageFormat format) {
    try {
        auto config = new BeamProfileConfig;
        config->width = width;
        config->height = height;
        config->format = format;
        config->dark_frame = nullptr;
        config->gain_map = nullptr;
        return config;
    } catch (...) {
        return nullptr;
    }
}

void beamprofile_destroy_config(BeamProfileConfig* config) {
    delete config;
}

void beamprofile_config_set_dark_frame(BeamProfileConfig* config, const void* dark_frame) {
    if (config) {
        config->dark_frame = dark_frame;
    }
}

void beamprofile_config_set_gain_map(BeamProfileConfig* config, const void* gain_map) {
    if (config) {
        config->gain_map = gain_map;
    }
}

BeamProfileMetrics* beamprofile_analyze(const void* image, const BeamProfileConfig* config) {
    if (!image || !config) {
        return nullptr;
    }

    try {
        // Convert C config to C++ config
        beamprofile::AnalysisConfig cpp_config{
            config->width,
            config->height,
            static_cast<beamprofile::ImageFormat>(config->format),
            config->dark_frame ? std::optional<const void*>(config->dark_frame) : std::nullopt,
            config->gain_map ? std::optional<const void*>(config->gain_map) : std::nullopt
        };

        // Analyze image using C++ API
        auto metrics = beamprofile::beam_analysis::analyze(image, cpp_config);

        // Convert C++ metrics to C metrics
        auto result = new BeamProfileMetrics;
        result->centroid_x = metrics.centroid_x;
        result->centroid_y = metrics.centroid_y;
        result->sigma_xx = metrics.sigma_xx;
        result->sigma_yy = metrics.sigma_yy;
        result->sigma_xy = metrics.sigma_xy;
        result->d4sigma_x = metrics.d4sigma_x;
        result->d4sigma_y = metrics.d4sigma_y;
        result->rotation = metrics.rotation;
        result->ellipticity = metrics.ellipticity;
        result->skewness_x = metrics.skewness_x;
        result->skewness_y = metrics.skewness_y;
        result->kurtosis_x = metrics.kurtosis_x;
        result->kurtosis_y = metrics.kurtosis_y;
        result->gaussian_amplitude = metrics.gaussian_amplitude;
        result->gaussian_center_x = metrics.gaussian_center_x;
        result->gaussian_center_y = metrics.gaussian_center_y;
        result->gaussian_sigma_x = metrics.gaussian_sigma_x;
        result->gaussian_sigma_y = metrics.gaussian_sigma_y;
        result->gaussian_theta = metrics.gaussian_theta;
        result->gaussian_fit_error = metrics.gaussian_fit_error;
        result->edge_steepness = metrics.edge_steepness;
        result->flatness = metrics.flatness;
        result->rms_uniformity = metrics.rms_uniformity;
        result->d86_radius = metrics.d86_radius;

        // Copy power in bucket data
        result->power_in_bucket_size = metrics.power_in_bucket.size();
        if (!metrics.power_in_bucket.empty()) {
            result->power_in_bucket = new double[metrics.power_in_bucket.size()];
            std::copy(metrics.power_in_bucket.begin(), metrics.power_in_bucket.end(),
                     result->power_in_bucket);
        } else {
            result->power_in_bucket = nullptr;
        }

        return result;
    } catch (...) {
        return nullptr;
    }
}

void beamprofile_destroy_metrics(BeamProfileMetrics* metrics) {
    if (metrics) {
        delete[] metrics->power_in_bucket;
        delete metrics;
    }
}

double beamprofile_calculate_m2(const double* z_positions, const double* beam_widths, int count) {
    if (!z_positions || !beam_widths || count < 3) {
        return -1.0;
    }

    try {
        std::vector<double> z_vec(z_positions, z_positions + count);
        std::vector<double> w_vec(beam_widths, beam_widths + count);
        return beamprofile::beam_analysis::calculateM2(z_vec, w_vec);
    } catch (...) {
        return -1.0;
    }
}

BeamProfileJitter* beamprofile_calculate_jitter(const double* x_positions, const double* y_positions, int count) {
    if (!x_positions || !y_positions || count < 1) {
        return nullptr;
    }

    try {
        std::vector<std::pair<double, double>> centroids;
        centroids.reserve(count);
        for (int i = 0; i < count; ++i) {
            centroids.emplace_back(x_positions[i], y_positions[i]);
        }

        auto [rms_jitter, psd] = beamprofile::beam_analysis::calculateJitter(centroids);

        auto result = new BeamProfileJitter;
        result->rms_jitter = rms_jitter;
        result->psd_size = psd.size();
        result->psd = new double[psd.size()];
        std::copy(psd.begin(), psd.end(), result->psd);

        return result;
    } catch (...) {
        return nullptr;
    }
}

void beamprofile_destroy_jitter(BeamProfileJitter* jitter) {
    if (jitter) {
        delete[] jitter->psd;
        delete jitter;
    }
}

} // extern "C" 