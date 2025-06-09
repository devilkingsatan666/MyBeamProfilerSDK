#pragma once

#include <cstdint>
#include <memory>
#include <vector>
#include <array>
#include <optional>

namespace beamprofile {

/**
 * @brief Supported image formats
 */
enum class ImageFormat {
    MONO8,    ///< 8-bit monochrome
    MONO12,   ///< 12-bit monochrome
    MONO16,   ///< 16-bit monochrome
    RGB8      ///< 24-bit RGB
};

/**
 * @brief Beam profile analysis results
 */
struct BeamMetrics {
    // Second moments (ISO 11146)
    double centroid_x;      ///< X centroid position
    double centroid_y;      ///< Y centroid position
    double sigma_xx;        ///< Second moment in x
    double sigma_yy;        ///< Second moment in y
    double sigma_xy;        ///< Cross moment
    double d4sigma_x;       ///< D4σ width in x
    double d4sigma_y;       ///< D4σ width in y
    double rotation;        ///< Beam rotation angle (radians)
    double ellipticity;     ///< Beam ellipticity

    // Higher moments
    double skewness_x;      ///< Skewness in x
    double skewness_y;      ///< Skewness in y
    double kurtosis_x;      ///< Kurtosis in x
    double kurtosis_y;      ///< Kurtosis in y

    // Gaussian fit parameters
    double gaussian_amplitude;  ///< Amplitude of Gaussian fit
    double gaussian_center_x;   ///< X center of Gaussian fit
    double gaussian_center_y;   ///< Y center of Gaussian fit
    double gaussian_sigma_x;    ///< X sigma of Gaussian fit
    double gaussian_sigma_y;    ///< Y sigma of Gaussian fit
    double gaussian_theta;      ///< Rotation angle of Gaussian fit (radians)
    double gaussian_fit_error;  ///< RMS error of Gaussian fit

    // ISO 13694 metrics
    double edge_steepness;      ///< Edge steepness
    double flatness;           ///< Flatness
    double rms_uniformity;     ///< RMS uniformity

    // Power in bucket metrics
    double d86_radius;         ///< D86 radius
    std::vector<double> power_in_bucket;  ///< Power in bucket at user-specified radii
};

/**
 * @brief Configuration for beam profile analysis
 */
struct AnalysisConfig {
    const int width;                    ///< Image width in pixels
    const int height;                   ///< Image height in pixels
    const ImageFormat format;           ///< Image format
    const std::optional<const void*> dark_frame;  ///< Optional dark frame for subtraction
    const std::optional<const void*> gain_map;    ///< Optional gain map for correction
};

/**
 * @brief Namespace containing pure functions for beam profile analysis
 */
namespace beam_analysis {

/**
 * @brief Convert image to float32 format
 * @param image Input image data
 * @param config Analysis configuration
 * @return Vector of float32 values
 */
std::vector<float> convertToFloat(const void* image, const AnalysisConfig& config);

/**
 * @brief Apply dark frame subtraction
 * @param image Input image data
 * @param dark_frame Dark frame data
 * @return Corrected image data
 */
std::vector<float> applyDarkFrame(const std::vector<float>& image, const void* dark_frame);

/**
 * @brief Apply gain map correction
 * @param image Input image data
 * @param gain_map Gain map data
 * @return Corrected image data
 */
std::vector<float> applyGainMap(const std::vector<float>& image, const void* gain_map);

/**
 * @brief Calculate second moments and higher moments
 * @param image Input image data
 * @param config Analysis configuration
 * @return Beam metrics with moment calculations
 */
BeamMetrics calculateMoments(const std::vector<float>& image, const AnalysisConfig& config);

/**
 * @brief Fit Gaussian to beam profile
 * @param image Input image data
 * @param config Analysis configuration
 * @param moments Previously calculated moments
 * @return Updated beam metrics with Gaussian fit
 */
BeamMetrics fitGaussian(const std::vector<float>& image, const AnalysisConfig& config, const BeamMetrics& moments);

/**
 * @brief Calculate ISO 13694 metrics
 * @param image Input image data
 * @param config Analysis configuration
 * @param moments Previously calculated moments
 * @return Updated beam metrics with ISO 13694 metrics
 */
BeamMetrics calculateISO13694Metrics(const std::vector<float>& image, const AnalysisConfig& config, const BeamMetrics& moments);

/**
 * @brief Analyze beam profile
 * @param image Input image data
 * @param config Analysis configuration
 * @return Complete beam metrics
 */
BeamMetrics analyze(const void* image, const AnalysisConfig& config);

/**
 * @brief Calculate M² from multiple measurements
 * @param z_positions Vector of z positions
 * @param beam_widths Vector of beam widths (D4σ)
 * @return M² value
 */
double calculateM2(const std::vector<double>& z_positions, const std::vector<double>& beam_widths);

/**
 * @brief Calculate pointing jitter metrics
 * @param centroids Vector of centroid positions over time
 * @return RMS jitter and PSD
 */
std::pair<double, std::vector<double>> calculateJitter(
    const std::vector<std::pair<double, double>>& centroids);

} // namespace beam_analysis

} // namespace beamprofile 