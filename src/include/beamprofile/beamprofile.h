/**
 * @file beamprofile.h
 * @brief C API for beam profile analysis
 * 
 * This header provides a C-compatible interface for analyzing laser beam profiles.
 * The API supports various image formats, dark frame subtraction, gain map correction,
 * and calculation of beam metrics according to ISO 11146 and ISO 13694 standards.
 * 
 * @author Beam Profiler SDK Team
 * @version 1.0.0
 * @date 2024
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Opaque handle for beam profiler instance
 */
typedef void* BeamProfilerHandle;

/**
 * @brief Supported image formats for beam profile analysis
 */
typedef enum {
    BEAMPROFILE_MONO8,    ///< 8-bit monochrome image (0-255)
    BEAMPROFILE_MONO12,   ///< 12-bit monochrome image (0-4095)
    BEAMPROFILE_MONO16,   ///< 16-bit monochrome image (0-65535)
    BEAMPROFILE_RGB8      ///< 24-bit RGB image (8 bits per channel)
} BeamProfileImageFormat;

/**
 * @brief Configuration structure for beam profile analysis
 */
typedef struct {
    int width;                    ///< Image width in pixels
    int height;                   ///< Image height in pixels
    BeamProfileImageFormat format;  ///< Image format
    const void* dark_frame;       ///< Optional dark frame for subtraction (same format as input image)
    const void* gain_map;         ///< Optional gain map for correction (float32, same dimensions as image)
} BeamProfileConfig;

/**
 * @brief Results structure for beam profile analysis
 */
typedef struct {
    // Second moments (ISO 11146)
    double centroid_x;      ///< X centroid position in pixels
    double centroid_y;      ///< Y centroid position in pixels
    double sigma_xx;        ///< Second moment in x direction (pixels²)
    double sigma_yy;        ///< Second moment in y direction (pixels²)
    double sigma_xy;        ///< Cross moment (pixels²)
    double d4sigma_x;       ///< D4σ width in x direction (pixels)
    double d4sigma_y;       ///< D4σ width in y direction (pixels)
    double rotation;        ///< Beam rotation angle in radians
    double ellipticity;     ///< Beam ellipticity (ratio of major to minor axis)

    // Higher moments
    double skewness_x;      ///< Skewness in x direction
    double skewness_y;      ///< Skewness in y direction
    double kurtosis_x;      ///< Kurtosis in x direction
    double kurtosis_y;      ///< Kurtosis in y direction

    // Gaussian fit parameters
    double gaussian_amplitude;  ///< Amplitude of Gaussian fit
    double gaussian_center_x;   ///< X center of Gaussian fit in pixels
    double gaussian_center_y;   ///< Y center of Gaussian fit in pixels
    double gaussian_sigma_x;    ///< X sigma of Gaussian fit in pixels
    double gaussian_sigma_y;    ///< Y sigma of Gaussian fit in pixels
    double gaussian_theta;      ///< Rotation angle of Gaussian fit in radians
    double gaussian_fit_error;  ///< RMS error of Gaussian fit

    // ISO 13694 metrics
    double edge_steepness;      ///< Edge steepness
    double flatness;           ///< Flatness
    double rms_uniformity;     ///< RMS uniformity
    double d86_radius;         ///< D86 radius in pixels

    // Power in bucket metrics
    int power_in_bucket_size;  ///< Size of power_in_bucket array
    double* power_in_bucket;   ///< Power in bucket at user-specified radii
} BeamProfileMetrics;

/**
 * @brief Results structure for pointing jitter analysis
 */
typedef struct {
    double rms_jitter;         ///< RMS pointing jitter in pixels
    int psd_size;             ///< Size of PSD array
    double* psd;              ///< Power spectral density of pointing jitter
} BeamProfileJitter;

/**
 * @brief Create a new beam profile analysis configuration
 */
BeamProfileConfig* beamprofile_create_config(int width, int height, BeamProfileImageFormat format);

/**
 * @brief Destroy a beam profile analysis configuration
 */
void beamprofile_destroy_config(BeamProfileConfig* config);

/**
 * @brief Set dark frame for subtraction
 */
void beamprofile_config_set_dark_frame(BeamProfileConfig* config, const void* dark_frame);

/**
 * @brief Set gain map for correction
 */
void beamprofile_config_set_gain_map(BeamProfileConfig* config, const void* gain_map);

/**
 * @brief Analyze beam profile
 */
BeamProfileMetrics* beamprofile_analyze(const void* image, const BeamProfileConfig* config);

/**
 * @brief Destroy beam profile analysis results
 */
void beamprofile_destroy_metrics(BeamProfileMetrics* metrics);

/**
 * @brief Calculate M² from multiple measurements
 */
double beamprofile_calculate_m2(const double* z_positions, const double* beam_widths, int count);

/**
 * @brief Calculate pointing jitter metrics
 */
BeamProfileJitter* beamprofile_calculate_jitter(const double* x_positions, const double* y_positions, int count);

/**
 * @brief Destroy pointing jitter results
 */
void beamprofile_destroy_jitter(BeamProfileJitter* jitter);

#ifdef __cplusplus
}
#endif 