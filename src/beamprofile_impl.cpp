#include "beamprofile/beamprofile.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <array>
#include <complex>
#include <immintrin.h>  // For AVX/SSE instructions
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace beamprofile {
namespace beam_analysis {

// Forward declaration of fft function
// Helper function for FFT
void fft(std::vector<std::complex<double>>& x) {
    const int n = static_cast<int>(x.size());
    if (n <= 1) return;
    
    // Bit reversal
    for (int i = 0; i < n; ++i) {
        int j = 0;
        for (int k = 0; k < static_cast<int>(std::log2(n)); ++k) {
            j = (j << 1) | (i >> k & 1);
        }
        if (j > i) {
            std::swap(x[i], x[j]);
        }
    }
    
    // Cooley-Tukey FFT
    for (int size = 2; size <= n; size *= 2) {
        double angle = -2 * M_PI / size;
        std::complex<double> w(1, 0);
        std::complex<double> wn(std::cos(angle), std::sin(angle));
        
        for (int i = 0; i < n; i += size) {
            std::complex<double> temp(1, 0);
            for (int j = 0; j < size / 2; ++j) {
                std::complex<double> u = x[i + j];
                std::complex<double> t = temp * x[i + j + size / 2];
                x[i + j] = u + t;
                x[i + j + size / 2] = u - t;
                temp *= wn;
            }
        }
    }
}

std::vector<float> convertToFloat(const void* image, const AnalysisConfig& config) {
    std::vector<float> output(config.width * config.height);
    
    #ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp for schedule(static)
        for (int i = 0; i < config.width * config.height; ++i) {
            switch (config.format) {
                case ImageFormat::MONO8:
                    output[i] = static_cast<float>(static_cast<const uint8_t*>(image)[i]);
                    break;
                case ImageFormat::MONO12:
                    output[i] = static_cast<float>(static_cast<const uint16_t*>(image)[i] & 0x0FFF);
                    break;
                case ImageFormat::MONO16:
                    output[i] = static_cast<float>(static_cast<const uint16_t*>(image)[i]);
                    break;
                case ImageFormat::RGB8:
                    const uint8_t* rgb = static_cast<const uint8_t*>(image) + i * 3;
                    output[i] = 0.299f * rgb[0] + 0.587f * rgb[1] + 0.114f * rgb[2];
                    break;
            }
        }
    }
    #else
    switch (config.format) {
        case ImageFormat::MONO8:
            for (int i = 0; i < config.width * config.height; ++i) {
                output[i] = static_cast<float>(static_cast<const uint8_t*>(image)[i]);
            }
            break;
        case ImageFormat::MONO12:
            for (int i = 0; i < config.width * config.height; ++i) {
                output[i] = static_cast<float>(static_cast<const uint16_t*>(image)[i] & 0x0FFF);
            }
            break;
        case ImageFormat::MONO16:
            for (int i = 0; i < config.width * config.height; ++i) {
                output[i] = static_cast<float>(static_cast<const uint16_t*>(image)[i]);
            }
            break;
        case ImageFormat::RGB8:
            for (int i = 0; i < config.width * config.height; ++i) {
                const uint8_t* rgb = static_cast<const uint8_t*>(image) + i * 3;
                output[i] = 0.299f * rgb[0] + 0.587f * rgb[1] + 0.114f * rgb[2];
            }
            break;
    }
    #endif
    
    return output;
}

std::vector<float> applyDarkFrame(const std::vector<float>& image, const void* dark_frame) {
    std::vector<float> result = image;  // Create a copy
    const uint8_t* df = static_cast<const uint8_t*>(dark_frame);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < image.size(); ++i) {
        result[i] -= static_cast<float>(df[i]);
    }
    #else
    for (int i = 0; i < image.size(); ++i) {
        result[i] -= static_cast<float>(df[i]);
    }
    #endif
    
    return result;
}

std::vector<float> applyGainMap(const std::vector<float>& image, const void* gain_map) {
    std::vector<float> result = image;  // Create a copy
    const float* gm = static_cast<const float*>(gain_map);
    
    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < image.size(); ++i) {
        result[i] *= gm[i];
    }
    #else
    for (int i = 0; i < image.size(); ++i) {
        result[i] *= gm[i];
    }
    #endif
    
    return result;
}

BeamMetrics calculateMoments(const std::vector<float>& image, const AnalysisConfig& config) {
    // Calculate total intensity using SIMD
    float total = 0.0f;
    #ifdef __AVX2__
    __m256 sum = _mm256_setzero_ps();
    int i = 0;
    for (; i + 8 <= config.width * config.height; i += 8) {
        __m256 data = _mm256_loadu_ps(&image[i]);
        sum = _mm256_add_ps(sum, data);
    }
    float sum_array[8];
    _mm256_storeu_ps(sum_array, sum);
    for (int j = 0; j < 8; ++j) {
        total += sum_array[j];
    }
    for (; i < config.width * config.height; ++i) {
        total += image[i];
    }
    #else
    for (int i = 0; i < config.width * config.height; ++i) {
        total += image[i];
    }
    #endif

    // Calculate centroid and second moments using OpenMP
    float sum_x = 0.0f, sum_y = 0.0f;
    float sum_xx = 0.0f, sum_yy = 0.0f, sum_xy = 0.0f;
    float sum_xxx = 0.0f, sum_yyy = 0.0f;
    float sum_xxxx = 0.0f, sum_yyyy = 0.0f;

    #ifdef _OPENMP
    #pragma omp parallel reduction(+:sum_x,sum_y,sum_xx,sum_yy,sum_xy,sum_xxx,sum_yyy,sum_xxxx,sum_yyyy)
    {
        #pragma omp for schedule(static)
        for (int y = 0; y < config.height; ++y) {
            for (int x = 0; x < config.width; ++x) {
                float intensity = image[y * config.width + x];
                sum_x += x * intensity;
                sum_y += y * intensity;
            }
        }
    }
    #else
    for (int y = 0; y < config.height; ++y) {
        for (int x = 0; x < config.width; ++x) {
            float intensity = image[y * config.width + x];
            sum_x += x * intensity;
            sum_y += y * intensity;
        }
    }
    #endif

    double centroid_x = sum_x / total;
    double centroid_y = sum_y / total;

    #ifdef _OPENMP
    #pragma omp parallel reduction(+:sum_xx,sum_yy,sum_xy,sum_xxx,sum_yyy,sum_xxxx,sum_yyyy)
    {
        #pragma omp for schedule(static)
        for (int y = 0; y < config.height; ++y) {
            for (int x = 0; x < config.width; ++x) {
                float intensity = image[y * config.width + x];
                float dx = x - centroid_x;
                float dy = y - centroid_y;
                
                sum_xx += dx * dx * intensity;
                sum_yy += dy * dy * intensity;
                sum_xy += dx * dy * intensity;
                sum_xxx += dx * dx * dx * intensity;
                sum_yyy += dy * dy * dy * intensity;
                sum_xxxx += dx * dx * dx * dx * intensity;
                sum_yyyy += dy * dy * dy * dy * intensity;
            }
        }
    }
    #else
    for (int y = 0; y < config.height; ++y) {
        for (int x = 0; x < config.width; ++x) {
            float intensity = image[y * config.width + x];
            float dx = x - centroid_x;
            float dy = y - centroid_y;
            
            sum_xx += dx * dx * intensity;
            sum_yy += dy * dy * intensity;
            sum_xy += dx * dy * intensity;
            sum_xxx += dx * dx * dx * intensity;
            sum_yyy += dy * dy * dy * intensity;
            sum_xxxx += dx * dx * dx * dx * intensity;
            sum_yyyy += dy * dy * dy * dy * intensity;
        }
    }
    #endif

    double sigma_xx = sum_xx / total;
    double sigma_yy = sum_yy / total;
    double sigma_xy = sum_xy / total;
    double d4sigma_x = 4.0 * std::sqrt(sigma_xx);
    double d4sigma_y = 4.0 * std::sqrt(sigma_yy);
    double rotation = 0.5 * std::atan2(2.0 * sigma_xy, sigma_xx - sigma_yy);
    double ellipticity = std::sqrt(std::pow(sigma_xx - sigma_yy, 2) + 4.0 * std::pow(sigma_xy, 2)) /
                        (sigma_xx + sigma_yy);
    
    double skewness_x = (sum_xxx / total) / std::pow(sigma_xx, 1.5);
    double skewness_y = (sum_yyy / total) / std::pow(sigma_yy, 1.5);
    double kurtosis_x = (sum_xxxx / total) / std::pow(sigma_xx, 2) - 3.0;
    double kurtosis_y = (sum_yyyy / total) / std::pow(sigma_yy, 2) - 3.0;

    return BeamMetrics{
        centroid_x, centroid_y,
        sigma_xx, sigma_yy, sigma_xy,
        d4sigma_x, d4sigma_y,
        rotation, ellipticity,
        skewness_x, skewness_y,
        kurtosis_x, kurtosis_y,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  // Gaussian fit parameters (to be filled later)
        0.0, 0.0, 0.0,  // ISO 13694 metrics (to be filled later)
        0.0, {}  // Power in bucket metrics (to be filled later)
    };
}

BeamMetrics fitGaussian(const std::vector<float>& image, const AnalysisConfig& config, const BeamMetrics& moments) {
    // Initial guess from second moments
    double p0[6] = {
        *std::max_element(image.begin(), image.end()),  // amplitude
        moments.centroid_x,                             // center x
        moments.centroid_y,                             // center y
        moments.d4sigma_x / 4.0,                        // sigma x
        moments.d4sigma_y / 4.0,                        // sigma y
        moments.rotation                                // rotation
    };

    // Levenberg-Marquardt parameters
    const int max_iterations = 100;
    const double lambda_init = 0.001;
    const double lambda_factor = 10.0;
    const double tolerance = 1e-6;
    
    double lambda = lambda_init;
    double chi2 = std::numeric_limits<double>::max();
    double p[6];
    std::copy(p0, p0 + 6, p);
    
    // Create coordinate grids
    std::vector<double> x_grid(config.width * config.height);
    std::vector<double> y_grid(config.width * config.height);
    for (int y = 0; y < config.height; ++y) {
        for (int x = 0; x < config.width; ++x) {
            x_grid[y * config.width + x] = static_cast<double>(x);
            y_grid[y * config.width + x] = static_cast<double>(y);
        }
    }
    
    // Levenberg-Marquardt iterations with OpenMP
    for (int iter = 0; iter < max_iterations; ++iter) {
        std::vector<double> residuals(config.width * config.height);
        std::vector<std::array<double, 6>> jacobian(config.width * config.height);
        double chi2_new = 0.0;

        #ifdef _OPENMP
        #pragma omp parallel reduction(+:chi2_new)
        {
            #pragma omp for schedule(static)
            for (int i = 0; i < config.width * config.height; ++i) {
                double x = x_grid[i];
                double y = y_grid[i];
                double dx = x - p[1];
                double dy = y - p[2];
                
                // Calculate Gaussian and its derivatives
                double a = (std::cos(p[5]) * std::cos(p[5])) / (2 * p[3] * p[3]) +
                          (std::sin(p[5]) * std::sin(p[5])) / (2 * p[4] * p[4]);
                double b = -(std::sin(2 * p[5])) / (4 * p[3] * p[3]) +
                          (std::sin(2 * p[5])) / (4 * p[4] * p[4]);
                double c = (std::sin(p[5]) * std::sin(p[5])) / (2 * p[3] * p[3]) +
                          (std::cos(p[5]) * std::cos(p[5])) / (2 * p[4] * p[4]);
                
                double g = p[0] * std::exp(-(a * dx * dx + 2 * b * dx * dy + c * dy * dy));
                residuals[i] = image[i] - g;
                chi2_new += residuals[i] * residuals[i];
                
                // Calculate Jacobian
                double dg_da = g / p[0];
                double dg_dx0 = g * (2 * a * dx + 2 * b * dy);
                double dg_dy0 = g * (2 * b * dx + 2 * c * dy);
                double dg_dsx = g * dx * dx * std::cos(p[5]) * std::cos(p[5]) / (p[3] * p[3] * p[3]);
                double dg_dsy = g * dy * dy * std::sin(p[5]) * std::sin(p[5]) / (p[4] * p[4] * p[4]);
                double dg_dt = g * (dx * dx - dy * dy) * std::sin(2 * p[5]) / (2 * p[3] * p[3]);
                
                jacobian[i] = {dg_da, dg_dx0, dg_dy0, dg_dsx, dg_dsy, dg_dt};
            }
        }
        #else
        // Original implementation without OpenMP
        for (int i = 0; i < config.width * config.height; ++i) {
            double x = x_grid[i];
            double y = y_grid[i];
            double dx = x - p[1];
            double dy = y - p[2];
            
            // Calculate Gaussian and its derivatives
            double a = (std::cos(p[5]) * std::cos(p[5])) / (2 * p[3] * p[3]) +
                      (std::sin(p[5]) * std::sin(p[5])) / (2 * p[4] * p[4]);
            double b = -(std::sin(2 * p[5])) / (4 * p[3] * p[3]) +
                      (std::sin(2 * p[5])) / (4 * p[4] * p[4]);
            double c = (std::sin(p[5]) * std::sin(p[5])) / (2 * p[3] * p[3]) +
                      (std::cos(p[5]) * std::cos(p[5])) / (2 * p[4] * p[4]);
            
            double g = p[0] * std::exp(-(a * dx * dx + 2 * b * dx * dy + c * dy * dy));
            residuals[i] = image[i] - g;
            chi2_new += residuals[i] * residuals[i];
            
            // Calculate Jacobian
            double dg_da = g / p[0];
            double dg_dx0 = g * (2 * a * dx + 2 * b * dy);
            double dg_dy0 = g * (2 * b * dx + 2 * c * dy);
            double dg_dsx = g * dx * dx * std::cos(p[5]) * std::cos(p[5]) / (p[3] * p[3] * p[3]);
            double dg_dsy = g * dy * dy * std::sin(p[5]) * std::sin(p[5]) / (p[4] * p[4] * p[4]);
            double dg_dt = g * (dx * dx - dy * dy) * std::sin(2 * p[5]) / (2 * p[3] * p[3]);
            
            jacobian[i] = {dg_da, dg_dx0, dg_dy0, dg_dsx, dg_dsy, dg_dt};
        }
        #endif

        // Check convergence
        if (std::abs(chi2_new - chi2) < tolerance) {
            break;
        }
        
        // Calculate Hessian and gradient
        std::array<std::array<double, 6>, 6> hessian = {};
        std::array<double, 6> gradient = {};
        
        for (int i = 0; i < config.width * config.height; ++i) {
            for (int j = 0; j < 6; ++j) {
                gradient[j] += residuals[i] * jacobian[i][j];
                for (int k = 0; k < 6; ++k) {
                    hessian[j][k] += jacobian[i][j] * jacobian[i][k];
                }
            }
        }
        
        // Add damping term
        for (int i = 0; i < 6; ++i) {
            hessian[i][i] *= (1.0 + lambda);
        }
        
        // Solve linear system for parameter update
        std::array<double, 6> dp;
        // Simple Gaussian elimination
        for (int i = 0; i < 6; ++i) {
            double max_val = std::abs(hessian[i][i]);
            int max_row = i;
            for (int j = i + 1; j < 6; ++j) {
                if (std::abs(hessian[j][i]) > max_val) {
                    max_val = std::abs(hessian[j][i]);
                    max_row = j;
                }
            }
            
            if (max_row != i) {
                std::swap(hessian[i], hessian[max_row]);
                std::swap(gradient[i], gradient[max_row]);
            }
            
            for (int j = i + 1; j < 6; ++j) {
                double factor = hessian[j][i] / hessian[i][i];
                gradient[j] -= factor * gradient[i];
                for (int k = i; k < 6; ++k) {
                    hessian[j][k] -= factor * hessian[i][k];
                }
            }
        }
        
        // Back substitution
        for (int i = 5; i >= 0; --i) {
            dp[i] = gradient[i];
            for (int j = i + 1; j < 6; ++j) {
                dp[i] -= hessian[i][j] * dp[j];
            }
            dp[i] /= hessian[i][i];
        }
        
        // Update parameters
        double p_new[6];
        for (int i = 0; i < 6; ++i) {
            p_new[i] = p[i] + dp[i];
        }
        
        // Check if new parameters are valid
        bool valid = true;
        if (p_new[0] <= 0.0 || p_new[3] <= 0.0 || p_new[4] <= 0.0) {
            valid = false;
        }
        
        if (valid && chi2_new < chi2) {
            std::copy(p_new, p_new + 6, p);
            chi2 = chi2_new;
            lambda /= lambda_factor;
        } else {
            lambda *= lambda_factor;
        }
    }
    
    // Store results
    double gaussian_amplitude = p[0];
    double gaussian_center_x = p[1];
    double gaussian_center_y = p[2];
    double gaussian_sigma_x = p[3];
    double gaussian_sigma_y = p[4];
    double gaussian_theta = p[5];
    double gaussian_fit_error = std::sqrt(chi2 / (config.width * config.height));

    return BeamMetrics{
        moments.centroid_x, moments.centroid_y,
        moments.sigma_xx, moments.sigma_yy, moments.sigma_xy,
        moments.d4sigma_x, moments.d4sigma_y,
        moments.rotation, moments.ellipticity,
        moments.skewness_x, moments.skewness_y,
        moments.kurtosis_x, moments.kurtosis_y,
        gaussian_amplitude, gaussian_center_x, gaussian_center_y,
        gaussian_sigma_x, gaussian_sigma_y, gaussian_theta,
        gaussian_fit_error,
        moments.edge_steepness, moments.flatness, moments.rms_uniformity,
        moments.d86_radius, moments.power_in_bucket
    };
}

BeamMetrics calculateISO13694Metrics(const std::vector<float>& image, const AnalysisConfig& config, const BeamMetrics& moments) {
    // Calculate edge steepness using OpenMP
    float max_gradient = 0.0f;
    #ifdef _OPENMP
    #pragma omp parallel
    {
        float local_max = 0.0f;
        #pragma omp for schedule(static)
        for (int y = 1; y < config.height - 1; ++y) {
            for (int x = 1; x < config.width - 1; ++x) {
                float dx = image[y * config.width + x + 1] - image[y * config.width + x - 1];
                float dy = image[(y + 1) * config.width + x] - image[(y - 1) * config.width + x];
                float gradient = std::sqrt(dx * dx + dy * dy);
                local_max = std::max(local_max, gradient);
            }
        }
        #pragma omp critical
        {
            max_gradient = std::max(max_gradient, local_max);
        }
    }
    #else
    // Original implementation without OpenMP
    for (int y = 1; y < config.height - 1; ++y) {
        for (int x = 1; x < config.width - 1; ++x) {
            float dx = image[y * config.width + x + 1] - image[y * config.width + x - 1];
            float dy = image[(y + 1) * config.width + x] - image[(y - 1) * config.width + x];
            float gradient = std::sqrt(dx * dx + dy * dy);
            max_gradient = std::max(max_gradient, gradient);
        }
    }
    #endif
    float edge_steepness = max_gradient;

    // Calculate flatness and RMS uniformity using SIMD
    float min_intensity = std::numeric_limits<float>::max();
    float max_intensity = std::numeric_limits<float>::lowest();
    float sum = 0.0f;
    float sum_sq = 0.0f;

    #ifdef __AVX2__
    __m256 min_vec = _mm256_set1_ps(std::numeric_limits<float>::max());
    __m256 max_vec = _mm256_set1_ps(std::numeric_limits<float>::lowest());
    __m256 sum_vec = _mm256_setzero_ps();
    __m256 sum_sq_vec = _mm256_setzero_ps();

    int i = 0;
    for (; i + 8 <= config.width * config.height; i += 8) {
        __m256 data = _mm256_loadu_ps(&image[i]);
        min_vec = _mm256_min_ps(min_vec, data);
        max_vec = _mm256_max_ps(max_vec, data);
        sum_vec = _mm256_add_ps(sum_vec, data);
        sum_sq_vec = _mm256_add_ps(sum_sq_vec, _mm256_mul_ps(data, data));
    }

    float min_array[8], max_array[8], sum_array[8], sum_sq_array[8];
    _mm256_storeu_ps(min_array, min_vec);
    _mm256_storeu_ps(max_array, max_vec);
    _mm256_storeu_ps(sum_array, sum_vec);
    _mm256_storeu_ps(sum_sq_array, sum_sq_vec);

    for (int j = 0; j < 8; ++j) {
        min_intensity = std::min(min_intensity, min_array[j]);
        max_intensity = std::max(max_intensity, max_array[j]);
        sum += sum_array[j];
        sum_sq += sum_sq_array[j];
    }

    for (; i < config.width * config.height; ++i) {
        min_intensity = std::min(min_intensity, image[i]);
        max_intensity = std::max(max_intensity, image[i]);
        sum += image[i];
        sum_sq += image[i] * image[i];
    }
    #else
    for (int i = 0; i < config.width * config.height; ++i) {
        min_intensity = std::min(min_intensity, image[i]);
        max_intensity = std::max(max_intensity, image[i]);
        sum += image[i];
        sum_sq += image[i] * image[i];
    }
    #endif

    float flatness = min_intensity / max_intensity;
    float mean = sum / (config.width * config.height);
    float variance = sum_sq / (config.width * config.height) - mean * mean;
    float rms_uniformity = std::sqrt(variance) / mean;

    // Calculate D86 radius
    std::vector<float> sorted_intensities(image.begin(), image.end());
    std::sort(sorted_intensities.begin(), sorted_intensities.end());
    float threshold = 0.86f * sum;
    float sum_above_threshold = 0.0f;
    for (float intensity : sorted_intensities) {
        if (intensity > threshold) {
            sum_above_threshold += intensity;
        }
    }
    float d86_radius = std::sqrt(sum_above_threshold / M_PI);

    return BeamMetrics{
        moments.centroid_x, moments.centroid_y,
        moments.sigma_xx, moments.sigma_yy, moments.sigma_xy,
        moments.d4sigma_x, moments.d4sigma_y,
        moments.rotation, moments.ellipticity,
        moments.skewness_x, moments.skewness_y,
        moments.kurtosis_x, moments.kurtosis_y,
        moments.gaussian_amplitude, moments.gaussian_center_x, moments.gaussian_center_y,
        moments.gaussian_sigma_x, moments.gaussian_sigma_y, moments.gaussian_theta,
        moments.gaussian_fit_error,
        edge_steepness, flatness, rms_uniformity,
        d86_radius, moments.power_in_bucket
    };
}

BeamMetrics analyze(const void* image, const AnalysisConfig& config) {
    // Convert image to float32
    auto float_image = convertToFloat(image, config);
    
    // Apply dark frame if available
    if (config.dark_frame) {
        float_image = applyDarkFrame(float_image, *config.dark_frame);
    }
    
    // Apply gain map if available
    if (config.gain_map) {
        float_image = applyGainMap(float_image, *config.gain_map);
    }
    
    // Calculate moments
    auto metrics = calculateMoments(float_image, config);
    
    // Fit Gaussian
    metrics = fitGaussian(float_image, config, metrics);
    
    // Calculate ISO 13694 metrics
    metrics = calculateISO13694Metrics(float_image, config, metrics);
    
    return metrics;
}

double calculateM2(const std::vector<double>& z_positions, const std::vector<double>& beam_widths) {
    if (z_positions.size() != beam_widths.size() || z_positions.size() < 3) {
        throw std::invalid_argument("Invalid input data for M² calculation");
    }
    
    // Center z positions around z0 (beam waist position)
    double z0 = 0.0;
    double min_width = beam_widths[0];
    for (size_t i = 0; i < z_positions.size(); ++i) {
        if (beam_widths[i] < min_width) {
            min_width = beam_widths[i];
            z0 = z_positions[i];
        }
    }
    
    std::vector<double> z_centered(z_positions.size());
    for (size_t i = 0; i < z_positions.size(); ++i) {
        z_centered[i] = z_positions[i] - z0;
    }
    
    // Fit quadratic function w²(z) = w₀² + (M²λ/πw₀)²(z-z₀)²
    // where w₀ is the beam waist radius and λ is the wavelength
    // We'll use a simplified form: w²(z) = a + b(z-z₀)²
    
    // Calculate sums for least squares fit
    double sum_z2 = 0.0, sum_z4 = 0.0;
    double sum_w2 = 0.0, sum_z2w2 = 0.0;
    for (size_t i = 0; i < z_positions.size(); ++i) {
        double z2 = z_centered[i] * z_centered[i];
        double w2 = beam_widths[i] * beam_widths[i];
        sum_z2 += z2;
        sum_z4 += z2 * z2;
        sum_w2 += w2;
        sum_z2w2 += z2 * w2;
    }
    
    // Solve normal equations
    double n = static_cast<double>(z_positions.size());
    double det = n * sum_z4 - sum_z2 * sum_z2;
    double a = (sum_w2 * sum_z4 - sum_z2 * sum_z2w2) / det;
    double b = (n * sum_z2w2 - sum_z2 * sum_w2) / det;
    
    // Calculate M²
    // From the fit: w²(z) = w₀² + (M²λ/πw₀)²(z-z₀)²
    // We have: w²(z) = a + b(z-z₀)²
    // Therefore: M² = πw₀√b/λ
    // For simplicity, we'll assume λ = 1 (normalized units)
    double w0 = std::sqrt(a);  // Beam waist radius
    double m2 = M_PI * w0 * std::sqrt(b);
    
    return m2;
}

std::pair<double, std::vector<double>> calculateJitter(
    const std::vector<std::pair<double, double>>& centroids) {
    if (centroids.empty()) {
        throw std::invalid_argument("Empty centroid data");
    }
    
    // Calculate RMS jitter
    double sum_x = 0.0, sum_y = 0.0;
    double sum_x2 = 0.0, sum_y2 = 0.0;
    for (const auto& [x, y] : centroids) {
        sum_x += x;
        sum_y += y;
        sum_x2 += x * x;
        sum_y2 += y * y;
    }
    double mean_x = sum_x / centroids.size();
    double mean_y = sum_y / centroids.size();
    double var_x = sum_x2 / centroids.size() - mean_x * mean_x;
    double var_y = sum_y2 / centroids.size() - mean_y * mean_y;
    double rms_jitter = std::sqrt(var_x + var_y);
    
    // Calculate PSD using Welch's method
    const int n = static_cast<int>(centroids.size());
    const int nfft = 1 << static_cast<int>(std::ceil(std::log2(n)));  // Next power of 2
    const int nseg = 4;  // Number of segments
    const int overlap = nfft / 2;  // 50% overlap
    const int npsd = nfft / 2 + 1;  // Number of PSD points
    
    std::vector<double> psd(npsd, 0.0);
    std::vector<double> window(nfft);
    
    // Create Hanning window
    for (int i = 0; i < nfft; ++i) {
        window[i] = 0.5 * (1.0 - std::cos(2.0 * M_PI * i / (nfft - 1)));
    }
    
    // Process each segment
    for (int seg = 0; seg < nseg; ++seg) {
        int start = seg * (nfft - overlap);
        if (start + nfft > n) break;
        
        // Extract and window segment
        std::vector<std::complex<double>> x_fft(nfft);
        std::vector<std::complex<double>> y_fft(nfft);
        
        for (int i = 0; i < nfft; ++i) {
            double x = centroids[start + i].first - mean_x;
            double y = centroids[start + i].second - mean_y;
            x_fft[i] = std::complex<double>(x * window[i], 0.0);
            y_fft[i] = std::complex<double>(y * window[i], 0.0);
        }
        
        // Compute FFT
        fft(x_fft);
        fft(y_fft);
        
        // Accumulate PSD
        for (int i = 0; i < npsd; ++i) {
            double power = std::norm(x_fft[i]) + std::norm(y_fft[i]);
            psd[i] += power;
        }
    }
    
    // Average PSD
    double scale = 1.0 / (nseg * nfft);
    for (int i = 0; i < npsd; ++i) {
        psd[i] *= scale;
    }
    
    return {rms_jitter, psd};
}

} // namespace beam_analysis
} // namespace beamprofile 