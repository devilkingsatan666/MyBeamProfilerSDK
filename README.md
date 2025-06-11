# Beam Profiler SDK

A high-performance C/C++ library for laser beam profile analysis, implementing ISO 11146 and ISO 13694 standards.

## Features

- Support for multiple image formats (8/12/16-bit monochrome, RGB)
- Dark frame subtraction and gain map correction
- Second-moment analysis (ISO 11146)
- Gaussian beam fitting
- M² calculation
- Pointing jitter analysis
- ISO 13694 metrics (edge steepness, flatness, RMS uniformity)
- SIMD/OpenMP optimizations
- C and C++ APIs

## Requirements

- C++17 or later
- CMake 3.15 or later
- vcpkg package manager
- Visual Studio 2019 or later (for Windows)
- GCC 7+ or Clang 6+ (for Linux/macOS)

## Building from Source

### Prerequisites

1. Install vcpkg if you haven't already:
```bash
git clone https://github.com/Microsoft/vcpkg.git
cd vcpkg
./bootstrap-vcpkg.bat  # On Windows
./bootstrap-vcpkg.sh   # On Linux/macOS
```

2. Integrate vcpkg with your system:
```bash
vcpkg integrate install
```

3. Install required libraries:
```bash
# Core libraries
vcpkg install libtiff:x64-windows      # Use :x64-linux or :x64-osx for other platforms
vcpkg install nlohmann-json:x64-windows # For configuration and test data

# Libraries required for testing
vcpkg install gtest:x64-windows        # For unit testing
```

### Building the Project

1. Clone the repository:
```bash
git clone https://github.com/devilkingsatan666/MyBeamProfilerSDK.git
cd beam_profiler_sdk
```

2. Create and navigate to the build directory:
```bash
mkdir build
cd build
```

3. Configure with CMake:
```bash
# For Windows
cmake .. -DCMAKE_TOOLCHAIN_FILE=[path_to_vcpkg]/scripts/buildsystems/vcpkg.cmake

# For Linux/macOS
cmake .. -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
```

4. Build the project:
```bash
cmake --build . --config Release
```

This will build:
- The C++ library (beamprofile_cpp)
- The C library (beamprofile_c)
- Example programs
- Test suite

### Running Tests

After building, you can run the tests:
```bash
ctest -C Release
```

### Project Structure

- `src/` - Source files for the core library
- `examples/` - Example programs demonstrating usage
- `tests/` - Test suite
- `build/` - Build output directory

## Usage

### C API Example

```c
#include <beamprofile/beamprofile.h>
#include <stdio.h>

int main() {
    // Create configuration
    BeamProfileConfig* config = beamprofile_create_config(640, 480, BEAMPROFILE_MONO8);
    if (!config) {
        fprintf(stderr, "Failed to create configuration\n");
        return 1;
    }

    // Load image data
    uint8_t* image_data = load_image("beam.tif");  // Your image loading code here

    // Analyze beam profile
    BeamProfileMetrics* metrics = beamprofile_analyze(image_data, config);
    if (!metrics) {
        fprintf(stderr, "Analysis failed\n");
        beamprofile_destroy_config(config);
        return 1;
    }

    // Print all available metrics
    printf("Basic Metrics:\n");
    printf("  Centroid: (%.2f, %.2f) pixels\n", metrics->centroid_x, metrics->centroid_y);
    printf("  D4σ Width: %.2f x %.2f pixels\n", metrics->d4sigma_x, metrics->d4sigma_y);
    printf("  Orientation: %.2f degrees\n", metrics->orientation);
    printf("  Peak Intensity: %.2f\n", metrics->peak_intensity);
    printf("  Total Power: %.2f\n", metrics->total_power);

    printf("\nShape Analysis:\n");
    printf("  Ellipticity: %.3f\n", metrics->ellipticity);
    printf("  Circularity: %.3f\n", metrics->circularity);
    printf("  Aspect Ratio: %.3f\n", metrics->aspect_ratio);

    printf("\nStatistical Moments:\n");
    printf("  Skewness: (%.3f, %.3f)\n", metrics->skewness_x, metrics->skewness_y);
    printf("  Kurtosis: (%.3f, %.3f)\n", metrics->kurtosis_x, metrics->kurtosis_y);

    printf("\nBeam Quality:\n");
    printf("  M²: (%.3f, %.3f)\n", metrics->m2_x, metrics->m2_y);
    printf("  Beam Quality Factor: %.3f\n", metrics->beam_quality);
    printf("  Divergence: (%.3f, %.3f) mrad\n", metrics->divergence_x, metrics->divergence_y);

    printf("\nPower Analysis:\n");
    printf("  Power in Bucket (86%%): %.2f%%\n", metrics->power_in_bucket_86 * 100.0);
    printf("  Effective Area: %.2f pixels²\n", metrics->effective_area);
    printf("  Flatness: %.3f\n", metrics->flatness);
    printf("  Roughness: %.3f\n", metrics->roughness);

    // Cleanup
    beamprofile_destroy_metrics(metrics);
    beamprofile_destroy_config(config);
    return 0;
}
```

### C++ API Example

```cpp
#include <beamprofile/beamprofile.hpp>
#include <iostream>
#include <iomanip>

int main() {
    // Create configuration
    beamprofile::AnalysisConfig config{
        640, 480,
        beamprofile::ImageFormat::MONO8
    };

    // Load image data
    std::vector<uint8_t> image_data = load_image("beam.tif");  // Your image loading code here

    // Analyze beam profile
    auto metrics = beamprofile::beam_analysis::analyze(image_data.data(), config);

    // Set output formatting
    std::cout << std::fixed << std::setprecision(3);

    // Print all available metrics
    std::cout << "Basic Metrics:\n";
    std::cout << "  Centroid: (" << metrics.centroid_x << ", " << metrics.centroid_y << ") pixels\n";
    std::cout << "  D4σ Width: " << metrics.d4sigma_x << " x " << metrics.d4sigma_y << " pixels\n";
    std::cout << "  Orientation: " << metrics.orientation << " degrees\n";
    std::cout << "  Peak Intensity: " << metrics.peak_intensity << "\n";
    std::cout << "  Total Power: " << metrics.total_power << "\n";

    std::cout << "\nShape Analysis:\n";
    std::cout << "  Ellipticity: " << metrics.ellipticity << "\n";
    std::cout << "  Circularity: " << metrics.circularity << "\n";
    std::cout << "  Aspect Ratio: " << metrics.aspect_ratio << "\n";

    std::cout << "\nStatistical Moments:\n";
    std::cout << "  Skewness: (" << metrics.skewness_x << ", " << metrics.skewness_y << ")\n";
    std::cout << "  Kurtosis: (" << metrics.kurtosis_x << ", " << metrics.kurtosis_y << ")\n";

    std::cout << "\nBeam Quality:\n";
    std::cout << "  M²: (" << metrics.m2_x << ", " << metrics.m2_y << ")\n";
    std::cout << "  Beam Quality Factor: " << metrics.beam_quality << "\n";
    std::cout << "  Divergence: (" << metrics.divergence_x << ", " << metrics.divergence_y << ") mrad\n";

    std::cout << "\nPower Analysis:\n";
    std::cout << "  Power in Bucket (86%): " << metrics.power_in_bucket_86 * 100.0 << "%\n";
    std::cout << "  Effective Area: " << metrics.effective_area << " pixels²\n";
    std::cout << "  Flatness: " << metrics.flatness << "\n";
    std::cout << "  Roughness: " << metrics.roughness << "\n";

    return 0;
}
```

## API Documentation

### C API

The C API provides a set of functions for beam profile analysis:

- `beamprofile_create_config()`: Create a new configuration
- `beamprofile_destroy_config()`: Destroy a configuration
- `beamprofile_set_dark_frame()`: Set dark frame for subtraction
- `beamprofile_set_gain_map()`: Set gain map for correction
- `beamprofile_analyze()`: Analyze beam profile
- `beamprofile_destroy_metrics()`: Destroy analysis results
- `beamprofile_calculate_m2()`: Calculate M² from multiple measurements
- `beamprofile_calculate_jitter()`: Calculate pointing jitter metrics
- `beamprofile_destroy_jitter()`: Destroy jitter analysis results

### C++ API

The C++ API provides a more modern interface with RAII and type safety:

- `beamprofile::AnalysisConfig`: Configuration structure
- `beamprofile::BeamMetrics`: Analysis results structure
- `beamprofile::beam_analysis::analyze()`: Main analysis function
- `beamprofile::beam_analysis::calculateM2()`: M² calculation
- `beamprofile::beam_analysis::calculateJitter()`: Jitter analysis

## Mathematical Background

### Beam Width and Center Calculations

The beam width and center calculations are performed according to ISO 11146 using the second moment method:

1. **Centroid (First Moment):**
   $$\bar{x} = \frac{\sum(x_i \cdot I_i)}{\sum I_i}$$
   $$\bar{y} = \frac{\sum(y_i \cdot I_i)}{\sum I_i}$$
   where $(x_i, y_i)$ are pixel coordinates and $I_i$ is the intensity at that pixel.

2. **Second Moments:**
   $$\sigma_x^2 = \frac{\sum((x_i - \bar{x})^2 \cdot I_i)}{\sum I_i}$$
   $$\sigma_y^2 = \frac{\sum((y_i - \bar{y})^2 \cdot I_i)}{\sum I_i}$$

3. **D4σ Beam Width:**
   $$D4\sigma_{x} = 4{\sqrt{\sigma_{x}^2}}$$
   $$D4\sigma_{y} = 4{\sqrt{\sigma_{y}^2}}$$
   This represents the width containing ~95% of the beam's power.

### Statistical Moments

1. **Skewness:**
   $$\gamma_x = \frac{\sum((x_i - \bar{x})^3 \cdot I_i)}{\sum I_i \cdot \sigma_x^3}$$
   $$\gamma_y = \frac{\sum((y_i - \bar{y})^3 \cdot I_i)}{\sum I_i \cdot \sigma_y^3}$$
   Measures beam asymmetry along X and Y axes.

2. **Kurtosis:**
   $$\beta_x = \frac{\sum((x_i - \bar{x})^4 \cdot I_i)}{\sum I_i \cdot \sigma_x^4}$$
   $$\beta_y = \frac{\sum((y_i - \bar{y})^4 \cdot I_i)}{\sum I_i \cdot \sigma_y^4}$$
   Measures how peaked or flat the beam profile is compared to a Gaussian.

### Power Metrics (ISO 13694)

1. **Power in Bucket:**
   $$P(r) = \iint_{x^2 + y^2 \leq r^2} I(x,y) \, dx \, dy$$
   Integrated power within radius $r$ from centroid.

2. **Fractional Power:**
   $$F(r) = \frac{P(r)}{P(\infty)}$$
   Ratio of power within radius $r$ to total power.

### Beam Quality (M²)

M² calculation requires measurements at multiple z positions:

1. **Beam Width vs. Position:**
   $$w^2(z) = w_0^2 + \theta^2(z - z_0)^2$$
   where:
   - $w_0$ is the waist size
   - $\theta$ is the far-field divergence
   - $z_0$ is the waist position

2. **M² Calculation:**
   $$M^2 = \frac{\pi w_0 \theta}{4\lambda}$$
   where $\lambda$ is the wavelength.

### Implementation Notes

- All calculations are performed using X/Y axes rather than principal axes for:
  - Direct correlation with physical setup
  - Computational efficiency
  - Industry standard compatibility
  - Simpler interpretation
- Intensity values are background-corrected before calculations
- Integration uses efficient numerical methods optimized for discrete pixel data
- Error propagation is considered in all calculations

## Performance

The library is optimized for performance using:

- SIMD instructions (AVX2) for vectorized operations
- OpenMP for parallel processing
- Efficient memory management
- Optimized algorithms for beam analysis

## Contributing

Contributions are welcome! Please read our [Contributing Guidelines](CONTRIBUTING.md) for details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributors

- Wu-Cheng, Chiang

## References

- ISO 11146 and ISO 13694 standards
- OpenCV for image processing inspiration
- libtiff for image I/O
- OpenMP for parallel processing 