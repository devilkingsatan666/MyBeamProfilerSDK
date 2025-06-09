# Contributing to Beam Profiler SDK

Thank you for your interest in contributing to the Beam Profiler SDK! This document provides guidelines and instructions for contributing.

## Code of Conduct

Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

## Development Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Run tests and ensure they pass
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## Pull Request Process

1. Update the README.md with details of changes if needed
2. Update the documentation if you've changed the API
3. Ensure all tests pass
4. The PR will be merged once you have the sign-off of at least one other developer

## Development Setup

### Prerequisites

- CMake 3.15 or later
- C++17 compatible compiler
- Git
- Python 3.7+ (for tests)
- libtiff development files

### Building

```bash
# Clone the repository
git clone https://github.com/yourusername/beam_profiler_sdk.git
cd beam_profiler_sdk

# Create build directory
mkdir build && cd build

# Configure with development options
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DBUILD_TESTS=ON \
      -DBUILD_EXAMPLES=ON \
      -DBUILD_DOCS=ON \
      ..

# Build
cmake --build .
```

### Running Tests

```bash
# Run C++ tests
ctest

# Run Python baseline tests
cd ../tests
python baseline_calculator.py
```

## Coding Standards

### C++ Code Style

- Follow the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html)
- Use clang-format for formatting
- Maximum line length: 100 characters
- Use spaces for indentation (4 spaces)
- Use meaningful variable and function names
- Add comments for complex logic

### C Code Style

- Follow the [Linux kernel coding style](https://www.kernel.org/doc/html/latest/process/coding-style.html)
- Use clang-format for formatting
- Maximum line length: 80 characters
- Use spaces for indentation (4 spaces)
- Use meaningful variable and function names
- Add comments for complex logic

### Documentation

- Use Doxygen-style comments for all public APIs
- Keep documentation up to date with code changes
- Include examples in documentation
- Document all parameters and return values
- Document error conditions and handling

### Testing

- Write unit tests for new features
- Maintain or improve test coverage
- Include test data with tests
- Document test requirements and setup

## Performance Guidelines

- Profile code changes for performance impact
- Use SIMD instructions where appropriate
- Optimize for cache locality
- Minimize memory allocations
- Use appropriate data structures
- Document performance characteristics

## Error Handling

- Use appropriate error codes
- Provide meaningful error messages
- Handle all error conditions
- Document error handling behavior
- Include error handling in tests

## Memory Management

- Use RAII in C++ code
- Properly free resources in C code
- Avoid memory leaks
- Document memory ownership
- Use smart pointers where appropriate

## Security

- Validate all input
- Avoid buffer overflows
- Use secure coding practices
- Document security considerations
- Include security tests

## Commit Messages

- Use present tense ("Add feature" not "Added feature")
- Use imperative mood ("Move cursor to..." not "Moves cursor to...")
- Limit the first line to 72 characters or less
- Reference issues and pull requests liberally
- Consider starting the commit message with an applicable emoji:
  - 🎨 `:art:` when improving the format/structure of the code
  - 🐎 `:racehorse:` when improving performance
  - 🚱 `:non-potable_water:` when plugging memory leaks
  - 📝 `:memo:` when writing docs
  - 🐛 `:bug:` when fixing a bug
  - 🔥 `:fire:` when removing code or files
  - 💚 `:green_heart:` when fixing the CI build
  - ✅ `:white_check_mark:` when adding tests
  - 🔒 `:lock:` when dealing with security
  - ⬆️ `:arrow_up:` when upgrading dependencies
  - ⬇️ `:arrow_down:` when downgrading dependencies

## License

By contributing to this project, you agree that your contributions will be licensed under the project's MIT License. 