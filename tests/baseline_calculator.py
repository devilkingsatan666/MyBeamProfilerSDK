#!/usr/bin/env python3

import os
import numpy as np
from scipy import ndimage, stats
from scipy.optimize import curve_fit
import tifffile
import json
from pathlib import Path

def gaussian_2d(xy, amplitude, x0, y0, sigma_x, sigma_y, theta):
    """2D Gaussian function for fitting"""
    x, y = xy
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    return amplitude * np.exp(-(a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2))

def calculate_metrics(image):
    """Calculate all beam profile metrics for a single image"""
    # Ensure image is float32
    image = image.astype(np.float32)
    
    # Calculate centroid and second moments
    y, x = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    total = np.sum(image)
    centroid_x = np.sum(x * image) / total
    centroid_y = np.sum(y * image) / total
    
    # Center coordinates
    x_centered = x - centroid_x
    y_centered = y - centroid_y
    
    # Second moments
    sigma_xx = np.sum(x_centered**2 * image) / total
    sigma_yy = np.sum(y_centered**2 * image) / total
    sigma_xy = np.sum(x_centered * y_centered * image) / total
    
    # D4σ widths
    d4sigma_x = 4 * np.sqrt(sigma_xx)
    d4sigma_y = 4 * np.sqrt(sigma_yy)
    
    # Rotation and ellipticity
    rotation = 0.5 * np.arctan2(2*sigma_xy, sigma_xx - sigma_yy)
    ellipticity = np.sqrt((sigma_xx - sigma_yy)**2 + 4*sigma_xy**2) / (sigma_xx + sigma_yy)
    
    # Higher moments
    skewness_x = np.sum(x_centered**3 * image) / (total * sigma_xx**1.5)
    skewness_y = np.sum(y_centered**3 * image) / (total * sigma_yy**1.5)
    kurtosis_x = np.sum(x_centered**4 * image) / (total * sigma_xx**2) - 3
    kurtosis_y = np.sum(y_centered**4 * image) / (total * sigma_yy**2) - 3
    
    # Gaussian fit
    p0 = [np.max(image), centroid_x, centroid_y, d4sigma_x/4, d4sigma_y/4, rotation]
    try:
        popt, pcov = curve_fit(gaussian_2d, (x, y), image.ravel(), p0=p0)
        gaussian_fit = gaussian_2d((x, y), *popt).reshape(image.shape)
        gaussian_fit_error = np.sqrt(np.mean((image - gaussian_fit)**2))
    except:
        popt = p0
        gaussian_fit_error = float('inf')
    
    # ISO 13694 metrics
    # Edge steepness (simplified)
    edge_steepness = np.max(np.abs(np.gradient(image)))
    
    # Flatness (simplified)
    flatness = np.min(image) / np.max(image)
    
    # RMS uniformity
    rms_uniformity = np.std(image) / np.mean(image)
    
    # D86 radius
    sorted_intensities = np.sort(image.ravel())
    threshold = 0.86 * np.sum(sorted_intensities)
    d86_radius = np.sqrt(np.sum(sorted_intensities[sorted_intensities > threshold]) / np.pi)
    
    return {
        'centroid_x': float(centroid_x),
        'centroid_y': float(centroid_y),
        'sigma_xx': float(sigma_xx),
        'sigma_yy': float(sigma_yy),
        'sigma_xy': float(sigma_xy),
        'd4sigma_x': float(d4sigma_x),
        'd4sigma_y': float(d4sigma_y),
        'rotation': float(rotation),
        'ellipticity': float(ellipticity),
        'skewness_x': float(skewness_x),
        'skewness_y': float(skewness_y),
        'kurtosis_x': float(kurtosis_x),
        'kurtosis_y': float(kurtosis_y),
        'gaussian_amplitude': float(popt[0]),
        'gaussian_center_x': float(popt[1]),
        'gaussian_center_y': float(popt[2]),
        'gaussian_sigma_x': float(popt[3]),
        'gaussian_sigma_y': float(popt[4]),
        'gaussian_theta': float(popt[5]),
        'gaussian_fit_error': float(gaussian_fit_error),
        'edge_steepness': float(edge_steepness),
        'flatness': float(flatness),
        'rms_uniformity': float(rms_uniformity),
        'd86_radius': float(d86_radius)
    }

def main():
    # Get all TIFF files in the BeamProfileData directory
    data_dir = Path('../BeamProfileData')
    tiff_files = list(data_dir.glob('*.tif'))
    
    # Calculate metrics for each file
    results = {}
    for tiff_file in tiff_files:
        print(f"Processing {tiff_file.name}...")
        try:
            image = tifffile.imread(tiff_file)
            # Convert to grayscale if RGB
            if len(image.shape) == 3:
                image = np.mean(image, axis=2)
            metrics = calculate_metrics(image)
            results[tiff_file.name] = metrics
        except Exception as e:
            print(f"Error processing {tiff_file.name}: {e}")
    
    # Save results to JSON
    output_file = 'baseline_metrics.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")

if __name__ == '__main__':
    main() 