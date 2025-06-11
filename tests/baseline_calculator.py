#!/usr/bin/env python3

import time
import numpy as np
from scipy import ndimage, stats
from scipy.optimize import curve_fit
import tifffile
import json
from pathlib import Path
import matplotlib.pyplot as plt

def gaussian_2d(xy, amplitude, x0, y0, sigma_x, sigma_y, theta):
    """2D Gaussian function for fitting"""
    x, y = xy
    if x.ndim == 1:  # If flattened arrays are passed
        x = x.reshape(-1, 1)  # Column vector
        y = y.reshape(-1, 1)  # Column vector
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    return amplitude * np.exp(-(a*(x-x0)**2 + 2*b*(x-x0)*(y-y0) + c*(y-y0)**2)).ravel()

def save_visualization(image, metrics, output_path):
    """Save visualization of the beam profile analysis results."""
    fig = plt.figure(figsize=(15, 5))
    
    # Plot 1: Original image with centroid and axes
    ax1 = fig.add_subplot(131)
    im1 = ax1.imshow(image, cmap='viridis')
    plt.colorbar(im1, ax=ax1, label='Intensity')
    
    # Plot centroid (removed label)
    ax1.plot(metrics['centroid_x'], metrics['centroid_y'], 'r+', markersize=10)
    
    # Plot major/minor axes
    theta = metrics['rotation']
    major_axis = metrics['d4sigma_x']
    minor_axis = metrics['d4sigma_y']
    
    t = np.linspace(0, 2*np.pi, 100)
    x = metrics['centroid_x'] + major_axis/2 * np.cos(t) * np.cos(theta) - minor_axis/2 * np.sin(t) * np.sin(theta)
    y = metrics['centroid_y'] + major_axis/2 * np.cos(t) * np.sin(theta) + minor_axis/2 * np.sin(t) * np.cos(theta)
    ax1.plot(x, y, 'g-')  # removed label
    
    # Plot ROI rectangle
    roi = metrics['roi']
    roi_rect = plt.Rectangle((roi['x_min'], roi['y_min']), 
                           roi['x_max'] - roi['x_min'], 
                           roi['y_max'] - roi['y_min'],
                           fill=False, color='white', linestyle='--')  # removed label
    ax1.add_patch(roi_rect)
    ax1.set_title('Beam Profile')
    
    # Plot 2: 3D surface plot of ROI
    ax2 = fig.add_subplot(132, projection='3d')
    roi_image = image[roi['y_min']:roi['y_max'], roi['x_min']:roi['x_max']]
    y, x = np.mgrid[roi['y_min']:roi['y_max'], roi['x_min']:roi['x_max']]
    surf = ax2.plot_surface(x, y, roi_image, cmap='viridis')
    plt.colorbar(surf, ax=ax2, label='Intensity')
    ax2.set_title('ROI 3D Intensity Distribution')
    
    # Plot 3: Gaussian fit (ROI only)
    ax3 = fig.add_subplot(133)
    y_fit, x_fit = np.mgrid[roi['y_min']:roi['y_max'], roi['x_min']:roi['x_max']]
    gaussian_fit = gaussian_2d((x_fit, y_fit), 
                             metrics['gaussian_amplitude'],
                             metrics['gaussian_center_x'],
                             metrics['gaussian_center_y'],
                             metrics['gaussian_sigma_x'],
                             metrics['gaussian_sigma_y'],
                             metrics['gaussian_theta'])
    im3 = ax3.imshow(gaussian_fit.reshape(roi_image.shape), cmap='viridis',
                     extent=[roi['x_min'], roi['x_max'], roi['y_max'], roi['y_min']])
    plt.colorbar(im3, ax=ax3, label='Intensity')
    ax3.set_title('Gaussian Fit (ROI)')
    
    # Add text with key metrics and legend
    metrics_text = (f"Ellipticity: {metrics['ellipticity']:.3f}\n"
                   f"RMS Uniformity: {metrics['rms_uniformity']:.3f}\n"
                   f"Fit Error: {metrics['gaussian_fit_error']:.3f}")
    
    # Create a new axes for the combined legend and metrics
    legend_ax = fig.add_axes([0.02, 0.02, 0.3, 0.1])  # [left, bottom, width, height]
    legend_ax.axis('off')
    
    # Add metrics text
    legend_ax.text(0, 0.5, metrics_text, 
                  fontsize=8, 
                  bbox=dict(facecolor='white', alpha=0.8),
                  verticalalignment='center')
    
    # Add legend items
    legend_items = [
        plt.Line2D([0], [0], marker='+', color='r', label='Centroid', markersize=10, linestyle='none'),
        plt.Line2D([0], [0], color='g', label='D4σ ellipse'),
        plt.Line2D([0], [0], color='white', label='ROI', linestyle='--')
    ]
    legend = legend_ax.legend(handles=legend_items, 
                            loc='center left',
                            bbox_to_anchor=(0.35, 0.5),
                            fontsize=8)
    legend.get_frame().set_alpha(0.8)
    legend.get_frame().set_facecolor('white')
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def calculate_metrics(image):
    """Calculate all beam profile metrics for a single image"""
    # Ensure image is float32
    image = image.astype(np.float32)
    
    # Start measure execution time for centroid and moments...
    start_time_without_gaussian_fit = time.time()
    
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
    
    # End measure execution time for centroid and moments...
    end_time_without_gaussian_fit = time.time()
    
    # Define ROI for Gaussian fitting
    # Use 3 times D4σ to ensure we capture the main beam while excluding background
    roi_half_width_x = int(np.ceil(d4sigma_x * 1.5))  # 3σ total width
    roi_half_width_y = int(np.ceil(d4sigma_y * 1.5))  # 3σ total width
    
    # Calculate ROI boundaries
    roi_x_min = max(0, int(centroid_x - roi_half_width_x))
    roi_x_max = min(image.shape[1], int(centroid_x + roi_half_width_x))
    roi_y_min = max(0, int(centroid_y - roi_half_width_y))
    roi_y_max = min(image.shape[0], int(centroid_y + roi_half_width_y))
    
    # Extract ROI
    roi_image = image[roi_y_min:roi_y_max, roi_x_min:roi_x_max]
    roi_y, roi_x = np.mgrid[roi_y_min:roi_y_max, roi_x_min:roi_x_max]
    
    # Prepare coordinates for fitting
    roi_x_flat = roi_x.ravel()
    roi_y_flat = roi_y.ravel()
    
    # Initial guess for parameters
    p0 = popt = [np.max(roi_image), centroid_x, centroid_y, d4sigma_x/4, d4sigma_y/4, rotation]
    try:
        popt, pcov = curve_fit(gaussian_2d, (roi_x_flat, roi_y_flat), roi_image.ravel(), 
                             p0=p0, 
                             method='lm',  # Explicitly specify Levenberg-Marquardt
                             maxfev=1000)  # Increase max iterations
        
        # Calculate fit error only on ROI
        gaussian_fit_roi = gaussian_2d((roi_x, roi_y), *popt).reshape(roi_image.shape)
        gaussian_fit_error = np.sqrt(np.mean((roi_image - gaussian_fit_roi)**2))
    except ZeroDivisionError:
        # print(f"Division by zero encountered during fitting")
        gaussian_fit_error = float('inf')
    except Exception as e:
        # print(f"Unexpected fitting error: {str(e)}")
        gaussian_fit_error = float('nan')
    
    # Store ROI information for visualization
    roi_info = {
        'x_min': roi_x_min,
        'x_max': roi_x_max,
        'y_min': roi_y_min,
        'y_max': roi_y_max
    }
    
    # End measure execution time for Guassian fit
    end_time_gaussian_fit = time.time()
    print(f"Time taken for Others: {end_time_without_gaussian_fit - start_time_without_gaussian_fit} seconds")
    print(f"Time taken for Gaussian fit: {end_time_gaussian_fit - end_time_without_gaussian_fit} seconds")
    
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
        'd86_radius': float(d86_radius),
        'roi': roi_info
    }

def main():
    # Create verify directory if it doesn't exist
    verify_dir = Path('../data/verify')
    verify_dir.mkdir(exist_ok=True)
    
    # Get all TIFF files in the BeamProfileData directory
    data_dir_mosa = Path('../data/mosa')    
    data_dir_oap = Path('../data/oap')
    
    tiff_files = list(data_dir_oap.glob('*.tif')) + list(data_dir_mosa.glob('*.tif'))
    
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
            
            # Check if gaussian_fit_error is infinity
            if metrics['gaussian_fit_error'] == float('inf'):
                print(f"Gaussian fit error is infinity for {tiff_file.name}")
            
            # Save visualization
            output_path = verify_dir / f"{tiff_file.stem}_analysis.png"
            save_visualization(image, metrics, output_path)
            print(f"Visualization saved to {output_path}")
            
        except Exception as e:
            print(f"Error processing {tiff_file.name}: {e}")
    
    # Save results to JSON
    output_file = 'baseline_metrics.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {output_file}")

if __name__ == '__main__':
    main() 