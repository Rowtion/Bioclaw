#!/usr/bin/env python3
"""
Time Series Anomaly Detection using Aeon

This script demonstrates anomaly detection in time series using the STOMP
algorithm and other methods.

Usage:
    python anomaly_detection.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from aeon.anomaly_detection import STOMP, DWT_MLEAD, KMeansAD


def generate_synthetic_data(n_points=1000, anomaly_rate=0.05, random_state=42):
    """
    Generate synthetic time series with injected anomalies.
    
    Args:
        n_points: Number of time points
        anomaly_rate: Proportion of anomalous points
        random_state: Random seed
        
    Returns:
        time_series: The generated time series
        true_anomalies: Boolean array indicating true anomalies
    """
    np.random.seed(random_state)
    
    # Generate base signal (sine wave with trend)
    t = np.arange(n_points)
    trend = 0.01 * t
    seasonal = 10 * np.sin(2 * np.pi * t / 100)
    noise = np.random.normal(0, 1, n_points)
    
    time_series = trend + seasonal + noise
    
    # Inject anomalies
    n_anomalies = int(n_points * anomaly_rate)
    anomaly_indices = np.random.choice(n_points, n_anomalies, replace=False)
    
    true_anomalies = np.zeros(n_points, dtype=bool)
    true_anomalies[anomaly_indices] = True
    
    # Add spike anomalies
    spike_anomalies = anomaly_indices[:n_anomalies//2]
    time_series[spike_anomalies] += np.random.choice([-1, 1], len(spike_anomalies)) * \
                                     np.random.uniform(15, 25, len(spike_anomalies))
    
    # Add level shift anomalies
    shift_anomalies = anomaly_indices[n_anomalies//2:]
    for idx in shift_anomalies:
        if idx < n_points - 10:
            time_series[idx:idx+10] += np.random.choice([-1, 1]) * 10
    
    return time_series, true_anomalies


def stomp_anomaly_detection():
    """
    Demonstrate anomaly detection using STOMP (Scalable Time series
    Ordered search Matrix Profile).
    """
    print("\n" + "="*60)
    print("STOMP Anomaly Detection")
    print("="*60)
    
    # Generate data
    print("\nGenerating synthetic time series...")
    time_series, true_anomalies = generate_synthetic_data(n_points=1000)
    
    # Apply STOMP anomaly detection
    window_size = 50
    print(f"Applying STOMP with window size {window_size}...")
    
    detector = STOMP(window_size=window_size)
    anomaly_scores = detector.fit_predict(time_series)
    
    # Determine threshold (95th percentile)
    threshold = np.percentile(anomaly_scores, 95)
    predicted_anomalies = anomaly_scores > threshold
    
    # Calculate metrics
    true_positives = np.sum(predicted_anomalies & true_anomalies)
    false_positives = np.sum(predicted_anomalies & ~true_anomalies)
    false_negatives = np.sum(~predicted_anomalies & true_anomalies)
    
    precision = true_positives / (true_positives + false_positives) if (true_positives + false_positives) > 0 else 0
    recall = true_positives / (true_positives + false_negatives) if (true_positives + false_negatives) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
    
    print(f"\nAnomaly Detection Results:")
    print(f"  Threshold (95th percentile): {threshold:.2f}")
    print(f"  True anomalies: {np.sum(true_anomalies)}")
    print(f"  Detected anomalies: {np.sum(predicted_anomalies)}")
    print(f"  Precision: {precision:.3f}")
    print(f"  Recall: {recall:.3f}")
    print(f"  F1 Score: {f1:.3f}")
    
    # Plot results
    fig, axes = plt.subplots(3, 1, figsize=(14, 10))
    
    # Plot 1: Time series with anomalies
    axes[0].plot(time_series, label='Time Series', color='blue', alpha=0.7)
    anomaly_indices = np.where(true_anomalies)[0]
    axes[0].scatter(anomaly_indices, time_series[anomaly_indices], 
                    color='red', label='True Anomalies', s=50, zorder=5)
    axes[0].set_title('Time Series with True Anomalies')
    axes[0].set_xlabel('Time')
    axes[0].set_ylabel('Value')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Anomaly scores
    axes[1].plot(anomaly_scores, label='Anomaly Score', color='orange')
    axes[1].axhline(threshold, color='red', linestyle='--', label=f'Threshold ({threshold:.2f})')
    axes[1].set_title('STOMP Anomaly Scores')
    axes[1].set_xlabel('Time')
    axes[1].set_ylabel('Score')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # Plot 3: Predicted anomalies
    axes[2].plot(time_series, label='Time Series', color='blue', alpha=0.7)
    pred_indices = np.where(predicted_anomalies)[0]
    axes[2].scatter(pred_indices, time_series[pred_indices], 
                    color='green', label='Detected Anomalies', s=50, zorder=5)
    axes[2].set_title('Detected Anomalies')
    axes[2].set_xlabel('Time')
    axes[2].set_ylabel('Value')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('anomaly_detection_results.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved to 'anomaly_detection_results.png'")
    plt.close()


def segmentation_example():
    """
    Demonstrate time series segmentation using change point detection.
    """
    print("\n" + "="*60)
    print("Time Series Segmentation Example")
    print("="*60)
    
    from aeon.segmentation import ClaSPSegmenter
    
    # Generate piecewise stationary data
    np.random.seed(42)
    n_points = 500
    
    # Create segments with different characteristics
    segment1 = np.random.normal(0, 1, 100)
    segment2 = np.random.normal(5, 2, 100)
    segment3 = np.random.normal(-3, 0.5, 100)
    segment4 = np.random.normal(2, 1.5, 100)
    segment5 = np.random.normal(8, 1, 100)
    
    time_series = np.concatenate([segment1, segment2, segment3, segment4, segment5])
    true_change_points = [100, 200, 300, 400]
    
    print(f"Generated time series with {len(true_change_points)} change points")
    print(f"True change points: {true_change_points}")
    
    # Apply ClaSP segmentation
    print("\nApplying ClaSP segmentation...")
    segmenter = ClaSPSegmenter()
    predicted_change_points = segmenter.fit_predict(time_series)
    
    print(f"Detected change points: {predicted_change_points}")
    
    # Plot segmentation
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(time_series, label='Time Series', color='blue')
    
    # Mark true change points
    for cp in true_change_points:
        ax.axvline(cp, color='green', linestyle='--', alpha=0.7, label='True CP' if cp == true_change_points[0] else '')
    
    # Mark predicted change points
    for cp in predicted_change_points:
        ax.axvline(cp, color='red', linestyle=':', alpha=0.7, label='Predicted CP' if cp == predicted_change_points[0] else '')
    
    ax.set_title('Time Series Segmentation - Change Point Detection')
    ax.set_xlabel('Time')
    ax.set_ylabel('Value')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('segmentation_results.png', dpi=150, bbox_inches='tight')
    print("Plot saved to 'segmentation_results.png'")
    plt.close()


def multivariate_anomaly_detection():
    """
    Demonstrate anomaly detection on multivariate time series.
    """
    print("\n" + "="*60)
    print("Multivariate Anomaly Detection")
    print("="*60)
    
    # Generate multivariate data
    np.random.seed(42)
    n_points = 500
    n_channels = 3
    
    # Generate correlated channels
    t = np.arange(n_points)
    channel1 = np.sin(2 * np.pi * t / 50) + np.random.normal(0, 0.1, n_points)
    channel2 = np.cos(2 * np.pi * t / 50) + np.random.normal(0, 0.1, n_points)
    channel3 = 0.5 * (channel1 + channel2) + np.random.normal(0, 0.1, n_points)
    
    # Stack into multivariate series
    time_series = np.column_stack([channel1, channel2, channel3])
    
    # Inject anomalies
    anomaly_indices = [150, 300, 450]
    for idx in anomaly_indices:
        time_series[idx:idx+10, :] += np.random.uniform(3, 5, (10, n_channels))
    
    print(f"Generated {n_channels}-channel time series")
    print(f"Injected {len(anomaly_indices)} anomalous segments")
    
    # Detect anomalies on first channel
    print("\nApplying STOMP on first channel...")
    detector = STOMP(window_size=30)
    anomaly_scores = detector.fit_predict(time_series[:, 0])
    
    threshold = np.percentile(anomaly_scores, 95)
    detected = anomaly_scores > threshold
    
    print(f"Detected {np.sum(detected)} anomalous points")
    
    # Plot results
    fig, axes = plt.subplots(n_channels + 1, 1, figsize=(14, 10))
    
    for i in range(n_channels):
        axes[i].plot(time_series[:, i], label=f'Channel {i+1}')
        axes[i].set_ylabel(f'Ch {i+1}')
        axes[i].grid(True, alpha=0.3)
    
    axes[-1].plot(anomaly_scores, label='Anomaly Score', color='orange')
    axes[-1].axhline(threshold, color='red', linestyle='--', label='Threshold')
    axes[-1].set_ylabel('Score')
    axes[-1].set_xlabel('Time')
    axes[-1].legend()
    axes[-1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('multivariate_anomaly_detection.png', dpi=150, bbox_inches='tight')
    print("Plot saved to 'multivariate_anomaly_detection.png'")
    plt.close()


def main():
    """Run all anomaly detection examples."""
    print("Aeon Time Series Anomaly Detection Examples")
    print("=" * 60)
    
    try:
        # Example 1: STOMP anomaly detection
        stomp_anomaly_detection()
        
        # Example 2: Segmentation
        segmentation_example()
        
        # Example 3: Multivariate detection
        multivariate_anomaly_detection()
        
    except Exception as e:
        print(f"\nError occurred: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*60)
    print("Anomaly detection examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
