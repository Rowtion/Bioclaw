#!/usr/bin/env python3
"""
Time Series Classification using Aeon

This script demonstrates time series classification using ROCKET and MiniRocket
transformers with various classifiers.

Usage:
    python time_series_classification.py
"""

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import RidgeClassifierCV
from sklearn.metrics import accuracy_score, classification_report
from sklearn.model_selection import train_test_split

# Aeon imports
from aeon.classification.convolution_based import RocketClassifier, MiniRocketClassifier
from aeon.datasets import load_classification
from aeon.transformations.collection.convolution_based import (
    MiniRocket,
    MultiRocket,
    Rocket,
)


def load_example_data(dataset_name="GunPoint"):
    """
    Load example time series classification dataset.
    
    Args:
        dataset_name: Name of the UCR/UEA dataset
        
    Returns:
        X_train, y_train, X_test, y_test
    """
    print(f"Loading dataset: {dataset_name}")
    X_train, y_train = load_classification(dataset_name, split="train")
    X_test, y_test = load_classification(dataset_name, split="test")
    
    print(f"Training set shape: {X_train.shape}")
    print(f"Test set shape: {X_test.shape}")
    print(f"Number of classes: {len(np.unique(y_train))}")
    
    return X_train, y_train, X_test, y_test


def rocket_classification_example():
    """
    Demonstrate ROCKET classifier for time series classification.
    ROCKET is fast and achieves state-of-the-art accuracy.
    """
    print("\n" + "="*60)
    print("ROCKET Classifier Example")
    print("="*60)
    
    # Load data
    X_train, y_train, X_test, y_test = load_example_data("GunPoint")
    
    # Create and train ROCKET classifier
    # n_kernels: number of random convolutional kernels (default 10000)
    classifier = RocketClassifier(n_kernels=10000, random_state=42)
    
    print("\nTraining ROCKET classifier...")
    classifier.fit(X_train, y_train)
    
    # Make predictions
    print("Making predictions...")
    y_pred = classifier.predict(X_test)
    
    # Evaluate
    accuracy = accuracy_score(y_test, y_pred)
    print(f"\nAccuracy: {accuracy:.4f}")
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))


def minirocket_classification_example():
    """
    Demonstrate MiniRocket classifier - faster version of ROCKET.
    """
    print("\n" + "="*60)
    print("MiniRocket Classifier Example")
    print("="*60)
    
    # Load data
    X_train, y_train, X_test, y_test = load_example_data("GunPoint")
    
    # Create and train MiniRocket classifier
    # Much faster than ROCKET with comparable accuracy
    classifier = MiniRocketClassifier(random_state=42)
    
    print("\nTraining MiniRocket classifier...")
    classifier.fit(X_train, y_train)
    
    # Make predictions
    print("Making predictions...")
    y_pred = classifier.predict(X_test)
    
    # Evaluate
    accuracy = accuracy_score(y_test, y_pred)
    print(f"\nAccuracy: {accuracy:.4f}")
    print("\nClassification Report:")
    print(classification_report(y_test, y_pred))


def rocket_feature_extraction_example():
    """
    Demonstrate using ROCKET for feature extraction followed by
    traditional machine learning classifiers.
    """
    print("\n" + "="*60)
    print("ROCKET Feature Extraction + Ridge Classifier")
    print("="*60)
    
    # Load data
    X_train, y_train, X_test, y_test = load_example_data("GunPoint")
    
    # Transform data using ROCKET
    print("\nExtracting ROCKET features...")
    rocket = Rocket(n_kernels=10000, random_state=42)
    X_train_features = rocket.fit_transform(X_train)
    X_test_features = rocket.transform(X_test)
    
    print(f"Original shape: {X_train.shape}")
    print(f"Feature shape: {X_train_features.shape}")
    
    # Train Ridge classifier on features
    print("\nTraining Ridge classifier on ROCKET features...")
    classifier = RidgeClassifierCV(alphas=np.logspace(-3, 3, 10))
    classifier.fit(X_train_features, y_train)
    
    # Predict and evaluate
    y_pred = classifier.predict(X_test_features)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"\nAccuracy: {accuracy:.4f}")


def multivariate_classification_example():
    """
    Demonstrate classification on multivariate time series.
    """
    print("\n" + "="*60)
    print("Multivariate Time Series Classification")
    print("="*60)
    
    # Create synthetic multivariate data
    n_samples = 100
    n_channels = 3
    n_timepoints = 50
    
    # Generate synthetic data
    np.random.seed(42)
    X = np.random.randn(n_samples, n_channels, n_timepoints)
    
    # Create labels based on simple pattern
    y = np.array([
        "Class_A" if np.sum(x[0]) > 0 else "Class_B"
        for x in X
    ])
    
    print(f"Data shape: {X.shape}")
    print(f"Channels: {n_channels}, Timepoints: {n_timepoints}")
    print(f"Classes: {np.unique(y)}")
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    # Use MiniRocket for multivariate data
    print("\nTraining MiniRocket classifier...")
    classifier = MiniRocketClassifier(random_state=42)
    classifier.fit(X_train, y_train)
    
    y_pred = classifier.predict(X_test)
    accuracy = accuracy_score(y_test, y_pred)
    print(f"\nAccuracy: {accuracy:.4f}")


def main():
    """Run all classification examples."""
    print("Aeon Time Series Classification Examples")
    print("=" * 60)
    
    try:
        # Example 1: ROCKET classifier
        rocket_classification_example()
        
        # Example 2: MiniRocket classifier
        minirocket_classification_example()
        
        # Example 3: ROCKET + custom classifier
        rocket_feature_extraction_example()
        
        # Example 4: Multivariate data
        multivariate_classification_example()
        
    except Exception as e:
        print(f"\nError occurred: {e}")
        print("Note: Some examples require internet connection to download datasets.")
    
    print("\n" + "="*60)
    print("Classification examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
