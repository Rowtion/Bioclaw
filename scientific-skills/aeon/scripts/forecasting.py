#!/usr/bin/env python3
"""
Time Series Forecasting using Aeon

This script demonstrates forecasting future values using various
forecasting algorithms including ARIMA, exponential smoothing, and
machine learning approaches.

Usage:
    python forecasting.py
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from aeon.datasets import load_airline
from aeon.forecasting.arima import ARIMA
from aeon.forecasting.exp_smoothing import ExponentialSmoothing
from aeon.forecasting.naive import NaiveForecaster
from aeon.forecasting.trend import TrendForecaster
from aeon.performance_metrics.forecasting import (
    mean_absolute_error,
    mean_absolute_percentage_error,
    mean_squared_error,
)


def load_airline_data():
    """
    Load the classic Box-Jenkins airline dataset.
    
    Returns:
        y: Time series of airline passengers
    """
    y = load_airline()
    print(f"Loaded airline dataset: {len(y)} observations")
    print(f"Frequency: {y.index.freq}")
    return y


def arima_forecasting_example():
    """
    Demonstrate ARIMA forecasting on the airline dataset.
    """
    print("\n" + "="*60)
    print("ARIMA Forecasting Example")
    print("="*60)
    
    # Load data
    y = load_airline_data()
    
    # Split into train and test
    train_size = int(len(y) * 0.8)
    y_train = y.iloc[:train_size]
    y_test = y.iloc[train_size:]
    
    print(f"\nTraining set size: {len(y_train)}")
    print(f"Test set size: {len(y_test)}")
    
    # Fit ARIMA model
    # order=(p, d, q) where:
    #   p: autoregressive order
    #   d: differencing order
    #   q: moving average order
    print("\nFitting ARIMA(1,1,1) model...")
    forecaster = ARIMA(order=(1, 1, 1))
    forecaster.fit(y_train)
    
    # Forecast
    fh = np.arange(1, len(y_test) + 1)  # Forecast horizon
    y_pred = forecaster.predict(fh)
    
    # Calculate metrics
    mae = mean_absolute_error(y_test, y_pred)
    mse = mean_squared_error(y_test, y_pred)
    mape = mean_absolute_percentage_error(y_test, y_pred)
    
    print(f"\nForecast Accuracy:")
    print(f"  MAE: {mae:.2f}")
    print(f"  MSE: {mse:.2f}")
    print(f"  MAPE: {mape:.2%}")
    
    # Plot results
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(y_train.index, y_train.values, label='Training', color='blue')
    ax.plot(y_test.index, y_test.values, label='Actual', color='green')
    ax.plot(y_test.index, y_pred.values, label='Forecast', color='red', linestyle='--')
    ax.set_title('ARIMA Forecast - Airline Passengers')
    ax.set_xlabel('Time')
    ax.set_ylabel('Passengers')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('arima_forecast.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved to 'arima_forecast.png'")
    plt.close()


def exponential_smoothing_example():
    """
    Demonstrate Exponential Smoothing forecasting.
    """
    print("\n" + "="*60)
    print("Exponential Smoothing Forecasting Example")
    print("="*60)
    
    # Load data
    y = load_airline_data()
    
    # Split data
    train_size = int(len(y) * 0.8)
    y_train = y.iloc[:train_size]
    y_test = y.iloc[train_size:]
    
    # Fit Exponential Smoothing with trend and seasonality
    print("\nFitting Exponential Smoothing model...")
    forecaster = ExponentialSmoothing(trend='add', seasonal='multiplicative', sp=12)
    forecaster.fit(y_train)
    
    # Forecast
    fh = np.arange(1, len(y_test) + 1)
    y_pred = forecaster.predict(fh)
    
    # Calculate metrics
    mae = mean_absolute_error(y_test, y_pred)
    mape = mean_absolute_percentage_error(y_test, y_pred)
    
    print(f"\nForecast Accuracy:")
    print(f"  MAE: {mae:.2f}")
    print(f"  MAPE: {mape:.2%}")
    
    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(y_train.index, y_train.values, label='Training', color='blue')
    ax.plot(y_test.index, y_test.values, label='Actual', color='green')
    ax.plot(y_test.index, y_pred.values, label='Forecast', color='purple', linestyle='--')
    ax.set_title('Exponential Smoothing Forecast - Airline Passengers')
    ax.set_xlabel('Time')
    ax.set_ylabel('Passengers')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('exp_smoothing_forecast.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved to 'exp_smoothing_forecast.png'")
    plt.close()


def compare_forecasters():
    """
    Compare multiple forecasting methods.
    """
    print("\n" + "="*60)
    print("Comparing Multiple Forecasting Methods")
    print("="*60)
    
    # Load data
    y = load_airline_data()
    
    # Split data
    train_size = int(len(y) * 0.8)
    y_train = y.iloc[:train_size]
    y_test = y.iloc[train_size:]
    fh = np.arange(1, len(y_test) + 1)
    
    # Define forecasters
    forecasters = {
        'Naive': NaiveForecaster(strategy='last'),
        'Trend': TrendForecaster(),
        'ARIMA(1,1,1)': ARIMA(order=(1, 1, 1)),
        'ExpSmoothing': ExponentialSmoothing(trend='add', seasonal='multiplicative', sp=12),
    }
    
    results = {}
    predictions = {}
    
    print("\nTraining and evaluating forecasters...")
    for name, forecaster in forecasters.items():
        print(f"\n  {name}...")
        forecaster.fit(y_train)
        y_pred = forecaster.predict(fh)
        
        mae = mean_absolute_error(y_test, y_pred)
        mape = mean_absolute_percentage_error(y_test, y_pred)
        
        results[name] = {'MAE': mae, 'MAPE': mape}
        predictions[name] = y_pred
    
    # Print comparison
    print("\n" + "-"*60)
    print("Comparison Results:")
    print("-"*60)
    print(f"{'Method':<20} {'MAE':<15} {'MAPE':<15}")
    print("-"*60)
    for name, metrics in results.items():
        print(f"{name:<20} {metrics['MAE']:<15.2f} {metrics['MAPE']:<15.2%}")
    print("-"*60)
    
    # Plot comparison
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.plot(y_train.index, y_train.values, label='Training', color='black', linewidth=2)
    ax.plot(y_test.index, y_test.values, label='Actual', color='green', linewidth=2)
    
    colors = ['blue', 'orange', 'purple', 'red']
    for i, (name, y_pred) in enumerate(predictions.items()):
        ax.plot(y_test.index, y_pred.values, label=name, 
                linestyle='--', color=colors[i % len(colors)], alpha=0.8)
    
    ax.set_title('Forecasting Methods Comparison - Airline Passengers')
    ax.set_xlabel('Time')
    ax.set_ylabel('Passengers')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('forecast_comparison.png', dpi=150, bbox_inches='tight')
    print("\nComparison plot saved to 'forecast_comparison.png'")
    plt.close()


def synthetic_forecasting_example():
    """
    Demonstrate forecasting on synthetic data with trend and seasonality.
    """
    print("\n" + "="*60)
    print("Synthetic Data Forecasting Example")
    print("="*60)
    
    # Generate synthetic data
    np.random.seed(42)
    n_points = 200
    t = np.arange(n_points)
    
    # Components
    trend = 0.05 * t
    seasonality = 10 * np.sin(2 * np.pi * t / 12)
    noise = np.random.normal(0, 2, n_points)
    
    y = pd.Series(trend + seasonality + noise + 100)
    y.index = pd.date_range(start='2020-01', periods=n_points, freq='M')
    
    print(f"Generated synthetic time series: {n_points} observations")
    print(f"  - Linear trend")
    print(f"  - Seasonal component (period=12)")
    print(f"  - Random noise")
    
    # Split data
    train_size = int(len(y) * 0.85)
    y_train = y.iloc[:train_size]
    y_test = y.iloc[train_size:]
    
    # Fit and forecast
    forecaster = ARIMA(order=(2, 1, 2))
    forecaster.fit(y_train)
    
    fh = np.arange(1, len(y_test) + 1)
    y_pred = forecaster.predict(fh)
    
    # Calculate confidence intervals (approximate)
    residuals = y_train - forecaster.predict(fh=np.arange(1, len(y_train) + 1))
    std_residuals = np.std(residuals)
    
    upper = y_pred + 1.96 * std_residuals
    lower = y_pred - 1.96 * std_residuals
    
    # Plot with confidence intervals
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(y_train.index, y_train.values, label='Training', color='blue')
    ax.plot(y_test.index, y_test.values, label='Actual', color='green')
    ax.plot(y_test.index, y_pred.values, label='Forecast', color='red', linestyle='--')
    ax.fill_between(y_test.index, lower.values, upper.values, alpha=0.2, color='red', label='95% CI')
    
    ax.set_title('ARIMA Forecast with Confidence Intervals')
    ax.set_xlabel('Time')
    ax.set_ylabel('Value')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('synthetic_forecast.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved to 'synthetic_forecast.png'")
    plt.close()
    
    # Metrics
    mape = mean_absolute_percentage_error(y_test, y_pred)
    print(f"\nForecast Accuracy: MAPE = {mape:.2%}")


def main():
    """Run all forecasting examples."""
    print("Aeon Time Series Forecasting Examples")
    print("=" * 60)
    
    try:
        # Example 1: ARIMA forecasting
        arima_forecasting_example()
        
        # Example 2: Exponential Smoothing
        exponential_smoothing_example()
        
        # Example 3: Compare multiple forecasters
        compare_forecasters()
        
        # Example 4: Synthetic data
        synthetic_forecasting_example()
        
    except Exception as e:
        print(f"\nError occurred: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "="*60)
    print("Forecasting examples completed!")
    print("="*60)


if __name__ == "__main__":
    main()
