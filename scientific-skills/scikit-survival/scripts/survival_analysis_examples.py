#!/usr/bin/env python3
"""
scikit-survival: Survival Analysis Examples
============================================
Practical scripts for time-to-event modeling and survival analysis.
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# ==============================================================================
# Example 1: Basic Cox Proportional Hazards Model
# ==============================================================================

def example_cox_ph_basic():
    """Basic Cox Proportional Hazards survival analysis workflow."""
    from sksurv.datasets import load_breast_cancer
    from sksurv.linear_model import CoxPHSurvivalAnalysis
    from sksurv.metrics import concordance_index_ipcw
    
    # Load breast cancer survival dataset
    X, y = load_breast_cancer()
    
    print("=" * 60)
    print("Example 1: Cox Proportional Hazards Model")
    print("=" * 60)
    print(f"Dataset shape: {X.shape}")
    print(f"Features: {list(X.columns[:5])}...")
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    # Standardize features (important for Cox models)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Fit Cox PH model
    estimator = CoxPHSurvivalAnalysis()
    estimator.fit(X_train_scaled, y_train)
    
    # Predict risk scores (higher = higher risk of event)
    risk_scores = estimator.predict(X_test_scaled)
    
    # Evaluate with Uno's C-index (robust to censoring)
    c_index = concordance_index_ipcw(y_train, y_test, risk_scores)[0]
    
    print(f"\nC-index (Uno's): {c_index:.3f}")
    print(f"Model coefficients shape: {estimator.coef_.shape}")
    
    # Show top risk factors
    feature_importance = pd.DataFrame({
        'feature': X.columns,
        'coefficient': estimator.coef_
    }).sort_values('coefficient', key=abs, ascending=False)
    
    print("\nTop 5 Risk Factors:")
    print(feature_importance.head())
    
    return estimator, c_index


# ==============================================================================
# Example 2: Kaplan-Meier Survival Curves
# ==============================================================================

def example_kaplan_meier():
    """Generate and plot Kaplan-Meier survival curves."""
    from sksurv.datasets import load_gbsg2
    from sksurv.nonparametric import kaplan_meier_estimator
    import matplotlib.pyplot as plt
    
    print("\n" + "=" * 60)
    print("Example 2: Kaplan-Meier Survival Curves")
    print("=" * 60)
    
    # Load dataset
    X, y = load_gbsg2()
    
    # Overall survival curve
    time, survival_prob = kaplan_meier_estimator(y["event"], y["time"])
    
    print(f"Number of patients: {len(y)}")
    print(f"Events occurred: {y['event'].sum()}")
    print(f"Censored: {len(y) - y['event'].sum()}")
    print(f"Median survival time: {np.median(y['time']):.1f} months")
    
    # Plot survival curve by treatment group
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for treatment in X['horTh'].unique():
        mask = X['horTh'] == treatment
        time_t, surv_t = kaplan_meier_estimator(
            y["event"][mask], y["time"][mask]
        )
        ax.step(time_t, surv_t, where="post", label=f"Hormonal: {treatment}")
    
    ax.set_xlabel("Time (months)")
    ax.set_ylabel("Survival Probability")
    ax.set_title("Kaplan-Meier Survival Curves by Treatment")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('kaplan_meier_curves.png', dpi=150)
    print("\nPlot saved: kaplan_meier_curves.png")
    plt.close()


# ==============================================================================
# Example 3: Random Survival Forest
# ==============================================================================

def example_random_survival_forest():
    """Random Survival Forest for non-linear survival relationships."""
    from sksurv.datasets import load_veterans_lung_cancer
    from sksurv.ensemble import RandomSurvivalForest
    from sksurv.metrics import concordance_index_ipcw
    
    print("\n" + "=" * 60)
    print("Example 3: Random Survival Forest")
    print("=" * 60)
    
    # Load veterans lung cancer dataset
    X, y = load_veterans_lung_cancer()
    
    print(f"Dataset shape: {X.shape}")
    
    # Handle categorical variables
    X = pd.get_dummies(X, drop_first=True)
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.25, random_state=42
    )
    
    # Fit Random Survival Forest
    rsf = RandomSurvivalForest(
        n_estimators=100,
        min_samples_split=10,
        min_samples_leaf=5,
        random_state=42
    )
    rsf.fit(X_train, y_train)
    
    # Predict risk scores
    risk_scores = rsf.predict(X_test)
    c_index = concordance_index_ipcw(y_train, y_test, risk_scores)[0]
    
    print(f"C-index: {c_index:.3f}")
    print(f"Number of trees: {rsf.n_estimators}")
    
    return rsf, c_index


# ==============================================================================
# Example 4: Penalized Cox (Elastic Net) for High-Dimensional Data
# ==============================================================================

def example_coxnet():
    """Penalized Cox regression with L1/L2 regularization."""
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sklearn.model_selection import GridSearchCV
    from sksurv.metrics import as_concordance_index_ipcw_scorer
    
    print("\n" + "=" * 60)
    print("Example 4: Penalized Cox (CoxNet) with Cross-Validation")
    print("=" * 60)
    
    # Load data
    X, y = load_breast_cancer()
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    # Standardize
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # CoxNet with L1 regularization (Lasso-like)
    estimator = CoxnetSurvivalAnalysis(l1_ratio=0.9)
    
    # Grid search for alpha
    param_grid = {'alphas': [[0.01], [0.05], [0.1], [0.5]]}
    
    cv = GridSearchCV(
        estimator, 
        param_grid,
        scoring=as_concordance_index_ipcw_scorer(),
        cv=3
    )
    cv.fit(X_train_scaled, y_train)
    
    print(f"Best alpha: {cv.best_params_}")
    print(f"Best CV C-index: {cv.best_score_:.3f}")
    
    # Identify selected features
    best_model = cv.best_estimator_
    selected = np.where(best_model.coef_ != 0)[0]
    print(f"Selected features: {len(selected)} / {X.shape[1]}")
    
    return cv.best_estimator_


# ==============================================================================
# Example 5: Survival Prediction and Survival Functions
# ==============================================================================

def example_survival_predictions():
    """Generate survival probability predictions over time."""
    from sksurv.ensemble import GradientBoostingSurvivalAnalysis
    from sksurv.datasets import load_whas500
    import matplotlib.pyplot as plt
    
    print("\n" + "=" * 60)
    print("Example 5: Survival Function Predictions")
    print("=" * 60)
    
    # Load WHAS500 dataset
    X, y = load_whas500()
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    # Gradient Boosting Survival Model
    gbs = GradientBoostingSurvivalAnalysis(
        n_estimators=100,
        learning_rate=0.1,
        max_depth=3,
        random_state=42
    )
    gbs.fit(X_train, y_train)
    
    # Predict survival functions for test samples
    surv_funcs = gbs.predict_survival_function(X_test[:5])
    
    # Plot survival functions
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for i, surv_func in enumerate(surv_funcs):
        ax.step(surv_func.x, surv_func.y, where="post", label=f"Patient {i+1}")
    
    ax.set_xlabel("Time")
    ax.set_ylabel("Survival Probability")
    ax.set_title("Predicted Survival Functions (Gradient Boosting)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('survival_functions.png', dpi=150)
    print("\nPlot saved: survival_functions.png")
    plt.close()
    
    # Print survival probabilities at specific times
    times = np.percentile(y['time'], [25, 50, 75])
    print(f"\nSurvival probabilities at key time points:")
    for t in times:
        surv_at_t = [fn(t) for fn in surv_funcs]
        print(f"  Time {t:.1f}: mean survival = {np.mean(surv_at_t):.3f}")


# ==============================================================================
# Example 6: Model Comparison
# ==============================================================================

def example_model_comparison():
    """Compare multiple survival models."""
    from sksurv.linear_model import CoxPHSurvivalAnalysis
    from sksurv.ensemble import RandomSurvivalForest, GradientBoostingSurvivalAnalysis
    from sksurv.svm import FastSurvivalSVM
    from sksurv.metrics import concordance_index_ipcw
    
    print("\n" + "=" * 60)
    print("Example 6: Survival Model Comparison")
    print("=" * 60)
    
    # Load and prepare data
    X, y = load_breast_cancer()
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Define models to compare
    models = {
        'Cox PH': CoxPHSurvivalAnalysis(),
        'Random Survival Forest': RandomSurvivalForest(
            n_estimators=100, random_state=42
        ),
        'Gradient Boosting': GradientBoostingSurvivalAnalysis(
            random_state=42
        ),
        'Survival SVM': FastSurvivalSVM(random_state=42)
    }
    
    results = {}
    
    for name, model in models.items():
        print(f"\nTraining {name}...")
        
        # Fit model
        if name == 'Random Survival Forest':
            # RSF handles dataframes directly
            model.fit(X_train, y_train)
            risk_scores = model.predict(X_test)
        else:
            model.fit(X_train_scaled, y_train)
            risk_scores = model.predict(X_test_scaled)
        
        # Evaluate
        c_index = concordance_index_ipcw(y_train, y_test, risk_scores)[0]
        results[name] = c_index
        print(f"  C-index: {c_index:.3f}")
    
    # Summary
    print("\n" + "-" * 40)
    print("Model Comparison Summary:")
    print("-" * 40)
    for name, c_idx in sorted(results.items(), key=lambda x: x[1], reverse=True):
        print(f"{name:25s}: {c_idx:.3f}")
    
    return results


# ==============================================================================
# Main Execution
# ==============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("scikit-survival: Survival Analysis Examples")
    print("=" * 70)
    
    try:
        # Run all examples
        example_cox_ph_basic()
        example_kaplan_meier()
        example_random_survival_forest()
        example_coxnet()
        example_survival_predictions()
        example_model_comparison()
        
        print("\n" + "=" * 70)
        print("All examples completed successfully!")
        print("=" * 70)
        
    except ImportError as e:
        print(f"\nImport Error: {e}")
        print("Please install required packages:")
        print("  pip install scikit-survival scikit-learn pandas matplotlib")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
