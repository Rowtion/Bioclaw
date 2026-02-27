#!/usr/bin/env python3
"""
SHAP: Model Interpretability Examples
======================================
SHapley Additive exPlanations for explaining ML model predictions.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# Example 1: Basic SHAP with Tree-based Models
# ==============================================================================

def example_tree_explainer():
    """SHAP explanations for Random Forest and XGBoost models."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.datasets import load_breast_cancer
    import shap
    
    print("=" * 60)
    print("Example 1: TreeExplainer for Random Forest")
    print("=" * 60)
    
    # Load data
    data = load_breast_cancer()
    X = pd.DataFrame(data.data, columns=data.feature_names)
    y = data.target
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    print(f"Dataset: {X.shape}")
    print(f"Features: {list(X.columns[:3])}...")
    
    # Train Random Forest
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    print(f"\nModel accuracy: {model.score(X_test, y_test):.3f}")
    
    # Create TreeExplainer
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(X_test)
    
    print(f"SHAP values shape: {shap_values.values.shape}")
    
    # Global importance - Beeswarm plot
    plt.figure(figsize=(10, 8))
    shap.plots.beeswarm(shap_values, max_display=10, show=False)
    plt.title('SHAP Beeswarm: Feature Importance & Value Impact')
    plt.tight_layout()
    plt.savefig('shap_beeswarm.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_beeswarm.png")
    plt.close()
    
    # Bar plot summary
    plt.figure(figsize=(10, 6))
    shap.plots.bar(shap_values, max_display=10, show=False)
    plt.title('SHAP Bar: Mean Absolute Feature Importance')
    plt.tight_layout()
    plt.savefig('shap_bar.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_bar.png")
    plt.close()
    
    return explainer, shap_values


# ==============================================================================
# Example 2: Individual Prediction Explanation
# ==============================================================================

def example_individual_explanations():
    """Explain individual predictions with waterfall and force plots."""
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.datasets import load_diabetes
    import shap
    
    print("\n" + "=" * 60)
    print("Example 2: Individual Prediction Explanations")
    print("=" * 60)
    
    # Use diabetes dataset (regression converted to classification)
    data = load_diabetes()
    X = pd.DataFrame(data.data, columns=data.feature_names)
    y = (data.target > np.median(data.target)).astype(int)  # Binary classification
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    # Train model
    model = GradientBoostingClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    # SHAP explanations
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(X_test)
    
    # Waterfall plot for first prediction
    plt.figure(figsize=(10, 6))
    shap.plots.waterfall(shap_values[0], max_display=10, show=False)
    plt.title(f'Waterfall: Prediction Explanation (Base: {explainer.expected_value:.3f})')
    plt.tight_layout()
    plt.savefig('shap_waterfall.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_waterfall.png")
    plt.close()
    
    # Force plot for first prediction
    plt.figure(figsize=(12, 4))
    shap.plots.force(
        shap_values[0],
        matplotlib=True,
        show=False
    )
    plt.title('Force Plot: Feature Contributions')
    plt.tight_layout()
    plt.savefig('shap_force.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_force.png")
    plt.close()
    
    # Explain multiple predictions
    for i in range(min(3, len(X_test))):
        prediction = model.predict(X_test.iloc[i:i+1])[0]
        probability = model.predict_proba(X_test.iloc[i:i+1])[0]
        print(f"\nSample {i+1}:")
        print(f"  Prediction: {prediction}")
        print(f"  Probability: {probability}")
        print(f"  SHAP base value: {shap_values.base_values[i]:.3f}")
        print(f"  SHAP sum: {shap_values.values[i].sum():.3f}")
    
    return explainer, shap_values


# ==============================================================================
# Example 3: Feature Relationships with Scatter Plots
# ==============================================================================

def example_feature_relationships():
    """Analyze feature-prediction relationships and interactions."""
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import train_test_split
    from sklearn.datasets import fetch_california_housing
    import shap
    
    print("\n" + "=" * 60)
    print("Example 3: Feature Relationships and Interactions")
    print("=" * 60)
    
    # Load California housing data
    data = fetch_california_housing()
    X = pd.DataFrame(data.data, columns=data.feature_names)
    y = data.target
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    # Train model
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    print(f"R² score: {model.score(X_test, y_test):.3f}")
    
    # SHAP values
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(X_test)
    
    # Scatter plot for most important feature
    feature_names = X.columns
    importance = np.abs(shap_values.values).mean(axis=0)
    top_feature = feature_names[np.argmax(importance)]
    
    print(f"\nMost important feature: {top_feature}")
    
    # SHAP scatter plot
    plt.figure(figsize=(10, 6))
    shap.plots.scatter(shap_values[:, top_feature], show=False)
    plt.title(f'SHAP Scatter: {top_feature} Impact')
    plt.tight_layout()
    plt.savefig('shap_scatter.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_scatter.png")
    plt.close()
    
    # Scatter with color (interaction)
    second_feature = feature_names[np.argsort(importance)[-2]]
    plt.figure(figsize=(10, 6))
    shap.plots.scatter(
        shap_values[:, top_feature],
        color=shap_values[:, second_feature],
        show=False
    )
    plt.title(f'SHAP Scatter: {top_feature} colored by {second_feature}')
    plt.tight_layout()
    plt.savefig('shap_scatter_interaction.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_scatter_interaction.png")
    plt.close()
    
    # Heatmap of SHAP values across samples
    plt.figure(figsize=(12, 8))
    shap.plots.heatmap(shap_values[:50], max_display=10, show=False)  # First 50 samples
    plt.title('SHAP Heatmap: Feature Contributions Across Samples')
    plt.tight_layout()
    plt.savefig('shap_heatmap.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_heatmap.png")
    plt.close()
    
    return explainer, shap_values


# ==============================================================================
# Example 4: SHAP with Linear Models
# ==============================================================================

def example_linear_explainer():
    """SHAP explanations for linear/logistic regression models."""
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import train_test_split
    from sklearn.datasets import load_breast_cancer
    import shap
    
    print("\n" + "=" * 60)
    print("Example 4: LinearExplainer for Logistic Regression")
    print("=" * 60)
    
    # Load data
    data = load_breast_cancer()
    X = pd.DataFrame(data.data, columns=data.feature_names)
    y = data.target
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    # Standardize features (important for linear models)
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train logistic regression
    model = LogisticRegression(max_iter=1000, random_state=42)
    model.fit(X_train_scaled, y_train)
    
    print(f"Accuracy: {model.score(X_test_scaled, y_test):.3f}")
    
    # LinearExplainer
    explainer = shap.LinearExplainer(model, X_train_scaled)
    shap_values = explainer(X_test_scaled)
    
    # Summary plot
    plt.figure(figsize=(10, 8))
    shap.summary_plot(
        shap_values, X_test_scaled,
        feature_names=data.feature_names,
        max_display=10, show=False
    )
    plt.title('Linear SHAP: Feature Contributions')
    plt.tight_layout()
    plt.savefig('shap_linear_summary.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_linear_summary.png")
    plt.close()
    
    return explainer, shap_values


# ==============================================================================
# Example 5: SHAP with KernelExplainer (Model-Agnostic)
# ==============================================================================

def example_kernel_explainer():
    """Model-agnostic SHAP explanations using KernelExplainer."""
    from sklearn.svm import SVC
    from sklearn.model_selection import train_test_split
    from sklearn.datasets import load_iris
    from sklearn.preprocessing import StandardScaler
    import shap
    
    print("\n" + "=" * 60)
    print("Example 5: KernelExplainer (Model-Agnostic)")
    print("=" * 60)
    
    # Load iris dataset (simpler for SVM)
    data = load_iris()
    X = pd.DataFrame(data.data, columns=data.feature_names)
    y = data.target
    
    # Binary classification for simplicity
    mask = y < 2
    X = X[mask]
    y = y[mask]
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )
    
    # Scale features
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)
    
    # Train SVM
    model = SVC(probability=True, kernel='rbf', random_state=42)
    model.fit(X_train_scaled, y_train)
    
    print(f"SVM Accuracy: {model.score(X_test_scaled, y_test):.3f}")
    
    # KernelExplainer (model-agnostic but slower)
    # Use background data (subset of training)
    background = X_train_scaled[:50]
    
    explainer = shap.KernelExplainer(model.predict_proba, background)
    
    # Explain a few test samples (KernelExplainer is slower)
    shap_values = explainer.shap_values(X_test_scaled[:5])
    
    print(f"SHAP values shape: {np.array(shap_values).shape}")
    print("(KernelExplainer is slower but works with any model)")
    
    # Summary for class 1
    plt.figure(figsize=(10, 6))
    shap.summary_plot(
        shap_values[1], X_test_scaled[:5],
        feature_names=data.feature_names,
        show=False
    )
    plt.title('Kernel SHAP: SVM Explanations (Class 1)')
    plt.tight_layout()
    plt.savefig('shap_kernel_summary.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_kernel_summary.png")
    plt.close()
    
    return explainer, shap_values


# ==============================================================================
# Example 6: SHAP for Model Debugging
# ==============================================================================

def example_model_debugging():
    """Use SHAP to debug and validate model behavior."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.datasets import make_classification
    import shap
    
    print("\n" + "=" * 60)
    print("Example 6: Model Debugging with SHAP")
    print("=" * 60)
    
    # Create synthetic dataset with known structure
    X, y = make_classification(
        n_samples=1000, n_features=20,
        n_informative=5, n_redundant=5,
        n_repeated=0, n_classes=2,
        random_state=42
    )
    
    feature_names = [f'Feature_{i}' for i in range(20)]
    X = pd.DataFrame(X, columns=feature_names)
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    # Train model
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    print(f"Accuracy: {model.score(X_test, y_test):.3f}")
    
    # SHAP analysis
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(X_test)
    
    # Check feature importance
    importance = np.abs(shap_values.values).mean(axis=0)
    top_features = np.argsort(importance)[-5:]
    
    print("\nTop 5 features by SHAP importance:")
    for idx in reversed(top_features):
        print(f"  {feature_names[idx]}: {importance[idx]:.4f}")
    
    # Expected: Features 0-4 should be most important (they're the informative ones)
    print("\nExpected important features: Feature_0 to Feature_4")
    print("If model learned correctly, these should have high SHAP values")
    
    # Plot debugging summary
    plt.figure(figsize=(10, 8))
    shap.plots.beeswarm(shap_values, max_display=15, show=False)
    plt.title('Model Debugging: Feature Importance Check')
    plt.tight_layout()
    plt.savefig('shap_debugging.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_debugging.png")
    plt.close()
    
    return explainer, shap_values


# ==============================================================================
# Example 7: Cohort Comparison with SHAP
# ==============================================================================

def example_cohort_comparison():
    """Compare SHAP values across different data cohorts."""
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.datasets import load_breast_cancer
    import shap
    
    print("\n" + "=" * 60)
    print("Example 7: Cohort Comparison")
    print("=" * 60)
    
    # Load data
    data = load_breast_cancer()
    X = pd.DataFrame(data.data, columns=data.feature_names)
    y = data.target
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.3, random_state=42
    )
    
    # Train model
    model = RandomForestClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    # SHAP values
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(X_test)
    
    # Split test set into two cohorts based on mean radius
    median_radius = X_test['mean radius'].median()
    cohort1_mask = X_test['mean radius'] <= median_radius
    cohort2_mask = X_test['mean radius'] > median_radius
    
    print(f"Cohort 1 (small radius): {cohort1_mask.sum()} samples")
    print(f"Cohort 2 (large radius): {cohort2_mask.sum()} samples")
    
    # Compare bar plots
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    plt.sca(axes[0])
    shap.plots.bar(shap_values[cohort1_mask], max_display=8, show=False)
    axes[0].set_title('Cohort 1: Small Tumor Radius')
    
    plt.sca(axes[1])
    shap.plots.bar(shap_values[cohort2_mask], max_display=8, show=False)
    axes[1].set_title('Cohort 2: Large Tumor Radius')
    
    plt.tight_layout()
    plt.savefig('shap_cohort_comparison.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_cohort_comparison.png")
    plt.close()
    
    return explainer, shap_values


# ==============================================================================
# Example 8: SHAP Decision Plot
# ==============================================================================

def example_decision_plot():
    """Visualize decision paths with SHAP decision plots."""
    from sklearn.ensemble import GradientBoostingClassifier
    from sklearn.model_selection import train_test_split
    from sklearn.datasets import load_breast_cancer
    import shap
    
    print("\n" + "=" * 60)
    print("Example 8: SHAP Decision Plot")
    print("=" * 60)
    
    # Load data
    data = load_breast_cancer()
    X = pd.DataFrame(data.data, columns=data.feature_names)
    y = data.target
    
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42
    )
    
    # Train model
    model = GradientBoostingClassifier(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    # SHAP
    explainer = shap.TreeExplainer(model)
    shap_values = explainer(X_test)
    
    # Decision plot for a few samples
    plt.figure(figsize=(10, 10))
    shap.decision_plot(
        explainer.expected_value,
        shap_values.values[:5],
        X_test.iloc[:5],
        feature_names=data.feature_names,
        show=False
    )
    plt.title('SHAP Decision Plot: Prediction Paths')
    plt.tight_layout()
    plt.savefig('shap_decision.png', dpi=150, bbox_inches='tight')
    print("Saved: shap_decision.png")
    plt.close()
    
    return explainer, shap_values


# ==============================================================================
# Main Execution
# ==============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("SHAP: Model Interpretability Examples")
    print("=" * 70)
    
    try:
        example_tree_explainer()
        example_individual_explanations()
        example_feature_relationships()
        example_linear_explainer()
        example_kernel_explainer()
        example_model_debugging()
        example_cohort_comparison()
        example_decision_plot()
        
        print("\n" + "=" * 70)
        print("All examples completed successfully!")
        print("Generated plots:")
        print("  - shap_beeswarm.png")
        print("  - shap_bar.png")
        print("  - shap_waterfall.png")
        print("  - shap_force.png")
        print("  - shap_scatter.png")
        print("  - shap_scatter_interaction.png")
        print("  - shap_heatmap.png")
        print("  - shap_linear_summary.png")
        print("  - shap_kernel_summary.png")
        print("  - shap_debugging.png")
        print("  - shap_cohort_comparison.png")
        print("  - shap_decision.png")
        print("=" * 70)
        
    except ImportError as e:
        print(f"\nImport Error: {e}")
        print("Please install required packages:")
        print("  pip install shap scikit-learn")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
