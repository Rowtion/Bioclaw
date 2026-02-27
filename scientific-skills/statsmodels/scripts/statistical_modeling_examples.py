#!/usr/bin/env python3
"""
Statsmodels: Statistical Modeling Examples
==========================================
Statistical models, econometrics, and rigorous inference.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ==============================================================================
# Example 1: Ordinary Least Squares (OLS) Regression
# ==============================================================================

def example_ols_regression():
    """Basic OLS regression with diagnostics."""
    import statsmodels.api as sm
    from statsmodels.stats.diagnostic import het_breuschpagan
    from statsmodels.stats.outliers_influence import variance_inflation_factor
    
    print("=" * 60)
    print("Example 1: OLS Linear Regression")
    print("=" * 60)
    
    # Generate synthetic data
    np.random.seed(42)
    n = 200
    X1 = np.random.randn(n)
    X2 = np.random.randn(n)
    X3 = 0.5 * X1 + 0.3 * X2 + np.random.randn(n) * 0.5  # Correlated with X1, X2
    
    # True model: y = 2 + 3*X1 - 2*X2 + 0.5*X3 + error
    y = 2 + 3*X1 - 2*X2 + 0.5*X3 + np.random.randn(n) * 2
    
    # Create DataFrame
    df = pd.DataFrame({'X1': X1, 'X2': X2, 'X3': X3, 'y': y})
    
    # Prepare X with constant (intercept)
    X = sm.add_constant(df[['X1', 'X2', 'X3']])
    
    # Fit OLS model
    model = sm.OLS(df['y'], X)
    results = model.fit()
    
    # Print comprehensive results
    print(results.summary())
    
    # Extract key metrics
    print(f"\nKey Results:")
    print(f"  R-squared: {results.rsquared:.4f}")
    print(f"  Adjusted R-squared: {results.rsquared_adj:.4f}")
    print(f"  F-statistic: {results.fvalue:.2f} (p={results.f_pvalue:.4f})")
    print(f"  AIC: {results.aic:.2f}, BIC: {results.bic:.2f}")
    
    # Coefficient table
    print(f"\nCoefficients:")
    coef_df = pd.DataFrame({
        'coef': results.params,
        'std_err': results.bse,
        't': results.tvalues,
        'p_value': results.pvalues,
        'ci_lower': results.conf_int()[0],
        'ci_upper': results.conf_int()[1]
    })
    print(coef_df)
    
    # Heteroskedasticity test (Breusch-Pagan)
    bp_test = het_breuschpagan(results.resid, X)
    print(f"\nBreusch-Pagan test for heteroskedasticity:")
    print(f"  LM statistic: {bp_test[0]:.4f}")
    print(f"  p-value: {bp_test[1]:.4f}")
    
    # VIF for multicollinearity
    vif_data = pd.DataFrame()
    vif_data["Variable"] = X.columns
    vif_data["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
    print(f"\nVariance Inflation Factors (VIF):")
    print(vif_data)
    
    # Residual plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Residuals vs Fitted
    axes[0, 0].scatter(results.fittedvalues, results.resid, alpha=0.6)
    axes[0, 0].axhline(y=0, color='r', linestyle='--')
    axes[0, 0].set_xlabel('Fitted values')
    axes[0, 0].set_ylabel('Residuals')
    axes[0, 0].set_title('Residuals vs Fitted')
    
    # Q-Q plot
    from scipy import stats
    stats.probplot(results.resid, dist="norm", plot=axes[0, 1])
    axes[0, 1].set_title('Normal Q-Q Plot')
    
    # Scale-Location
    axes[1, 0].scatter(results.fittedvalues, np.sqrt(np.abs(results.resid)), alpha=0.6)
    axes[1, 0].set_xlabel('Fitted values')
    axes[1, 0].set_ylabel('√|Residuals|')
    axes[1, 0].set_title('Scale-Location')
    
    # Residuals vs Leverage
    axes[1, 1].scatter(results.get_influence().hat_matrix_diag, results.resid, alpha=0.6)
    axes[1, 1].set_xlabel('Leverage')
    axes[1, 1].set_ylabel('Residuals')
    axes[1, 1].set_title('Residuals vs Leverage')
    
    plt.tight_layout()
    plt.savefig('statsmodels_ols_diagnostics.png', dpi=150)
    print("\nSaved: statsmodels_ols_diagnostics.png")
    plt.close()
    
    return results


# ==============================================================================
# Example 2: Logistic Regression
# ==============================================================================

def example_logistic_regression():
    """Binary logistic regression with odds ratios."""
    from statsmodels.discrete.discrete_model import Logit
    from sklearn.datasets import make_classification
    from sklearn.metrics import classification_report, roc_auc_score, roc_curve
    
    print("\n" + "=" * 60)
    print("Example 2: Logistic Regression")
    print("=" * 60)
    
    # Generate binary classification data
    X, y = make_classification(
        n_samples=500, n_features=5,
        n_informative=3, n_redundant=1,
        n_classes=2, random_state=42
    )
    
    feature_names = [f'Feature_{i}' for i in range(5)]
    X = pd.DataFrame(X, columns=feature_names)
    X = sm.add_constant(X)
    
    # Fit logistic regression
    model = Logit(y, X)
    results = model.fit()
    
    print(results.summary())
    
    # Odds ratios
    odds_ratios = np.exp(results.params)
    or_ci = np.exp(results.conf_int())
    
    print(f"\nOdds Ratios:")
    or_df = pd.DataFrame({
        'OR': odds_ratios,
        'CI_lower': or_ci[0],
        'CI_upper': or_ci[1]
    })
    print(or_df)
    
    # Model evaluation
    y_pred_prob = results.predict(X)
    y_pred = (y_pred_prob > 0.5).astype(int)
    
    print(f"\nClassification Metrics:")
    print(f"  AUC: {roc_auc_score(y, y_pred_prob):.4f}")
    print(f"  Accuracy: {(y_pred == y).mean():.4f}")
    
    # Marginal effects
    margeff = results.get_margeff()
    print(f"\nMarginal Effects (at means):")
    print(margeff.summary())
    
    # ROC Curve
    fpr, tpr, thresholds = roc_curve(y, y_pred_prob)
    
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, label=f'ROC Curve (AUC = {roc_auc_score(y, y_pred_prob):.3f})')
    plt.plot([0, 1], [0, 1], 'k--', label='Random')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve - Logistic Regression')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('statsmodels_logistic_roc.png', dpi=150)
    print("Saved: statsmodels_logistic_roc.png")
    plt.close()
    
    return results


# ==============================================================================
# Example 3: Poisson Regression (Count Data)
# ==============================================================================

def example_poisson_regression():
    """Poisson regression for count data."""
    import statsmodels.api as sm
    from statsmodels.discrete.count_model import NegativeBinomial
    
    print("\n" + "=" * 60)
    print("Example 3: Poisson Regression (Count Data)")
    print("=" * 60)
    
    # Generate count data
    np.random.seed(42)
    n = 300
    
    X1 = np.random.randn(n)
    X2 = np.random.randn(n)
    
    # Poisson rate parameter
    log_lambda = 1 + 0.5 * X1 - 0.3 * X2
    lambda_ = np.exp(log_lambda)
    
    # Generate counts
    y = np.random.poisson(lambda_)
    
    df = pd.DataFrame({'X1': X1, 'X2': X2, 'y': y})
    X = sm.add_constant(df[['X1', 'X2']])
    
    # Fit Poisson model
    poisson_model = sm.GLM(df['y'], X, family=sm.families.Poisson())
    poisson_results = poisson_model.fit()
    
    print("Poisson Regression Results:")
    print(poisson_results.summary())
    
    # Rate ratios
    rate_ratios = np.exp(poisson_results.params)
    print(f"\nRate Ratios (incidence rate ratios):")
    for name, rr in rate_ratios.items():
        print(f"  {name}: {rr:.4f}")
    
    # Check overdispersion
    overdispersion = poisson_results.pearson_chi2 / poisson_results.df_resid
    print(f"\nOverdispersion ratio: {overdispersion:.4f}")
    print(f"  (> 1.5 suggests Negative Binomial may be better)")
    
    # If overdispersed, fit Negative Binomial
    if overdispersion > 1.5:
        print("\nFitting Negative Binomial (due to overdispersion)...")
        nb_model = NegativeBinomial(df['y'], X)
        nb_results = nb_model.fit()
        print(nb_results.summary())
        
        # Compare models
        print(f"\nModel Comparison:")
        print(f"  Poisson AIC: {poisson_results.aic:.2f}")
        print(f"  NegBinom AIC: {nb_results.aic:.2f}")
    
    return poisson_results


# ==============================================================================
# Example 4: Time Series Analysis (ARIMA)
# ==============================================================================

def example_arima_timeseries():
    """ARIMA modeling for time series."""
    from statsmodels.tsa.arima.model import ARIMA
    from statsmodels.tsa.stattools import adfuller
    from statsmodels.graphics.tsaplots import plot_acf, plot_pacf
    from statsmodels.stats.diagnostic import acorr_ljungbox
    
    print("\n" + "=" * 60)
    print("Example 4: ARIMA Time Series Analysis")
    print("=" * 60)
    
    # Generate synthetic time series with trend and seasonality
    np.random.seed(42)
    n = 200
    t = np.arange(n)
    
    trend = 0.02 * t
    seasonal = 10 * np.sin(2 * np.pi * t / 12)
    noise = np.random.randn(n) * 2
    
    y = 100 + trend + seasonal + noise
    
    # Create pandas Series with date index
    dates = pd.date_range(start='2020-01-01', periods=n, freq='M')
    series = pd.Series(y, index=dates)
    
    print(f"Time series length: {len(series)}")
    
    # Check stationarity (ADF test)
    adf_result = adfuller(series)
    print(f"\nADF Test for Stationarity:")
    print(f"  Test statistic: {adf_result[0]:.4f}")
    print(f"  p-value: {adf_result[1]:.4f}")
    print(f"  Critical values: {adf_result[4]}")
    
    if adf_result[1] > 0.05:
        print("  Series is non-stationary (p > 0.05)")
    
    # Plot ACF and PACF
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    plot_acf(series, lags=40, ax=axes[0])
    axes[0].set_title('Autocorrelation Function (ACF)')
    
    plot_pacf(series, lags=40, ax=axes[1])
    axes[1].set_title('Partial Autocorrelation Function (PACF)')
    
    plt.tight_layout()
    plt.savefig('statsmodels_acf_pacf.png', dpi=150)
    print("Saved: statsmodels_acf_pacf.png")
    plt.close()
    
    # Fit ARIMA model
    print("\nFitting ARIMA(1,1,1) model...")
    model = ARIMA(series, order=(1, 1, 1))
    results = model.fit()
    
    print(results.summary())
    
    # Ljung-Box test for residual autocorrelation
    lb_test = acorr_ljungbox(results.resid, lags=10)
    print(f"\nLjung-Box Test (residual autocorrelation):")
    print(lb_test)
    
    # Forecast
    forecast_steps = 12
    forecast = results.get_forecast(steps=forecast_steps)
    forecast_df = forecast.summary_frame()
    
    print(f"\n{forecast_steps}-step Forecast:")
    print(forecast_df)
    
    # Plot series with forecast
    plt.figure(figsize=(12, 6))
    plt.plot(series.index, series, label='Observed', color='blue')
    
    # Forecast dates
    forecast_dates = pd.date_range(
        start=series.index[-1] + pd.DateOffset(months=1),
        periods=forecast_steps, freq='M'
    )
    
    plt.plot(forecast_dates, forecast_df['mean'], label='Forecast', color='red')
    plt.fill_between(
        forecast_dates,
        forecast_df['mean_ci_lower'],
        forecast_df['mean_ci_upper'],
        color='red', alpha=0.2, label='95% CI'
    )
    
    plt.xlabel('Date')
    plt.ylabel('Value')
    plt.title('ARIMA Forecast')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('statsmodels_arima_forecast.png', dpi=150)
    print("Saved: statsmodels_arima_forecast.png")
    plt.close()
    
    return results


# ==============================================================================
# Example 5: Formula API (R-style)
# ==============================================================================

def example_formula_api():
    """Using statsmodels formula API for convenient model specification."""
    import statsmodels.formula.api as smf
    
    print("\n" + "=" * 60)
    print("Example 5: Formula API (R-style)")
    print("=" * 60)
    
    # Load example dataset
    df = sm.datasets.get_rdataset('mtcars').data
    print(f"Dataset: mtcars ({df.shape})")
    print(df.head())
    
    # OLS with formula
    model1 = smf.ols('mpg ~ wt + hp', data=df).fit()
    print("\nModel 1: mpg ~ wt + hp")
    print(model1.summary().tables[1])
    
    # With transformations
    model2 = smf.ols('mpg ~ np.log(wt) + hp + I(hp**2)', data=df).fit()
    print("\nModel 2: mpg ~ log(wt) + hp + hp^2")
    print(model2.summary().tables[1])
    
    # With interactions
    model3 = smf.ols('mpg ~ wt * hp', data=df).fit()  # wt + hp + wt:hp
    print("\nModel 3: mpg ~ wt * hp (main effects + interaction)")
    print(model3.summary().tables[1])
    
    # Categorical variables
    # Create categorical version
    df['cyl_cat'] = df['cyl'].astype('category')
    model4 = smf.ols('mpg ~ wt + C(cyl_cat)', data=df).fit()
    print("\nModel 4: mpg ~ wt + C(cyl_cat) (categorical)")
    print(model4.summary().tables[1])
    
    # Logistic regression with formula
    df['high_mpg'] = (df['mpg'] > df['mpg'].median()).astype(int)
    logit_model = smf.logit('high_mpg ~ wt + hp', data=df).fit()
    print("\nLogistic: high_mpg ~ wt + hp")
    print(logit_model.summary().tables[1])
    
    return model1, model2, model3, model4


# ==============================================================================
# Example 6: Robust Regression and Outlier Detection
# ==============================================================================

def example_robust_regression():
    """Robust regression methods and outlier analysis."""
    import statsmodels.api as sm
    from statsmodels.robust.robust_linear_model import RLM
    
    print("\n" + "=" * 60)
    print("Example 6: Robust Regression")
    print("=" * 60)
    
    # Generate data with outliers
    np.random.seed(42)
    n = 100
    X = np.linspace(0, 10, n)
    y_true = 2 + 3 * X
    y_clean = y_true + np.random.randn(n) * 2
    
    # Add outliers
    outlier_idx = [10, 25, 50, 75, 90]
    y_outliers = y_clean.copy()
    y_outliers[outlier_idx] += [20, -15, 25, -20, 18]
    
    X_const = sm.add_constant(X)
    
    # Standard OLS
    ols_model = sm.OLS(y_outliers, X_const)
    ols_results = ols_model.fit()
    
    # Robust regression (Huber T norm)
    huber_model = RLM(y_outliers, X_const, M=sm.robust.norms.HuberT())
    huber_results = huber_model.fit()
    
    # Robust regression (Tukey Biweight)
    tukey_model = RLM(y_outliers, X_const, M=sm.robust.norms.TukeyBiweight())
    tukey_results = tukey_model.fit()
    
    print("Comparison: OLS vs Robust Regression")
    print("-" * 50)
    print(f"{'Parameter':<15} {'OLS':<12} {'Huber':<12} {'Tukey':<12}")
    print("-" * 50)
    print(f"{'Intercept':<15} {ols_results.params[0]:<12.4f} {huber_results.params[0]:<12.4f} {tukey_results.params[0]:<12.4f}")
    print(f"{'Slope':<15} {ols_results.params[1]:<12.4f} {huber_results.params[1]:<12.4f} {tukey_results.params[1]:<12.4f}")
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.scatter(X, y_outliers, alpha=0.5, label='Data (with outliers)')
    plt.scatter(X[outlier_idx], y_outliers[outlier_idx], 
                color='red', s=100, marker='x', label='Outliers')
    
    plt.plot(X, y_true, 'g-', label='True relationship', linewidth=2)
    plt.plot(X, ols_results.fittedvalues, 'b--', label='OLS', linewidth=2)
    plt.plot(X, huber_results.fittedvalues, 'r:', label='Huber Robust', linewidth=2)
    plt.plot(X, tukey_results.fittedvalues, 'm-.', label='Tukey Robust', linewidth=2)
    
    plt.xlabel('X')
    plt.ylabel('y')
    plt.title('Robust Regression: Handling Outliers')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('statsmodels_robust_regression.png', dpi=150)
    print("\nSaved: statsmodels_robust_regression.png")
    plt.close()
    
    return ols_results, huber_results, tukey_results


# ==============================================================================
# Example 7: ANOVA and Statistical Tests
# ==============================================================================

def example_anova_tests():
    """ANOVA and statistical hypothesis testing."""
    from statsmodels.stats.anova import anova_lm
    from scipy import stats
    
    print("\n" + "=" * 60)
    print("Example 7: ANOVA and Statistical Tests")
    print("=" * 60)
    
    # Generate data for 3 groups
    np.random.seed(42)
    group_a = np.random.normal(10, 2, 30)
    group_b = np.random.normal(12, 2, 30)
    group_c = np.random.normal(11, 2, 30)
    
    # Create DataFrame
    df = pd.DataFrame({
        'value': np.concatenate([group_a, group_b, group_c]),
        'group': ['A']*30 + ['B']*30 + ['C']*30
    })
    
    print(f"Group means:")
    print(df.groupby('group')['value'].mean())
    
    # One-way ANOVA using formula
    import statsmodels.formula.api as smf
    model = smf.ols('value ~ C(group)', data=df).fit()
    anova_table = anova_lm(model)
    
    print(f"\nOne-way ANOVA:")
    print(anova_table)
    
    # F-test interpretation
    f_stat = anova_table['F'][0]
    p_value = anova_table['PR(>F)'][0]
    print(f"\nF-statistic: {f_stat:.4f}")
    print(f"p-value: {p_value:.4f}")
    
    if p_value < 0.05:
        print("Significant difference between groups (p < 0.05)")
    
    # Post-hoc: Pairwise t-tests
    print(f"\nPairwise t-tests (post-hoc):")
    groups = ['A', 'B', 'C']
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            g1 = df[df['group'] == groups[i]]['value']
            g2 = df[df['group'] == groups[j]]['value']
            t_stat, p_val = stats.ttest_ind(g1, g2)
            print(f"  {groups[i]} vs {groups[j]}: t={t_stat:.3f}, p={p_val:.4f}")
    
    # Visualize
    plt.figure(figsize=(10, 6))
    df.boxplot(column='value', by='group', ax=plt.gca())
    plt.title('One-way ANOVA: Group Comparisons')
    plt.suptitle('')  # Remove default title
    plt.savefig('statsmodels_anova.png', dpi=150)
    print("Saved: statsmodels_anova.png")
    plt.close()
    
    return anova_table


# ==============================================================================
# Main Execution
# ==============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("Statsmodels: Statistical Modeling Examples")
    print("=" * 70)
    
    try:
        example_ols_regression()
        example_logistic_regression()
        example_poisson_regression()
        example_arima_timeseries()
        example_formula_api()
        example_robust_regression()
        example_anova_tests()
        
        print("\n" + "=" * 70)
        print("All examples completed successfully!")
        print("Generated plots:")
        print("  - statsmodels_ols_diagnostics.png")
        print("  - statsmodels_logistic_roc.png")
        print("  - statsmodels_acf_pacf.png")
        print("  - statsmodels_arima_forecast.png")
        print("  - statsmodels_robust_regression.png")
        print("  - statsmodels_anova.png")
        print("=" * 70)
        
    except ImportError as e:
        print(f"\nImport Error: {e}")
        print("Please install required packages:")
        print("  pip install statsmodels scipy matplotlib pandas")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
