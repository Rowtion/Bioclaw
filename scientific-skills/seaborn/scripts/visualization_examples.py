#!/usr/bin/env python3
"""
Seaborn: Statistical Visualization Examples
============================================
Publication-quality statistical graphics with pandas integration.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Set style for publication-quality plots
sns.set_theme(style='whitegrid', context='paper', font_scale=1.1)

# ==============================================================================
# Example 1: Distribution Plots
# ==============================================================================

def example_distributions():
    """Various distribution visualization techniques."""
    print("=" * 60)
    print("Example 1: Distribution Plots")
    print("=" * 60)
    
    # Load example dataset
    df = sns.load_dataset('tips')
    print(f"Dataset shape: {df.shape}")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # 1. Histogram
    sns.histplot(data=df, x='total_bill', bins=20, kde=True, ax=axes[0, 0])
    axes[0, 0].set_title('Histogram with KDE')
    
    # 2. Box plot
    sns.boxplot(data=df, x='day', y='total_bill', ax=axes[0, 1])
    axes[0, 1].set_title('Box Plot by Day')
    
    # 3. Violin plot
    sns.violinplot(data=df, x='day', y='total_bill', split=True, hue='sex', ax=axes[0, 2])
    axes[0, 2].set_title('Violin Plot (Split by Sex)')
    
    # 4. ECDF plot
    sns.ecdfplot(data=df, x='total_bill', hue='time', ax=axes[1, 0])
    axes[1, 0].set_title('ECDF by Time')
    
    # 5. Strip plot with jitter
    sns.stripplot(data=df, x='day', y='tip', hue='smoker', dodge=True, ax=axes[1, 1])
    axes[1, 1].set_title('Strip Plot (Jittered)')
    
    # 6. Swarm plot
    sns.swarmplot(data=df, x='day', y='tip', size=4, ax=axes[1, 2])
    axes[1, 2].set_title('Swarm Plot')
    
    plt.tight_layout()
    plt.savefig('seaborn_distributions.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_distributions.png")
    plt.close()


# ==============================================================================
# Example 2: Relational Plots (Scatter and Line)
# ==============================================================================

def example_relational():
    """Scatter plots, line plots, and relational visualizations."""
    print("\n" + "=" * 60)
    print("Example 2: Relational Plots")
    print("=" * 60)
    
    # Load datasets
    tips = sns.load_dataset('tips')
    fmri = sns.load_dataset('fmri')
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Scatter with multiple encodings
    sns.scatterplot(
        data=tips, x='total_bill', y='tip',
        hue='time', size='size', style='sex',
        sizes=(50, 200), ax=axes[0, 0]
    )
    axes[0, 0].set_title('Scatter: Multiple Encodings')
    
    # 2. Line plot with confidence intervals
    sns.lineplot(
        data=fmri, x='timepoint', y='signal',
        hue='region', style='event', ci=95,
        ax=axes[0, 1]
    )
    axes[0, 1].set_title('Line Plot with 95% CI')
    
    # 3. Regression plot
    sns.regplot(
        data=tips, x='total_bill', y='tip',
        scatter_kws={'alpha': 0.5}, line_kws={'color': 'red'},
        ax=axes[1, 0]
    )
    axes[1, 0].set_title('Linear Regression with Confidence Interval')
    
    # 4. Residual plot
    sns.residplot(
        data=tips, x='total_bill', y='tip',
        scatter_kws={'alpha': 0.5}, ax=axes[1, 1]
    )
    axes[1, 1].set_title('Residual Plot')
    axes[1, 1].axhline(y=0, color='r', linestyle='--')
    
    plt.tight_layout()
    plt.savefig('seaborn_relational.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_relational.png")
    plt.close()


# ==============================================================================
# Example 3: Categorical Plots
# ==============================================================================

def example_categorical():
    """Categorical data visualization techniques."""
    print("\n" + "=" * 60)
    print("Example 3: Categorical Plots")
    print("=" * 60)
    
    df = sns.load_dataset('titanic')
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Bar plot with confidence intervals
    sns.barplot(
        data=df, x='class', y='survived',
        hue='sex', ci=95, ax=axes[0, 0]
    )
    axes[0, 0].set_title('Survival Rate by Class and Sex')
    axes[0, 0].set_ylabel('Survival Probability')
    
    # 2. Point plot
    sns.pointplot(
        data=df, x='class', y='survived',
        hue='sex', dodge=True, capsize=0.1, ax=axes[0, 1]
    )
    axes[0, 1].set_title('Point Plot: Survival Trends')
    
    # 3. Count plot
    sns.countplot(
        data=df, x='class', hue='survived',
        palette='Set2', ax=axes[1, 0]
    )
    axes[1, 0].set_title('Count by Class and Survival')
    
    # 4. Boxen plot (enhanced box plot for large datasets)
    sns.boxenplot(
        data=sns.load_dataset('diamonds'),
        x='cut', y='price', ax=axes[1, 1]
    )
    axes[1, 1].set_title('Boxen Plot: Diamond Prices')
    
    plt.tight_layout()
    plt.savefig('seaborn_categorical.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_categorical.png")
    plt.close()


# ==============================================================================
# Example 4: Matrix and Heatmap Plots
# ==============================================================================

def example_matrices():
    """Heatmaps and matrix visualizations."""
    print("\n" + "=" * 60)
    print("Example 4: Matrix and Heatmap Plots")
    print("=" * 60)
    
    # Load flights dataset
    flights = sns.load_dataset('flights')
    flights_pivot = flights.pivot(index='month', columns='year', values='passengers')
    
    # Load correlation dataset
    iris = sns.load_dataset('iris')
    corr = iris.select_dtypes(include=[np.number]).corr()
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Correlation heatmap
    sns.heatmap(
        corr, annot=True, fmt='.2f',
        cmap='coolwarm', center=0,
        square=True, linewidths=0.5, ax=axes[0, 0]
    )
    axes[0, 0].set_title('Correlation Matrix')
    
    # 2. Clustermap (separate figure)
    plt.figure(figsize=(8, 8))
    sns.clustermap(
        flights_pivot, cmap='YlOrRd',
        figsize=(8, 8), annot=False
    )
    plt.savefig('seaborn_clustermap.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: seaborn_clustermap.png")
    
    # 3. Flight passengers heatmap
    sns.heatmap(
        flights_pivot, cmap='YlOrRd',
        annot=True, fmt='d', linewidths=0.5, ax=axes[0, 1]
    )
    axes[0, 1].set_title('Flight Passengers by Month and Year')
    
    # 4. Masked heatmap (upper triangle)
    mask = np.triu(np.ones_like(corr, dtype=bool))
    sns.heatmap(
        corr, mask=mask, annot=True, fmt='.2f',
        cmap='coolwarm', center=0, ax=axes[1, 0]
    )
    axes[1, 0].set_title('Lower Triangle Correlation')
    
    # 5. Diverging heatmap with custom center
    np.random.seed(42)
    data = pd.DataFrame(
        np.random.randn(10, 10),
        columns=[f'Var{i}' for i in range(10)]
    )
    sns.heatmap(
        data, cmap='vlag', center=0,
        annot=True, fmt='.1f', ax=axes[1, 1]
    )
    axes[1, 1].set_title('Diverging Heatmap')
    
    plt.tight_layout()
    plt.savefig('seaborn_matrices.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_matrices.png")
    plt.close()


# ==============================================================================
# Example 5: Multi-Plot Grids (FacetGrid, PairGrid, JointGrid)
# ==============================================================================

def example_grids():
    """Multi-panel figure creation with grids."""
    print("\n" + "=" * 60)
    print("Example 5: Multi-Plot Grids")
    print("=" * 60)
    
    tips = sns.load_dataset('tips')
    iris = sns.load_dataset('iris')
    
    # 1. FacetGrid
    print("Creating FacetGrid...")
    g = sns.FacetGrid(tips, col='time', row='smoker', height=3, aspect=1.2)
    g.map(sns.scatterplot, 'total_bill', 'tip', hue='sex', data=tips)
    g.add_legend()
    g.fig.suptitle('FacetGrid: Tips by Time and Smoker Status', y=1.02)
    plt.savefig('seaborn_facetgrid.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_facetgrid.png")
    plt.close()
    
    # 2. PairGrid
    print("Creating PairGrid...")
    g = sns.PairGrid(iris, hue='species', height=2.5)
    g.map_upper(sns.scatterplot, alpha=0.6)
    g.map_lower(sns.kdeplot, fill=True)
    g.map_diag(sns.histplot, kde=True)
    g.add_legend()
    g.fig.suptitle('PairGrid: Iris Dataset', y=1.02)
    plt.savefig('seaborn_pairgrid.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_pairgrid.png")
    plt.close()
    
    # 3. JointGrid
    print("Creating JointGrid...")
    g = sns.JointGrid(data=tips, x='total_bill', y='tip', height=7)
    g.plot_joint(sns.scatterplot, hue=tips['time'], alpha=0.6)
    g.plot_marginals(sns.histplot, kde=True, color='steelblue')
    g.fig.suptitle('JointGrid: Bill vs Tip Distribution', y=1.02)
    plt.savefig('seaborn_jointgrid.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_jointgrid.png")
    plt.close()


# ==============================================================================
# Example 6: Figure-Level Functions (relplot, displot, catplot, lmplot)
# ==============================================================================

def example_figure_level():
    """High-level figure functions with automatic faceting."""
    print("\n" + "=" * 60)
    print("Example 6: Figure-Level Functions")
    print("=" * 60)
    
    tips = sns.load_dataset('tips')
    
    # 1. relplot (relational plots with faceting)
    print("Creating relplot...")
    g = sns.relplot(
        data=tips, x='total_bill', y='tip',
        col='time', row='smoker', hue='sex',
        kind='scatter', height=3, aspect=1.2
    )
    g.fig.suptitle('relplot: Faceted Scatter Plot', y=1.02)
    plt.savefig('seaborn_relplot.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_relplot.png")
    plt.close()
    
    # 2. displot (distribution plots)
    print("Creating displot...")
    g = sns.displot(
        data=tips, x='total_bill', col='time',
        kind='kde', fill=True, height=4, aspect=1.2
    )
    g.fig.suptitle('displot: KDE by Time', y=1.02)
    plt.savefig('seaborn_displot.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_displot.png")
    plt.close()
    
    # 3. catplot (categorical plots)
    print("Creating catplot...")
    g = sns.catplot(
        data=tips, x='day', y='total_bill',
        col='smoker', hue='sex',
        kind='box', height=4, aspect=1
    )
    g.fig.suptitle('catplot: Box Plots by Day', y=1.02)
    plt.savefig('seaborn_catplot.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_catplot.png")
    plt.close()
    
    # 4. lmplot (regression plots)
    print("Creating lmplot...")
    g = sns.lmplot(
        data=tips, x='total_bill', y='tip',
        col='time', hue='smoker',
        height=4, aspect=1.2
    )
    g.fig.suptitle('lmplot: Regression by Time and Smoking', y=1.02)
    plt.savefig('seaborn_lmplot.png', dpi=150, bbox_inches='tight')
    print("Saved: seaborn_lmplot.png")
    plt.close()


# ==============================================================================
# Example 7: Objects Interface (Modern API)
# ==============================================================================

def example_objects_interface():
    """Modern seaborn.objects interface for declarative plotting."""
    print("\n" + "=" * 60)
    print("Example 7: Objects Interface (Modern API)")
    print("=" * 60)
    
    try:
        from seaborn import objects as so
        
        tips = sns.load_dataset('tips')
        
        # Create figure with objects interface
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        
        # 1. Basic scatter with objects
        p = (
            so.Plot(tips, x='total_bill', y='tip', color='time')
            .add(so.Dot())
            .on(axes[0, 0])
        )
        p.show()
        axes[0, 0].set_title('Objects: Scatter with Color')
        
        # 2. With regression line
        p = (
            so.Plot(tips, x='total_bill', y='tip')
            .add(so.Dot(alpha=0.5))
            .add(so.Line(), so.PolyFit(order=1))
            .on(axes[0, 1])
        )
        p.show()
        axes[0, 1].set_title('Objects: Scatter with Fit')
        
        # 3. Dodged bars
        p = (
            so.Plot(tips, x='day', y='tip', color='sex')
            .add(so.Bar(), so.Agg('mean'), so.Dodge())
            .on(axes[1, 0])
        )
        p.show()
        axes[1, 0].set_title('Objects: Dodged Bar Means')
        
        # 4. Faceted
        p = (
            so.Plot(tips, x='total_bill', y='tip')
            .add(so.Dot())
            .facet('time', 'smoker')
            .on(axes[1, 1])
        )
        p.show()
        axes[1, 1].set_title('Objects: Faceted')
        
        plt.tight_layout()
        plt.savefig('seaborn_objects.png', dpi=150, bbox_inches='tight')
        print("Saved: seaborn_objects.png")
        plt.close()
        
    except ImportError:
        print("seaborn.objects not available in this version")


# ==============================================================================
# Example 8: Publication-Quality Figure
# ==============================================================================

def example_publication_figure():
    """Create a publication-ready multi-panel figure."""
    print("\n" + "=" * 60)
    print("Example 8: Publication-Quality Figure")
    print("=" * 60)
    
    # Reset to paper style
    sns.set_theme(style='ticks', context='paper', font_scale=1.2)
    
    tips = sns.load_dataset('tips')
    
    fig = plt.figure(figsize=(14, 10))
    
    # Create grid spec for custom layout
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
    
    # Panel A: Main scatter
    ax1 = fig.add_subplot(gs[0, :2])
    sns.scatterplot(
        data=tips, x='total_bill', y='tip',
        hue='time', size='size', sizes=(50, 300),
        alpha=0.7, ax=ax1
    )
    ax1.set_title('A. Tip Amount vs Total Bill', fontweight='bold', loc='left')
    sns.despine(ax=ax1)
    
    # Panel B: Distribution by time
    ax2 = fig.add_subplot(gs[0, 2])
    sns.boxplot(
        data=tips, x='time', y='tip',
        palette='Set2', ax=ax2
    )
    ax2.set_title('B. Tip Distribution', fontweight='bold', loc='left')
    sns.despine(ax=ax2, trim=True)
    
    # Panel C: Day comparison
    ax3 = fig.add_subplot(gs[1, 0])
    sns.barplot(
        data=tips, x='day', y='tip',
        ci=95, palette='viridis', ax=ax3
    )
    ax3.set_title('C. Mean Tip by Day', fontweight='bold', loc='left')
    sns.despine(ax=ax3)
    
    # Panel D: Smoker analysis
    ax4 = fig.add_subplot(gs[1, 1])
    sns.violinplot(
        data=tips, x='smoker', y='total_bill',
        split=True, hue='sex', ax=ax4
    )
    ax4.set_title('D. Bill Distribution', fontweight='bold', loc='left')
    sns.despine(ax=ax4)
    
    # Panel E: Size relationship
    ax5 = fig.add_subplot(gs[1, 2])
    sns.pointplot(
        data=tips, x='size', y='tip',
        ci=95, color='steelblue', ax=ax5
    )
    ax5.set_title('E. Tip by Party Size', fontweight='bold', loc='left')
    sns.despine(ax=ax5)
    
    plt.suptitle('Restaurant Tipping Analysis', fontsize=16, fontweight='bold', y=0.98)
    plt.savefig('seaborn_publication.png', dpi=300, bbox_inches='tight')
    print("Saved: seaborn_publication.png (300 DPI)")
    plt.close()


# ==============================================================================
# Main Execution
# ==============================================================================

if __name__ == "__main__":
    print("\n" + "=" * 70)
    print("Seaborn: Statistical Visualization Examples")
    print("=" * 70)
    
    try:
        example_distributions()
        example_relational()
        example_categorical()
        example_matrices()
        example_grids()
        example_figure_level()
        example_objects_interface()
        example_publication_figure()
        
        print("\n" + "=" * 70)
        print("All examples completed successfully!")
        print("Generated plots:")
        print("  - seaborn_distributions.png")
        print("  - seaborn_relational.png")
        print("  - seaborn_categorical.png")
        print("  - seaborn_matrices.png")
        print("  - seaborn_clustermap.png")
        print("  - seaborn_facetgrid.png")
        print("  - seaborn_pairgrid.png")
        print("  - seaborn_jointgrid.png")
        print("  - seaborn_relplot.png")
        print("  - seaborn_displot.png")
        print("  - seaborn_catplot.png")
        print("  - seaborn_lmplot.png")
        print("  - seaborn_objects.png")
        print("  - seaborn_publication.png")
        print("=" * 70)
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
