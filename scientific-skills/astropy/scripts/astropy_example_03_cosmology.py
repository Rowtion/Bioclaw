#!/usr/bin/env python3
"""
astropy_example_03_cosmology.py
宇宙学计算

本脚本演示:
1. 使用内置宇宙学模型
2. 距离和红移计算
3. 宇宙年龄和回望时间
4. 体积计算
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import (
    Planck18, Planck15, WMAP9, 
    FlatLambdaCDM, LambdaCDM
)
from astropy import units as u


def builtin_cosmologies():
    """使用内置宇宙学模型"""
    print("="*60)
    print("内置宇宙学模型")
    print("="*60)
    
    # 常用的内置宇宙学
    cosmologies = {
        'Planck 2018': Planck18,
        'Planck 2015': Planck15,
        'WMAP 9-year': WMAP9,
    }
    
    print("\n基本参数:")
    print(f"{'模型':<15} {'H0':<10} {'Om0':<10} {'Ode0':<10} {'Age (Gyr)':<12}")
    print("-" * 60)
    
    for name, cosmo in cosmologies.items():
        print(f"{name:<15} {cosmo.H0.value:<10.2f} {cosmo.Om0:<10.3f} "
              f"{cosmo.Ode0:<10.3f} {cosmo.age(0).to(u.Gyr).value:<12.2f}")
    
    return Planck18  # 使用 Planck 2018 作为默认


def custom_cosmology():
    """创建自定义宇宙学模型"""
    print("\n" + "="*60)
    print("自定义宇宙学模型")
    print("="*60)
    
    # 平坦宇宙 (标准 LCDM)
    flat_lcdm = FlatLambdaCDM(
        H0=70 * u.km / u.s / u.Mpc,  # 哈勃常数
        Om0=0.3,  # 物质密度参数
        name='My Flat LCDM'
    )
    print(f"\n平坦 ΛCDM:")
    print(f"  H0 = {flat_lcdm.H0}")
    print(f"  Ωm = {flat_lcdm.Om0}")
    print(f"  ΩΛ = {flat_lcdm.Ode0}")
    print(f"  曲率 = {flat_lcdm.Ok0}")
    
    # 非平坦宇宙
    curved_lcdm = LambdaCDM(
        H0=70 * u.km / u.s / u.Mpc,
        Om0=0.3,
        Ode0=0.6,  # 非平坦: Ωm + ΩΛ ≠ 1
        name='Curved LCDM'
    )
    print(f"\n非平坦 ΛCDM:")
    print(f"  H0 = {curved_lcdm.H0}")
    print(f"  Ωm = {curved_lcdm.Om0}")
    print(f"  ΩΛ = {curved_lcdm.Ode0}")
    print(f"  曲率 Ωk = {curved_lcdm.Ok0:.3f}")
    
    return flat_lcdm


def distance_calculations(cosmo):
    """距离计算"""
    print("\n" + "="*60)
    print("宇宙学距离计算")
    print("="*60)
    
    # 测试红移
    z = 1.0
    
    print(f"\n在 z = {z} 处的各种距离:")
    print("-" * 50)
    
    # 共动距离
    d_comoving = cosmo.comoving_distance(z)
    print(f"共动距离:           {d_comoving:.2f}")
    
    # 光度距离
    d_luminosity = cosmo.luminosity_distance(z)
    print(f"光度距离:           {d_luminosity:.2f}")
    
    # 角直径距离
    d_angular = cosmo.angular_diameter_distance(z)
    print(f"角直径距离:         {d_angular:.2f}")
    
    # 视线共动距离
    d_line_of_sight = cosmo.comoving_distance(z)  # 对于平坦宇宙
    print(f"视线共动距离:       {d_line_of_sight:.2f}")
    
    # 横向共动距离
    d_transverse = cosmo.comoving_transverse_distance(z)
    print(f"横向共动距离:       {d_transverse:.2f}")
    
    # 距离模数
    distmod = cosmo.distmod(z)
    print(f"距离模数:           {distmod:.2f}")
    
    # 验证关系: d_L = d_A * (1+z)^2
    print(f"\n验证距离关系:")
    print(f"  d_L = d_A * (1+z)^2")
    print(f"  {d_luminosity:.2f} = {d_angular * (1+z)**2:.2f}")
    
    return z, d_luminosity, d_angular


def age_and_time(cosmo):
    """宇宙年龄和时间计算"""
    print("\n" + "="*60)
    print("宇宙年龄和回望时间")
    print("="*60)
    
    # 当前宇宙年龄
    age_now = cosmo.age(0)
    print(f"\n当前宇宙年龄: {age_now.to(u.Gyr):.2f}")
    
    # 不同红移处的宇宙年龄和回望时间
    redshifts = [0, 0.5, 1.0, 2.0, 5.0, 10.0]
    
    print(f"\n{'红移':<8} {'宇宙年龄 (Gyr)':<18} {'回望时间 (Gyr)':<18}")
    print("-" * 50)
    
    for z in redshifts:
        age_at_z = cosmo.age(z)
        lookback = cosmo.lookback_time(z)
        print(f"{z:<8.1f} {age_at_z.to(u.Gyr).value:<18.2f} {lookback.to(u.Gyr).value:<18.2f}")


def scale_factor(cosmo):
    """尺度因子演化"""
    print("\n" + "="*60)
    print("尺度因子演化")
    print("="*60)
    
    print(f"\n尺度因子 a = 1/(1+z):")
    print(f"{'红移 z':<10} {'尺度因子 a':<15} {'说明':<30}")
    print("-" * 60)
    
    scenarios = [
        (0, "现在"),
        (0.5, ""),
        (1.0, "宇宙年龄为现在的一半"),
        (2.0, ""),
        (9.0, "宇宙年龄约为现在的一成"),
        (99.0, "复合时期 (CMB)")
    ]
    
    for z, desc in scenarios:
        a = cosmo.scale_factor(z)
        print(f"{z:<10.1f} {a:<15.4f} {desc:<30}")
    
    # 临界密度演化
    print(f"\n临界密度演化:")
    for z in [0, 1, 2, 5]:
        rho_crit = cosmo.critical_density(z)
        print(f"  z={z}: ρ_crit = {rho_crit:.3e}")


def hubble_parameter(cosmo):
    """哈勃参数演化"""
    print("\n" + "="*60)
    print("哈勃参数演化")
    print("="*60)
    
    print(f"\n哈勃参数 H(z):")
    print(f"{'红移 z':<10} {'H(z) km/s/Mpc':<18} {'H(z)/H0':<12}")
    print("-" * 45)
    
    for z in [0, 0.5, 1.0, 2.0, 5.0, 10.0]:
        H = cosmo.H(z)
        H_ratio = H / cosmo.H0
        print(f"{z:<10.1f} {H.value:<18.2f} {H_ratio:<12.2f}")


def volume_calculations(cosmo):
    """体积计算"""
    print("\n" + "="*60)
    print("体积计算")
    print("="*60)
    
    # 共动体积
    z = 1.0
    V_comoving = cosmo.comoving_volume(z)
    print(f"\n到 z={z} 的共动体积:")
    print(f"  V = {V_comoving:.3e}")
    print(f"    = {V_comoving.to(u.Gpc**3):.2f}")
    
    # 微分共动体积 (单位红移、单位立体角)
    dV_dz = cosmo.differential_comoving_volume(z)
    print(f"\n微分共动体积 (z={z}):")
    print(f"  dV/dz/dΩ = {dV_dz:.3e}")
    
    # 球壳体积
    z1, z2 = 0.5, 1.0
    V_shell = cosmo.comoving_volume(z2) - cosmo.comoving_volume(z1)
    print(f"\n红移范围 {z1} < z < {z2} 的球壳体积:")
    print(f"  V_shell = {V_shell.to(u.Gpc**3):.2f}")
    
    # 可观测宇宙体积
    V_observable = cosmo.comoving_volume(np.inf)
    print(f"\n可观测宇宙共动体积:")
    print(f"  V_observable = {V_observable.to(u.Gpc**3):.2f}")


def inverse_calculations(cosmo):
    """反向计算 (从距离求红移)"""
    print("\n" + "="*60)
    print("反向计算")
    print("="*60)
    
    # 从光度距离求红移
    target_dL = 1000 * u.Mpc
    z = cosmo.z_at_value(cosmo.luminosity_distance, target_dL)
    print(f"\n光度距离 dL = {target_dL} 对应的红移:")
    print(f"  z = {z:.4f}")
    
    # 验证
    dL_check = cosmo.luminosity_distance(z)
    print(f"  验证: dL(z) = {dL_check:.2f}")
    
    # 从年龄求红移
    target_age = 5 * u.Gyr
    z_age = cosmo.z_at_value(cosmo.age, target_age)
    print(f"\n宇宙年龄 = {target_age} 对应的红移:")
    print(f"  z = {z_age:.4f}")
    
    age_check = cosmo.age(z_age)
    print(f"  验证: age(z) = {age_check:.2f}")


def visualize_cosmology(cosmo):
    """可视化宇宙学量"""
    print("\n" + "="*60)
    print("宇宙学可视化")
    print("="*60)
    
    # 红移范围
    z = np.linspace(0, 5, 200)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. 距离-红移关系
    ax1 = axes[0, 0]
    d_comoving = [cosmo.comoving_distance(zi).value for zi in z]
    d_luminosity = [cosmo.luminosity_distance(zi).value for zi in z]
    d_angular = [cosmo.angular_diameter_distance(zi).value for zi in z]
    
    ax1.plot(z, d_comoving, label='Comoving', linewidth=2)
    ax1.plot(z, d_luminosity, label='Luminosity', linewidth=2)
    ax1.plot(z, d_angular, label='Angular Diameter', linewidth=2)
    ax1.set_xlabel('Redshift z')
    ax1.set_ylabel('Distance (Mpc)')
    ax1.set_title('Distance-Redshift Relations')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. 宇宙年龄和回望时间
    ax2 = axes[0, 1]
    age = [cosmo.age(zi).to(u.Gyr).value for zi in z]
    lookback = [cosmo.lookback_time(zi).to(u.Gyr).value for zi in z]
    
    ax2.plot(z, age, label='Age of Universe', linewidth=2)
    ax2.plot(z, lookback, label='Lookback Time', linewidth=2)
    ax2.set_xlabel('Redshift z')
    ax2.set_ylabel('Time (Gyr)')
    ax2.set_title('Cosmic Time Evolution')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. 哈勃参数
    ax3 = axes[1, 0]
    H = [cosmo.H(zi).value for zi in z]
    ax3.plot(z, H, linewidth=2, color='purple')
    ax3.axhline(y=cosmo.H0.value, color='red', linestyle='--', 
                label=f'H0 = {cosmo.H0.value:.1f}')
    ax3.set_xlabel('Redshift z')
    ax3.set_ylabel('H(z) (km/s/Mpc)')
    ax3.set_title('Hubble Parameter Evolution')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # 4. 密度参数
    ax4 = axes[1, 1]
    Om = [cosmo.Om(zi) for zi in z]
    Ode = [cosmo.Ode(zi) for zi in z]
    
    ax4.plot(z, Om, label='Matter (Ωm)', linewidth=2)
    ax4.plot(z, Ode, label='Dark Energy (ΩΛ)', linewidth=2)
    ax4.set_xlabel('Redshift z')
    ax4.set_ylabel('Density Parameter')
    ax4.set_title('Density Evolution')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(0, 1.2)
    
    plt.tight_layout()
    plt.savefig('cosmology_plots.png', dpi=150, bbox_inches='tight')
    print("✓ 保存: cosmology_plots.png")
    plt.close()


def compare_cosmologies():
    """比较不同宇宙学模型"""
    print("\n" + "="*60)
    print("宇宙学模型比较")
    print("="*60)
    
    # 不同的哈勃常数
    cosmo1 = FlatLambdaCDM(H0=67, Om0=0.31, name='Planck-like')
    cosmo2 = FlatLambdaCDM(H0=73, Om0=0.27, name='SH0ES-like')
    
    z = np.linspace(0.01, 2, 100)
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # 距离比较
    ax1 = axes[0]
    for cosmo, label, color in [(cosmo1, 'H0=67', 'blue'), 
                                  (cosmo2, 'H0=73', 'red')]:
        dL = [cosmo.luminosity_distance(zi).value for zi in z]
        ax1.plot(z, dL, label=f'{label}, Ωm={cosmo.Om0}', 
                linewidth=2, color=color)
    ax1.set_xlabel('Redshift z')
    ax1.set_ylabel('Luminosity Distance (Mpc)')
    ax1.set_title('Distance Comparison')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 年龄比较
    ax2 = axes[1]
    for cosmo, label, color in [(cosmo1, 'H0=67', 'blue'), 
                                  (cosmo2, 'H0=73', 'red')]:
        age = [cosmo.age(zi).to(u.Gyr).value for zi in z]
        ax2.plot(z, age, label=f'{label}, Ωm={cosmo.Om0}', 
                linewidth=2, color=color)
    ax2.set_xlabel('Redshift z')
    ax2.set_ylabel('Age of Universe (Gyr)')
    ax2.set_title('Age Comparison')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('cosmology_comparison.png', dpi=150, bbox_inches='tight')
    print("✓ 保存: cosmology_comparison.png")
    plt.close()


def main():
    """主函数"""
    print("Astropy 宇宙学计算示例")
    
    # 1. 内置宇宙学
    cosmo = builtin_cosmologies()
    
    # 2. 自定义宇宙学
    custom_cosmo = custom_cosmology()
    
    # 3. 距离计算
    distance_calculations(cosmo)
    
    # 4. 年龄和时间
    age_and_time(cosmo)
    
    # 5. 尺度因子
    scale_factor(cosmo)
    
    # 6. 哈勃参数
    hubble_parameter(cosmo)
    
    # 7. 体积计算
    volume_calculations(cosmo)
    
    # 8. 反向计算
    inverse_calculations(cosmo)
    
    # 9. 可视化
    visualize_cosmology(cosmo)
    
    # 10. 模型比较
    compare_cosmologies()
    
    print("\n" + "="*60)
    print("示例完成!")
    print("="*60)


if __name__ == '__main__':
    main()
