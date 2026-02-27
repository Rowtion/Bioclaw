#!/usr/bin/env python3
"""
astropy_example_01_coordinates.py
天文坐标转换与计算

本脚本演示:
1. 创建天球坐标对象
2. 坐标系转换 (ICRS, Galactic, AltAz)
3. 计算角距离和位置角
4. 坐标匹配
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import (
    SkyCoord, EarthLocation, AltAz, 
    match_coordinates_sky, angular_separation
)
from astropy.time import Time


def create_coordinates():
    """创建各种坐标对象"""
    print("="*60)
    print("创建天文坐标")
    print("="*60)
    
    # 方法1: 使用角度单位
    coord1 = SkyCoord(ra=10.5*u.degree, dec=41.2*u.degree, frame='icrs')
    print(f"\n方法1 - 使用单位:")
    print(f"  坐标: RA={coord1.ra:.4f}, Dec={coord1.dec:.4f}")
    
    # 方法2: 使用时分秒格式
    coord2 = SkyCoord('05h23m34.5s', '-69d45m22s', frame='icrs')
    print(f"\n方法2 - 时分秒格式:")
    print(f"  坐标: RA={coord2.ra:.4f}, Dec={coord2.dec:.4f}")
    
    # 方法3: 批量创建坐标
    ras = [10.0, 15.5, 20.3, 25.7] * u.degree
    decs = [40.0, 42.5, -30.2, 15.8] * u.degree
    coords_batch = SkyCoord(ra=ras, dec=decs, frame='icrs')
    print(f"\n方法3 - 批量创建:")
    print(f"  创建了 {len(coords_batch)} 个坐标")
    
    return coord1, coord2, coords_batch


def coordinate_transformations():
    """坐标系转换示例"""
    print("\n" + "="*60)
    print("坐标系转换")
    print("="*60)
    
    # 创建一个 ICRS 坐标
    coord_icrs = SkyCoord(ra=266.4*u.degree, dec=-29.0*u.degree, frame='icrs')
    print(f"\n原始 ICRS 坐标:")
    print(f"  RA={coord_icrs.ra.to(u.hourangle):.4f}")
    print(f"  Dec={coord_icrs.dec:.4f}")
    
    # 转换为银河系坐标
    coord_gal = coord_icrs.galactic
    print(f"\n转换为银河系坐标:")
    print(f"  l={coord_gal.l:.4f}")
    print(f"  b={coord_gal.b:.4f}")
    
    # 转换为 FK5 坐标
    coord_fk5 = coord_icrs.fk5
    print(f"\n转换为 FK5 坐标:")
    print(f"  RA={coord_fk5.ra:.4f}")
    print(f"  Dec={coord_fk5.dec:.4f}")
    
    # 转换为地平坐标 (AltAz)
    # 需要观测位置和时间
    observing_location = EarthLocation(
        lat=31.0*u.deg, 
        lon=-110.0*u.deg, 
        height=2100*u.m
    )
    observing_time = Time('2024-06-15 23:00:00')
    
    aa_frame = AltAz(obstime=observing_time, location=observing_location)
    coord_altaz = coord_icrs.transform_to(aa_frame)
    
    print(f"\n转换为地平坐标 (AltAz):")
    print(f"  观测地点: 纬度={observing_location.lat}, 经度={observing_location.lon}")
    print(f"  观测时间: {observing_time}")
    print(f"  高度角 (Alt)={coord_altaz.alt:.4f}")
    print(f"  方位角 (Az)={coord_altaz.az:.4f}")
    
    # 检查目标是否在地平线上
    if coord_altaz.alt > 0*u.deg:
        print(f"  ✓ 目标在地平线上 (可观测)")
    else:
        print(f"  ✗ 目标在地平线下 (不可观测)")
    
    return coord_icrs


def calculate_distances():
    """计算角距离和位置角"""
    print("\n" + "="*60)
    print("角距离和位置角计算")
    print("="*60)
    
    # 定义两个坐标
    coord1 = SkyCoord(ra=10.0*u.degree, dec=41.0*u.degree, frame='icrs')
    coord2 = SkyCoord(ra=10.5*u.degree, dec=41.5*u.degree, frame='icrs')
    
    print(f"\n坐标1: RA={coord1.ra:.4f}, Dec={coord1.dec:.4f}")
    print(f"坐标2: RA={coord2.ra:.4f}, Dec={coord2.dec:.4f}")
    
    # 计算角距离
    separation = coord1.separation(coord2)
    print(f"\n角距离: {separation:.4f}")
    print(f"        {separation.to(u.arcmin):.2f}")
    print(f"        {separation.to(u.arcsec):.2f}")
    
    # 计算位置角 (从北向东测量)
    position_angle = coord1.position_angle(coord2)
    print(f"\n位置角: {position_angle.to(u.degree):.2f}")
    
    # 使用 angular_separation 函数
    sep = angular_separation(
        coord1.ra, coord1.dec,
        coord2.ra, coord2.dec
    )
    print(f"\n使用 angular_separation: {sep.to(u.arcmin):.2f}")
    
    return coord1, coord2


def catalog_matching():
    """星表匹配示例"""
    print("\n" + "="*60)
    print("星表匹配")
    print("="*60)
    
    # 创建模拟的星表1 (观测数据)
    np.random.seed(42)
    n_sources1 = 100
    ra1 = np.random.uniform(10, 11, n_sources1) * u.degree
    dec1 = np.random.uniform(40, 41, n_sources1) * u.degree
    catalog1 = SkyCoord(ra=ra1, dec=dec1, frame='icrs')
    
    # 创建模拟的星表2 (有轻微偏移，模拟不同观测)
    n_sources2 = 80
    ra2 = np.random.uniform(10, 11, n_sources2) * u.degree + 0.001*u.degree  # 小偏移
    dec2 = np.random.uniform(40, 41, n_sources2) * u.degree + 0.001*u.degree
    catalog2 = SkyCoord(ra=ra2, dec=dec2, frame='icrs')
    
    print(f"\n星表1: {len(catalog1)} 个源")
    print(f"星表2: {len(catalog2)} 个源")
    
    # 匹配星表
    idx, sep2d, dist3d = match_coordinates_sky(catalog1, catalog2)
    
    # 设置匹配阈值 (5 角秒)
    max_sep = 5 * u.arcsec
    matches = sep2d < max_sep
    
    print(f"\n匹配结果 (阈值 {max_sep}):")
    print(f"  匹配成功: {matches.sum()} 对")
    print(f"  匹配率: {matches.sum() / len(catalog1) * 100:.1f}%")
    
    # 显示匹配统计
    matched_sep = sep2d[matches]
    print(f"\n匹配距离统计:")
    print(f"  平均: {matched_sep.mean().to(u.arcsec):.3f}")
    print(f"  中位数: {np.median(matched_sep).to(u.arcsec):.3f}")
    print(f"  最大: {matched_sep.max().to(u.arcsec):.3f}")
    
    return catalog1, catalog2, idx, matches


def visualize_coordinates():
    """可视化坐标"""
    print("\n" + "="*60)
    print("坐标可视化")
    print("="*60)
    
    # 创建银河系坐标数据
    np.random.seed(42)
    n_stars = 500
    
    # 模拟银河系盘结构
    l_disk = np.random.uniform(0, 360, n_stars) * u.degree
    b_disk = np.random.normal(0, 5, n_stars) * u.degree
    
    # 模拟球状星团 (高银纬)
    l_cluster = np.random.uniform(0, 360, 50) * u.degree
    b_cluster = np.random.uniform(-90, 90, 50) * u.degree
    
    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # 图1: 银道坐标系
    ax1 = axes[0]
    ax1.scatter(l_disk.value, b_disk.value, s=1, alpha=0.5, label='Disk stars')
    ax1.scatter(l_cluster.value, b_cluster.value, s=10, c='red', alpha=0.7, label='Clusters')
    ax1.set_xlabel('Galactic Longitude (l) [deg]')
    ax1.set_ylabel('Galactic Latitude (b) [deg]')
    ax1.set_title('Galactic Coordinates')
    ax1.legend()
    ax1.set_xlim(0, 360)
    ax1.set_ylim(-90, 90)
    ax1.grid(True, alpha=0.3)
    
    # 图2: Aitoff 投影
    ax2 = axes[1]
    ax2 = plt.subplot(122, projection='aitoff')
    
    # 转换为 radians 用于 Aitoff 投影
    l_rad = l_disk.to(u.radian).value
    b_rad = b_disk.to(u.radian).value
    # 调整经度范围 (-pi to pi)
    l_rad = np.where(l_rad > np.pi, l_rad - 2*np.pi, l_rad)
    
    ax2.scatter(l_rad, b_rad, s=1, alpha=0.5)
    ax2.set_xlabel('Galactic Longitude')
    ax2.set_ylabel('Galactic Latitude')
    ax2.set_title('Aitoff Projection')
    ax2.grid(True)
    
    plt.tight_layout()
    plt.savefig('coordinate_visualization.png', dpi=150, bbox_inches='tight')
    print("✓ 保存: coordinate_visualization.png")
    plt.close()


def observing_planning():
    """观测规划示例"""
    print("\n" + "="*60)
    print("观测规划")
    print("="*60)
    
    # 观测站位置 (以智利的欧洲南方天文台为例)
    paranal = EarthLocation(
        lat=-24.6272*u.deg,
        lon=-70.4044*u.deg,
        height=2635*u.m
    )
    
    # 目标列表
    targets = {
        'Sirius': SkyCoord(ra='06h45m08.9s', dec='-16d42m58s', frame='icrs'),
        'Betelgeuse': SkyCoord(ra='05h55m10.3s', dec='+07d24m25s', frame='icrs'),
        'Andromeda': SkyCoord(ra='00h42m44.3s', dec='+41d16m09s', frame='icrs'),
    }
    
    # 观测时间范围
    times = Time('2024-06-20 20:00:00') + np.linspace(0, 8, 20) * u.hour
    
    print(f"\n观测站: Paranal (ESO)")
    print(f"日期: 2024-06-20")
    print(f"\n目标可见性:")
    print(f"{'目标':<15} {'最高高度角':<12} {'最佳观测时间':<15}")
    print("-" * 45)
    
    for name, target in targets.items():
        altazs = target.transform_to(AltAz(obstime=times, location=paranal))
        
        # 找到最高高度角
        max_alt_idx = np.argmax(altazs.alt)
        max_alt = altazs.alt[max_alt_idx]
        best_time = times[max_alt_idx]
        
        # 检查是否在地平线上
        if max_alt > 30*u.deg:
            status = "✓ 良好"
        elif max_alt > 0*u.deg:
            status = "△ 可观测"
        else:
            status = "✗ 不可见"
        
        print(f"{name:<15} {max_alt:.2f}      {best_time.strftime('%H:%M')} {status}")


def main():
    """主函数"""
    print("Astropy 坐标处理示例")
    
    # 1. 创建坐标
    create_coordinates()
    
    # 2. 坐标转换
    coordinate_transformations()
    
    # 3. 距离计算
    calculate_distances()
    
    # 4. 星表匹配
    catalog_matching()
    
    # 5. 可视化
    visualize_coordinates()
    
    # 6. 观测规划
    observing_planning()
    
    print("\n" + "="*60)
    print("示例完成!")
    print("="*60)


if __name__ == '__main__':
    main()
