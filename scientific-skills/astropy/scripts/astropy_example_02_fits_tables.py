#!/usr/bin/env python3
"""
astropy_example_02_fits_tables.py
FITS 文件和表格处理

本脚本演示:
1. 创建和读取 FITS 文件
2. 天文表格操作
3. 与 Astropy Tables 协作
4. 数据筛选和分析
"""

import numpy as np
from astropy.io import fits
from astropy.table import Table, QTable, join, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import matplotlib.pyplot as plt
from pathlib import Path


def create_fits_file(output_path='./sample_catalog.fits'):
    """创建示例 FITS 文件"""
    print("="*60)
    print("创建 FITS 文件")
    print("="*60)
    
    np.random.seed(42)
    n_sources = 1000
    
    # 创建模拟的星表数据
    data = {
        'ID': np.arange(1, n_sources + 1),
        'RA': np.random.uniform(10, 15, n_sources),  # 度
        'DEC': np.random.uniform(30, 35, n_sources),  # 度
        'MAG_R': np.random.normal(20, 1.5, n_sources),  # r波段星等
        'MAG_G': np.random.normal(20.5, 1.5, n_sources),  # g波段星等
        'REDSHIFT': np.random.exponential(0.3, n_sources),
        'CLASS': np.random.choice(['STAR', 'GALAXY', 'QSO'], n_sources, p=[0.3, 0.6, 0.1]),
    }
    
    # 创建 Astropy Table
    table = Table(data)
    
    # 添加单位
    table['RA'].unit = u.degree
    table['DEC'].unit = u.degree
    table['MAG_R'].unit = u.mag
    table['MAG_G'].unit = u.mag
    
    # 计算颜色
    table['G_R'] = table['MAG_G'] - table['MAG_R']
    
    # 添加描述性的元数据
    table.meta['CREATOR'] = 'Astropy Example Script'
    table.meta['DATE'] = Time.now().iso
    table.meta['COMMENT'] = 'Simulated astronomical catalog for demonstration'
    
    # 保存为 FITS 文件
    table.write(output_path, format='fits', overwrite=True)
    print(f"✓ 创建 FITS 文件: {output_path}")
    print(f"  包含 {len(table)} 个源, {len(table.columns)} 列")
    
    return table


def read_fits_file(file_path='./sample_catalog.fits'):
    """读取 FITS 文件"""
    print("\n" + "="*60)
    print("读取 FITS 文件")
    print("="*60)
    
    # 方法1: 使用 astropy.table
    table = Table.read(file_path, format='fits')
    print(f"\n方法1 - Table.read:")
    print(f"  行数: {len(table)}")
    print(f"  列: {list(table.columns)}")
    
    # 显示元数据
    print(f"\n  元数据:")
    for key, value in table.meta.items():
        print(f"    {key}: {value}")
    
    # 方法2: 使用 astropy.io.fits (底层访问)
    with fits.open(file_path) as hdul:
        print(f"\n方法2 - fits.open:")
        print(f"  HDU 列表:")
        hdul.info()
        
        # 访问数据
        data = hdul[1].data  # 第一个扩展 HDU 通常包含表格
        print(f"\n  前5行数据:")
        print(data[:5])
    
    return table


def table_operations(table):
    """表格操作"""
    print("\n" + "="*60)
    print("表格操作")
    print("="*60)
    
    # 1. 筛选数据
    print("\n1. 筛选星系且 r < 19:")
    galaxies_bright = table[(table['CLASS'] == 'GALAXY') & (table['MAG_R'] < 19)]
    print(f"   找到 {len(galaxies_bright)} 个源")
    
    # 2. 选择特定列
    print("\n2. 选择特定列:")
    subset = table[['ID', 'RA', 'DEC', 'REDSHIFT']][:5]
    print(subset)
    
    # 3. 排序
    print("\n3. 按红移排序 (前5):")
    sorted_table = table.sort_values('REDSHIFT', reverse=True)[:5]
    print(sorted_table[['ID', 'REDSHIFT', 'MAG_R']])
    
    # 4. 分组统计
    print("\n4. 按类型分组统计:")
    for class_name in np.unique(table['CLASS']):
        subset = table[table['CLASS'] == class_name]
        print(f"   {class_name}:")
        print(f"     数量: {len(subset)}")
        print(f"     平均 r 星等: {np.mean(subset['MAG_R']):.2f}")
        print(f"     平均红移: {np.mean(subset['REDSHIFT']):.3f}")
    
    # 5. 添加计算列
    print("\n5. 添加坐标对象列:")
    coords = SkyCoord(ra=table['RA'], dec=table['DEC'], unit=u.degree, frame='icrs')
    table['GAL_L'] = coords.galactic.l.degree
    table['GAL_B'] = coords.galactic.b.degree
    print("   已添加银河系坐标列 (GAL_L, GAL_B)")


def advanced_table_operations():
    """高级表格操作"""
    print("\n" + "="*60)
    print("高级表格操作")
    print("="*60)
    
    np.random.seed(42)
    
    # 创建两个示例表格
    table1 = Table({
        'ID': [1, 2, 3, 4, 5],
        'NAME': ['Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon'],
        'MAG': [10.5, 12.3, 11.8, 13.2, 9.8],
    })
    
    table2 = Table({
        'ID': [1, 2, 3, 4, 6],  # 注意: 6 是新的
        'REDSHIFT': [0.1, 0.5, 0.3, 0.8, 1.2],
        'TYPE': ['GAL', 'QSO', 'GAL', 'GAL', 'STAR'],
    })
    
    # 1. 表格连接 (Join)
    print("\n1. 表格连接 (内连接):")
    joined = join(table1, table2, keys='ID', join_type='inner')
    print(joined)
    
    print("\n   表格连接 (外连接):")
    joined_outer = join(table1, table2, keys='ID', join_type='outer')
    print(joined_outer)
    
    # 2. 表格堆叠 (Stack)
    print("\n2. 表格堆叠:")
    table3 = Table({
        'ID': [7, 8, 9],
        'NAME': ['Zeta', 'Eta', 'Theta'],
        'MAG': [14.1, 15.2, 11.5],
    })
    stacked = vstack([table1, table3])
    print(f"   堆叠后行数: {len(stacked)}")
    
    # 3. 使用 QTable (带单位的表格)
    print("\n3. QTable (带单位):")
    qtable = QTable({
        'wavelength': [400, 500, 600, 700] * u.nm,
        'flux': [1.2, 1.5, 1.3, 1.0] * u.erg / u.cm**2 / u.s / u.AA,
    })
    print(qtable)
    
    # 单位转换
    qtable['flux_Jy'] = qtable['flux'].to(u.Jy, 
        equivalencies=u.spectral_density(qtable['wavelength']))
    print("\n   转换后的流量 (Jansky):")
    print(qtable[['wavelength', 'flux_Jy']])


def visualize_catalog(table):
    """可视化星表"""
    print("\n" + "="*60)
    print("星表可视化")
    print("="*60)
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. 天空分布
    ax1 = axes[0, 0]
    scatter = ax1.scatter(table['RA'], table['DEC'], 
                         c=table['MAG_R'], cmap='viridis_r',
                         s=5, alpha=0.6)
    ax1.set_xlabel('RA (deg)')
    ax1.set_ylabel('DEC (deg)')
    ax1.set_title('Sky Distribution')
    ax1.invert_xaxis()  # 天文惯例: RA 从右向左
    plt.colorbar(scatter, ax=ax1, label='r magnitude')
    
    # 2. 颜色-星等图 (Color-Magnitude Diagram)
    ax2 = axes[0, 1]
    for class_name, color in zip(['STAR', 'GALAXY', 'QSO'], ['red', 'blue', 'green']):
        mask = table['CLASS'] == class_name
        ax2.scatter(table['G_R'][mask], table['MAG_R'][mask],
                   c=color, label=class_name, s=5, alpha=0.5)
    ax2.set_xlabel('g - r')
    ax2.set_ylabel('r magnitude')
    ax2.set_title('Color-Magnitude Diagram')
    ax2.invert_yaxis()  # 亮星等在下方
    ax2.legend()
    
    # 3. 红移分布
    ax3 = axes[1, 0]
    ax3.hist(table['REDSHIFT'], bins=30, alpha=0.7, color='steelblue', edgecolor='black')
    ax3.set_xlabel('Redshift')
    ax3.set_ylabel('Count')
    ax3.set_title('Redshift Distribution')
    
    # 4. 星等分布 (按类型)
    ax4 = axes[1, 1]
    for class_name in np.unique(table['CLASS']):
        mask = table['CLASS'] == class_name
        ax4.hist(table['MAG_R'][mask], bins=20, alpha=0.5, 
                label=class_name, edgecolor='black')
    ax4.set_xlabel('r magnitude')
    ax4.set_ylabel('Count')
    ax4.set_title('Magnitude Distribution by Class')
    ax4.legend()
    
    plt.tight_layout()
    plt.savefig('catalog_visualization.png', dpi=150, bbox_inches='tight')
    print("✓ 保存: catalog_visualization.png")
    plt.close()


def export_data(table, output_dir='./output'):
    """导出数据为不同格式"""
    print("\n" + "="*60)
    print("导出数据")
    print("="*60)
    
    Path(output_dir).mkdir(exist_ok=True)
    
    # 1. 保存为 FITS
    fits_path = Path(output_dir) / 'catalog_export.fits'
    table.write(fits_path, format='fits', overwrite=True)
    print(f"✓ FITS 格式: {fits_path}")
    
    # 2. 保存为 CSV
    csv_path = Path(output_dir) / 'catalog_export.csv'
    table.write(csv_path, format='csv', overwrite=True)
    print(f"✓ CSV 格式: {csv_path}")
    
    # 3. 保存为 ECSV (增强 CSV，保留元数据)
    ecsv_path = Path(output_dir) / 'catalog_export.ecsv'
    table.write(ecsv_path, format='ascii.ecsv', overwrite=True)
    print(f"✓ ECSV 格式: {ecsv_path}")
    
    # 4. 保存为 HDF5
    h5_path = Path(output_dir) / 'catalog_export.h5'
    table.write(h5_path, format='hdf5', path='data', overwrite=True)
    print(f"✓ HDF5 格式: {h5_path}")


def main():
    """主函数"""
    print("Astropy FITS 和表格处理示例")
    
    # 1. 创建 FITS 文件
    table = create_fits_file('./sample_catalog.fits')
    
    # 2. 读取 FITS 文件
    table = read_fits_file('./sample_catalog.fits')
    
    # 3. 表格操作
    table_operations(table)
    
    # 4. 高级表格操作
    advanced_table_operations()
    
    # 5. 可视化
    visualize_catalog(table)
    
    # 6. 导出数据
    export_data(table)
    
    print("\n" + "="*60)
    print("示例完成!")
    print("="*60)


if __name__ == '__main__':
    main()
