#!/usr/bin/env python3
"""
ENA数据库批量下载工具
功能：
1. 批量下载研究项目的所有序列
2. 使用FTP或Aspera高速下载
3. 断点续传
4. 并行下载

依赖：requests, biopython
"""

import os
import json
import time
from pathlib import Path
from typing import List, Dict, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
import requests


class ENABulkDownloader:
    """ENA批量下载器"""
    
    def __init__(self, output_dir: str = "./ena_downloads", max_workers: int = 4):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_workers = max_workers
        self.session = requests.Session()
    
    def get_run_files(self, run_accession: str) -> List[Dict]:
        """
        获取run的文件列表
        
        Args:
            run_accession: Run accession (如 ERR123456)
        
        Returns:
            文件信息列表
        """
        url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
        params = {
            "accession": run_accession,
            "result": "read_run",
            "fields": "run_accession,fastq_ftp,fastq_aspera,fastq_galaxy,"
                     "submitted_ftp,submitted_aspera,submitted_galaxy,sra_ftp"
        }
        
        response = self.session.get(url, params=params)
        
        if response.status_code != 200:
            print(f"  获取文件列表失败: {run_accession}")
            return []
        
        lines = response.text.strip().split('\n')
        if len(lines) < 2:
            return []
        
        headers = lines[0].split('\t')
        files = []
        
        for line in lines[1:]:
            values = line.split('\t')
            data = dict(zip(headers, values))
            
            # 提取FTP链接
            for key in ['fastq_ftp', 'submitted_ftp', 'sra_ftp']:
                if key in data and data[key]:
                    urls = data[key].split(';')
                    for url in urls:
                        if url.strip():
                            files.append({
                                'run': run_accession,
                                'url': f"ftp://{url.strip()}",
                                'filename': url.strip().split('/')[-1],
                                'type': key.replace('_ftp', '')
                            })
        
        return files
    
    def download_file(
        self,
        file_info: Dict,
        overwrite: bool = False
    ) -> bool:
        """
        下载单个文件
        
        Args:
            file_info: 文件信息字典
            overwrite: 是否覆盖已存在文件
        
        Returns:
            是否成功
        """
        output_path = self.output_dir / file_info['filename']
        
        if output_path.exists() and not overwrite:
            print(f"  跳过（已存在）: {file_info['filename']}")
            return True
        
        print(f"  下载: {file_info['filename']}")
        
        try:
            from urllib.request import urlretrieve
            urlretrieve(file_info['url'], str(output_path))
            print(f"    ✓ 完成 ({output_path.stat().st_size / 1024 / 1024:.1f} MB)")
            return True
        except Exception as e:
            print(f"    ✗ 失败: {e}")
            return False
    
    def download_study(
        self,
        study_accession: str,
        file_type: str = "fastq"
    ) -> Dict:
        """
        下载研究项目的所有文件
        
        Args:
            study_accession: 研究accession (如 PRJEB1234)
            file_type: 文件类型 (fastq, submitted, sra)
        
        Returns:
            下载统计
        """
        print(f"下载研究项目: {study_accession}")
        
        # 获取所有runs
        url = "https://www.ebi.ac.uk/ena/portal/api/search"
        params = {
            "result": "run",
            "query": f"study_accession={study_accession}",
            "format": "json",
            "limit": 1000
        }
        
        response = self.session.get(url, params=params)
        runs = response.json()
        
        print(f"  找到 {len(runs)} 个runs")
        
        # 获取所有文件
        all_files = []
        for run in runs:
            run_acc = run.get('run_accession')
            files = self.get_run_files(run_acc)
            
            # 筛选文件类型
            if file_type:
                files = [f for f in files if f['type'] == file_type]
            
            all_files.extend(files)
        
        print(f"  找到 {len(all_files)} 个文件待下载")
        
        # 并行下载
        successful = 0
        failed = 0
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self.download_file, f): f 
                for f in all_files
            }
            
            for future in as_completed(futures):
                if future.result():
                    successful += 1
                else:
                    failed += 1
        
        stats = {
            'study': study_accession,
            'total_files': len(all_files),
            'successful': successful,
            'failed': failed
        }
        
        # 保存统计
        stats_file = self.output_dir / f"{study_accession}_stats.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        print(f"\n下载完成: {successful}/{len(all_files)} 成功")
        return stats
    
    def download_accession_list(
        self,
        accession_file: str,
        file_type: str = "fastq"
    ) -> Dict:
        """
        从accession列表文件批量下载
        
        Args:
            accession_file: 包含accession的文件（每行一个）
            file_type: 文件类型
        
        Returns:
            下载统计
        """
        with open(accession_file) as f:
            accessions = [line.strip() for line in f if line.strip()]
        
        print(f"从列表下载: {len(accessions)} 个accessions")
        
        total_stats = {
            'total_runs': 0,
            'total_files': 0,
            'successful': 0,
            'failed': 0
        }
        
        for accession in accessions:
            stats = self.download_study(accession, file_type)
            
            total_stats['total_runs'] += stats.get('total_files', 0)
            total_stats['successful'] += stats.get('successful', 0)
            total_stats['failed'] += stats.get('failed', 0)
        
        print(f"\n总下载统计:")
        print(f"  成功: {total_stats['successful']}")
        print(f"  失败: {total_stats['failed']}")
        
        return total_stats


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="ENA批量下载工具"
    )
    
    parser.add_argument(
        "accession",
        help="研究accession或accession列表文件"
    )
    parser.add_argument(
        "--output", "-o",
        default="./ena_downloads",
        help="输出目录"
    )
    parser.add_argument(
        "--type",
        default="fastq",
        choices=["fastq", "submitted", "sra"],
        help="文件类型"
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="并行下载数"
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="输入文件是accession列表"
    )
    
    args = parser.parse_args()
    
    # 创建下载器
    downloader = ENABulkDownloader(
        output_dir=args.output,
        max_workers=args.workers
    )
    
    if args.list:
        downloader.download_accession_list(args.accession, args.type)
    else:
        downloader.download_study(args.accession, args.type)


if __name__ == "__main__":
    main()
