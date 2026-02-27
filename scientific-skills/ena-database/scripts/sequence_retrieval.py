#!/usr/bin/env python3
"""
ENA数据库序列检索脚本
功能：
1. 通过accession检索序列数据
2. 下载FASTA/FASTQ文件
3. 搜索特定条件的序列
4. 批量下载研究数据

依赖：requests, biopython
"""

import requests
import time
from pathlib import Path
from typing import List, Dict, Optional, Union
from Bio import Entrez

# ENA API endpoints
ENA_PORTAL_API = "https://www.ebi.ac.uk/ena/portal/api"
ENA_BROWSER_API = "https://www.ebi.ac.uk/ena/browser/api"


def search_ena(
    query: str,
    result_type: str = "sequence",
    format: str = "json",
    limit: int = 100
) -> List[Dict]:
    """
    搜索ENA数据库
    
    Args:
        query: 搜索查询（如 tax_tree(9606) 或 study_accession=PRJEB1234）
        result_type: 结果类型 (sequence, sample, study, run, assembly)
        format: 返回格式 (json, tsv, xml)
        limit: 最大返回数量
    
    Returns:
        搜索结果列表
    """
    url = f"{ENA_PORTAL_API}/search"
    params = {
        "result": result_type,
        "query": query,
        "format": format,
        "limit": limit
    }
    
    print(f"搜索ENA: {query}")
    print(f"  结果类型: {result_type}")
    
    response = requests.get(url, params=params)
    response.raise_for_status()
    
    if format == "json":
        results = response.json()
    else:
        results = response.text
    
    print(f"  找到 {len(results) if isinstance(results, list) else 'N/A'} 条记录")
    return results


def fetch_sequence_by_accession(
    accession: str,
    format: str = "fasta",
    output_file: Optional[str] = None
) -> str:
    """
    通过accession获取序列
    
    Args:
        accession: 序列accession号 (如 LR794134)
        format: 返回格式 (fasta, text, xml)
        output_file: 输出文件路径（可选）
    
    Returns:
        序列数据或文件路径
    """
    url = f"{ENA_BROWSER_API}/{format}/{accession}"
    
    print(f"获取序列: {accession}")
    
    response = requests.get(url)
    response.raise_for_status()
    
    sequence_data = response.text
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(sequence_data)
        print(f"  ✓ 已保存: {output_file}")
        return output_file
    
    return sequence_data


def fetch_study_data(
    study_accession: str,
    result_type: str = "sample",
    output_dir: str = "./ena_data"
) -> List[Dict]:
    """
    获取研究的所有相关数据
    
    Args:
        study_accession: 研究accession (如 PRJEB1234)
        result_type: 数据类型 (sample, run, experiment)
        output_dir: 输出目录
    
    Returns:
        数据记录列表
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"获取研究数据: {study_accession}")
    
    query = f"study_accession={study_accession}"
    results = search_ena(query, result_type=result_type, limit=1000)
    
    # 保存元数据
    import json
    metadata_file = Path(output_dir) / f"{study_accession}_{result_type}.json"
    with open(metadata_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"  ✓ 元数据已保存: {metadata_file}")
    print(f"  找到 {len(results)} 个{result_type}")
    
    return results


def download_fastq_files(
    run_accessions: List[str],
    output_dir: str = "./fastq_data",
    aspera: bool = False
) -> List[str]:
    """
    下载FASTQ文件
    
    Args:
        run_accessions: Run accession列表 (如 ERR123456)
        output_dir: 输出目录
        aspera: 是否使用Aspera高速下载
    
    Returns:
        下载的文件路径列表
    """
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"下载FASTQ文件: {len(run_accessions)} 个runs")
    
    downloaded_files = []
    
    for accession in run_accessions:
        # 获取文件位置
        url = f"{ENA_PORTAL_API}/filereport"
        params = {
            "accession": accession,
            "result": "read_run",
            "fields": "run_accession,fastq_ftp,fastq_aspera,submitted_ftp"
        }
        
        response = requests.get(url, params=params)
        
        if response.status_code == 200:
            lines = response.text.strip().split('\n')
            if len(lines) > 1:
                # 解析TSV
                headers = lines[0].split('\t')
                data = lines[1].split('\t')
                
                # 获取FTP链接
                ftp_urls = []
                for line in data:
                    if 'ftp.sra.ebi.ac.uk' in line:
                        ftp_urls.extend(line.split(';'))
                
                for ftp_url in ftp_urls:
                    if ftp_url.strip():
                        filename = ftp_url.split('/')[-1]
                        output_file = Path(output_dir) / filename
                        
                        if output_file.exists():
                            print(f"  跳过（已存在）: {filename}")
                            downloaded_files.append(str(output_file))
                            continue
                        
                        # 下载文件
                        print(f"  下载: {filename}")
                        try:
                            from urllib.request import urlretrieve
                            urlretrieve(f"ftp://{ftp_url}", str(output_file))
                            downloaded_files.append(str(output_file))
                            print(f"    ✓ 完成")
                        except Exception as e:
                            print(f"    ✗ 失败: {e}")
        
        # 限速
        time.sleep(0.5)
    
    print(f"\n下载完成: {len(downloaded_files)} 个文件")
    return downloaded_files


def get_taxonomy_info(taxon_id: Union[int, str]) -> Dict:
    """
    获取分类学信息
    
    Args:
        taxon_id: NCBI taxonomy ID
    
    Returns:
        分类学信息字典
    """
    url = f"https://www.ebi.ac.uk/ena/taxonomy/rest/tax-id/{taxon_id}"
    
    print(f"获取分类信息: taxon_id={taxon_id}")
    
    response = requests.get(url)
    response.raise_for_status()
    
    data = response.json()
    
    print(f"  物种: {data.get('scientificName', 'N/A')}")
    print(f"  等级: {data.get('rank', 'N/A')}")
    
    return data


def fetch_assembly(
    assembly_accession: str,
    format: str = "fasta",
    output_file: Optional[str] = None
) -> str:
    """
    获取基因组组装
    
    Args:
        assembly_accession: 组装accession (如 GCA_000001405)
        format: 格式 (fasta, xml)
        output_file: 输出文件路径
    
    Returns:
        组装数据或文件路径
    """
    # 搜索组装信息
    query = f"assembly_accession={assembly_accession}"
    results = search_ena(query, result_type="assembly", format="json")
    
    if not results:
        print(f"未找到组装: {assembly_accession}")
        return ""
    
    # 获取序列文件链接
    assembly = results[0]
    sequence_url = assembly.get('chromosome_sequence', '')
    
    if sequence_url and output_file:
        print(f"下载组装序列: {assembly_accession}")
        print(f"  从: {sequence_url}")
        
        # 下载
        try:
            from urllib.request import urlretrieve
            urlretrieve(sequence_url, output_file)
            print(f"  ✓ 已保存: {output_file}")
            return output_file
        except Exception as e:
            print(f"  ✗ 下载失败: {e}")
            return ""
    
    return str(results)


def main():
    """主函数"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="ENA数据库序列检索工具"
    )
    subparsers = parser.add_subparsers(dest="command", help="可用命令")
    
    # 搜索命令
    search_parser = subparsers.add_parser("search", help="搜索ENA")
    search_parser.add_argument("query", help="搜索查询")
    search_parser.add_argument("--type", default="sequence", 
                              choices=["sequence", "sample", "study", "run", "assembly"],
                              help="结果类型")
    search_parser.add_argument("--limit", type=int, default=100)
    search_parser.add_argument("--output", "-o", help="输出JSON文件")
    
    # 获取序列
    fetch_parser = subparsers.add_parser("fetch", help="获取序列")
    fetch_parser.add_argument("accession", help="序列accession")
    fetch_parser.add_argument("--format", default="fasta", choices=["fasta", "text", "xml"])
    fetch_parser.add_argument("--output", "-o", help="输出文件")
    
    # 获取研究数据
    study_parser = subparsers.add_parser("study", help="获取研究数据")
    study_parser.add_argument("accession", help="研究accession (如 PRJEB1234)")
    study_parser.add_argument("--type", default="sample",
                             choices=["sample", "run", "experiment"])
    study_parser.add_argument("--output", "-o", default="./ena_data")
    
    # 下载FASTQ
    fastq_parser = subparsers.add_parser("download-fastq", help="下载FASTQ文件")
    fastq_parser.add_argument("accessions", nargs="+", help="Run accessions")
    fastq_parser.add_argument("--output", "-o", default="./fastq_data")
    
    # 分类学信息
    tax_parser = subparsers.add_parser("taxonomy", help="获取分类学信息")
    tax_parser.add_argument("taxon_id", help="NCBI taxonomy ID")
    
    # 获取组装
    assembly_parser = subparsers.add_parser("assembly", help="获取基因组组装")
    assembly_parser.add_argument("accession", help="组装accession (如 GCA_000001405)")
    assembly_parser.add_argument("--output", "-o", help="输出文件")
    
    args = parser.parse_args()
    
    if args.command == "search":
        results = search_ena(args.query, args.type, limit=args.limit)
        
        if args.output:
            import json
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            print(f"结果已保存: {args.output}")
        else:
            for r in results[:10]:
                print(r)
    
    elif args.command == "fetch":
        result = fetch_sequence_by_accession(
            args.accession,
            args.format,
            args.output
        )
        if not args.output:
            print(result[:1000] + "..." if len(result) > 1000 else result)
    
    elif args.command == "study":
        fetch_study_data(args.accession, args.type, args.output)
    
    elif args.command == "download-fastq":
        download_fastq_files(args.accessions, args.output)
    
    elif args.command == "taxonomy":
        info = get_taxonomy_info(args.taxon_id)
        print(json.dumps(info, indent=2))
    
    elif args.command == "assembly":
        fetch_assembly(args.accession, output_file=args.output)
    
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
