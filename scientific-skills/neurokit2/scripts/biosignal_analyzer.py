#!/usr/bin/env python3
"""
NeuroKit2 生物信号处理工具
用于处理和分析生理信号 (ECG, EEG, EDA, RSP, EMG, EOG)
"""

import numpy as np
from typing import Optional, Dict, List, Union, Tuple, Any
from dataclasses import dataclass
from pathlib import Path
import json

# 尝试导入 neurokit2
try:
    import neurokit2 as nk
    NEUROKIT_AVAILABLE = True
except ImportError:
    NEUROKIT_AVAILABLE = False
    print("警告: neurokit2 未安装。运行: uv pip install neurokit2")

try:
    import pandas as pd
    PANDAS_AVAILABLE = True
except ImportError:
    PANDAS_AVAILABLE = False


@dataclass
class HRVResults:
    """心率变异性分析结果"""
    time_domain: Dict[str, float]
    frequency_domain: Optional[Dict[str, float]]
    nonlinear: Optional[Dict[str, float]]
    
    def to_dict(self) -> Dict:
        return {
            'time_domain': self.time_domain,
            'frequency_domain': self.frequency_domain,
            'nonlinear': self.nonlinear
        }


class ECGAnalyzer:
    """心电信号分析器"""
    
    def __init__(self, sampling_rate: int = 1000):
        """
        初始化 ECG 分析器
        
        Args:
            sampling_rate: 采样率 (Hz)
        """
        if not NEUROKIT_AVAILABLE:
            raise ImportError("neurokit2 未安装")
        
        self.sampling_rate = sampling_rate
    
    def process(self, ecg_signal: np.ndarray) -> Tuple[pd.DataFrame, Dict]:
        """
        处理 ECG 信号
        
        Args:
            ecg_signal: ECG 信号数组
        
        Returns:
            (处理后的信号 DataFrame, 信息字典)
        """
        signals, info = nk.ecg_process(ecg_signal, sampling_rate=self.sampling_rate)
        return signals, info
    
    def detect_peaks(self, ecg_signal: np.ndarray) -> np.ndarray:
        """
        检测 R 峰
        
        Args:
            ecg_signal: ECG 信号
        
        Returns:
            R 峰位置索引
        """
        signals, info = self.process(ecg_signal)
        return info['ECG_R_Peaks']
    
    def calculate_heart_rate(self, ecg_signal: np.ndarray, 
                            window_seconds: int = 10) -> np.ndarray:
        """
        计算心率
        
        Args:
            ecg_signal: ECG 信号
            window_seconds: 窗口大小 (秒)
        
        Returns:
            心率时间序列
        """
        signals, info = self.process(ecg_signal)
        
        # 计算逐搏心率
        r_peaks = info['ECG_R_Peaks']
        rr_intervals = np.diff(r_peaks) / self.sampling_rate
        heart_rates = 60 / rr_intervals
        
        return heart_rates
    
    def analyze_hrv(self, ecg_signal: np.ndarray) -> HRVResults:
        """
        分析心率变异性 (HRV)
        
        Args:
            ecg_signal: ECG 信号
        
        Returns:
            HRV 分析结果
        """
        signals, info = self.process(ecg_signal)
        peaks = info['ECG_R_Peaks']
        
        # 时域分析
        hrv_time = nk.hrv_time(peaks, sampling_rate=self.sampling_rate)
        time_dict = {col: float(hrv_time[col].values[0]) for col in hrv_time.columns}
        
        # 频域分析
        try:
            hrv_freq = nk.hrv_frequency(peaks, sampling_rate=self.sampling_rate)
            freq_dict = {col: float(hrv_freq[col].values[0]) for col in hrv_freq.columns}
        except:
            freq_dict = None
        
        # 非线性分析
        try:
            hrv_nonlinear = nk.hrv_nonlinear(peaks, sampling_rate=self.sampling_rate)
            nonlinear_dict = {col: float(hrv_nonlinear[col].values[0]) 
                            for col in hrv_nonlinear.columns}
        except:
            nonlinear_dict = None
        
        return HRVResults(
            time_domain=time_dict,
            frequency_domain=freq_dict,
            nonlinear=nonlinear_dict
        )
    
    def plot_ecg(self, ecg_signal: np.ndarray, save_path: Optional[str] = None):
        """
        绘制 ECG 信号
        
        Args:
            ecg_signal: ECG 信号
            save_path: 保存路径
        """
        signals, info = self.process(ecg_signal)
        nk.ecg_plot(signals, info)
        
        if save_path:
            import matplotlib.pyplot as plt
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"ECG 图已保存到: {save_path}")


class EEGAnalyzer:
    """脑电信号分析器"""
    
    def __init__(self, sampling_rate: int = 250):
        """
        初始化 EEG 分析器
        
        Args:
            sampling_rate: 采样率 (Hz)
        """
        if not NEUROKIT_AVAILABLE:
            raise ImportError("neurokit2 未安装")
        
        self.sampling_rate = sampling_rate
    
    def calculate_power_bands(self, eeg_signal: np.ndarray, 
                              channels: Optional[List[str]] = None) -> Dict:
        """
        计算 EEG 频段功率
        
        Args:
            eeg_signal: EEG 信号 (单通道或多通道)
            channels: 通道名称列表
        
        Returns:
            频段功率字典
        """
        if channels is None:
            channels = ['EEG']
        
        power = nk.eeg_power(eeg_signal, sampling_rate=self.sampling_rate, 
                            channels=channels)
        
        return {col: float(power[col].values[0]) for col in power.columns}
    
    def get_band_powers(self, eeg_signal: np.ndarray) -> Dict[str, float]:
        """
        获取主要频段功率
        
        Args:
            eeg_signal: EEG 信号
        
        Returns:
            频段功率 (Delta, Theta, Alpha, Beta, Gamma)
        """
        # 计算功率谱密度
        psd = nk.signal_psd(eeg_signal, sampling_rate=self.sampling_rate)
        
        # 定义频段
        bands = {
            'Delta': (0.5, 4),
            'Theta': (4, 8),
            'Alpha': (8, 13),
            'Beta': (13, 30),
            'Gamma': (30, 100)
        }
        
        band_powers = {}
        for band_name, (low, high) in bands.items():
            # 提取频段内的功率
            mask = (psd['Frequency'] >= low) & (psd['Frequency'] <= high)
            band_powers[band_name] = float(psd.loc[mask, 'Power'].mean())
        
        return band_powers


class EDAAnalyzer:
    """皮肤电活动分析器"""
    
    def __init__(self, sampling_rate: int = 100):
        """
        初始化 EDA 分析器
        
        Args:
            sampling_rate: 采样率 (Hz)
        """
        if not NEUROKIT_AVAILABLE:
            raise ImportError("neurokit2 未安装")
        
        self.sampling_rate = sampling_rate
    
    def process(self, eda_signal: np.ndarray) -> Tuple[pd.DataFrame, Dict]:
        """
        处理 EDA 信号
        
        Args:
            eda_signal: EDA 信号
        
        Returns:
            (处理后的信号, 信息字典)
        """
        signals, info = nk.eda_process(eda_signal, sampling_rate=self.sampling_rate)
        return signals, info
    
    def detect_scr(self, eda_signal: np.ndarray) -> Dict:
        """
        检测皮肤电反应 (SCR)
        
        Args:
            eda_signal: EDA 信号
        
        Returns:
            SCR 检测结果
        """
        signals, info = self.process(eda_signal)
        
        return {
            'scr_peaks': info.get('SCR_Peaks', []),
            'scr_amplitude': info.get('SCR_Amplitude', []),
            'scr_rise_time': info.get('SCR_RiseTime', []),
            'scr_recovery': info.get('SCR_Recovery', [])
        }
    
    def get_sympathetic_index(self, eda_signal: np.ndarray) -> float:
        """
        计算交感神经指数
        
        Args:
            eda_signal: EDA 信号
        
        Returns:
            交感神经指数
        """
        signals, info = self.process(eda_signal)
        
        try:
            sympathetic = nk.eda_sympathetic(signals, sampling_rate=self.sampling_rate)
            return float(symphatetic.values[0])
        except:
            return np.nan


class RSPAnalyzer:
    """呼吸信号分析器"""
    
    def __init__(self, sampling_rate: int = 100):
        """
        初始化 RSP 分析器
        
        Args:
            sampling_rate: 采样率 (Hz)
        """
        if not NEUROKIT_AVAILABLE:
            raise ImportError("neurokit2 未安装")
        
        self.sampling_rate = sampling_rate
    
    def process(self, rsp_signal: np.ndarray) -> Tuple[pd.DataFrame, Dict]:
        """
        处理呼吸信号
        
        Args:
            rsp_signal: 呼吸信号
        
        Returns:
            (处理后的信号, 信息字典)
        """
        signals, info = nk.rsp_process(rsp_signal, sampling_rate=self.sampling_rate)
        return signals, info
    
    def get_breathing_rate(self, rsp_signal: np.ndarray) -> float:
        """
        获取呼吸频率
        
        Args:
            rsp_signal: 呼吸信号
        
        Returns:
            呼吸频率 (次/分钟)
        """
        signals, info = self.process(rsp_signal)
        
        # 从峰值计算呼吸率
        peaks = info.get('RSP_Peaks', [])
        if len(peaks) < 2:
            return np.nan
        
        duration = len(rsp_signal) / self.sampling_rate
        rate = len(peaks) / duration * 60
        
        return rate


class BiosignalSimulator:
    """生物信号模拟器"""
    
    @staticmethod
    def simulate_ecg(duration: int = 10, 
                     sampling_rate: int = 1000,
                     heart_rate: int = 70,
                     noise: float = 0.1) -> np.ndarray:
        """
        模拟 ECG 信号
        
        Args:
            duration: 持续时间 (秒)
            sampling_rate: 采样率
            heart_rate: 心率 (bpm)
            noise: 噪声水平
        
        Returns:
            模拟的 ECG 信号
        """
        if not NEUROKIT_AVAILABLE:
            raise ImportError("neurokit2 未安装")
        
        return nk.ecg_simulate(duration=duration, 
                              sampling_rate=sampling_rate,
                              heart_rate=heart_rate,
                              noise=noise)
    
    @staticmethod
    def simulate_eda(duration: int = 10,
                     sampling_rate: int = 100,
                     n_scr: int = 3) -> np.ndarray:
        """
        模拟 EDA 信号
        
        Args:
            duration: 持续时间 (秒)
            sampling_rate: 采样率
            n_scr: SCR 数量
        
        Returns:
            模拟的 EDA 信号
        """
        if not NEUROKIT_AVAILABLE:
            raise ImportError("neurokit2 未安装")
        
        return nk.eda_simulate(duration=duration,
                              sampling_rate=sampling_rate,
                              n_scr=n_scr)
    
    @staticmethod
    def simulate_rsp(duration: int = 10,
                     sampling_rate: int = 100,
                     respiratory_rate: int = 15) -> np.ndarray:
        """
        模拟呼吸信号
        
        Args:
            duration: 持续时间 (秒)
            sampling_rate: 采样率
            respiratory_rate: 呼吸频率
        
        Returns:
            模拟的呼吸信号
        """
        if not NEUROKIT_AVAILABLE:
            raise ImportError("neurokit2 未安装")
        
        return nk.rsp_simulate(duration=duration,
                              sampling_rate=sampling_rate,
                              respiratory_rate=respiratory_rate)


def demo():
    """演示功能"""
    print("=" * 60)
    print("NeuroKit2 生物信号处理演示")
    print("=" * 60)
    
    if not NEUROKIT_AVAILABLE:
        print("错误: neurokit2 未安装")
        print("请运行: uv pip install neurokit2")
        return
    
    # 1. 模拟 ECG 信号
    print("\n1. ECG 信号分析:")
    ecg_sim = BiosignalSimulator()
    ecg_signal = ecg_sim.simulate_ecg(duration=60, heart_rate=70)
    print(f"   信号长度: {len(ecg_signal)} 点")
    
    # 分析 ECG
    ecg_analyzer = ECGAnalyzer(sampling_rate=1000)
    signals, info = ecg_analyzer.process(ecg_signal)
    print(f"   检测到 R 峰: {len(info['ECG_R_Peaks'])} 个")
    
    # HRV 分析
    hrv = ecg_analyzer.analyze_hrv(ecg_signal)
    print(f"   SDNN: {hrv.time_domain.get('HRV_SDNN', 'N/A'):.2f} ms")
    print(f"   RMSSD: {hrv.time_domain.get('HRV_RMSSD', 'N/A'):.2f} ms")
    
    # 2. EDA 信号
    print("\n2. EDA 信号分析:")
    eda_signal = BiosignalSimulator.simulate_eda(duration=30, n_scr=5)
    eda_analyzer = EDAAnalyzer(sampling_rate=100)
    scr_results = eda_analyzer.detect_scr(eda_signal)
    print(f"   信号长度: {len(eda_signal)} 点")
    print(f"   检测到 SCR: {len(scr_results['scr_peaks'])} 个")
    
    # 3. 呼吸信号
    print("\n3. 呼吸信号分析:")
    rsp_signal = BiosignalSimulator.simulate_rsp(duration=60, respiratory_rate=15)
    rsp_analyzer = RSPAnalyzer(sampling_rate=100)
    br = rsp_analyzer.get_breathing_rate(rsp_signal)
    print(f"   信号长度: {len(rsp_signal)} 点")
    print(f"   呼吸频率: {br:.1f} 次/分钟")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    demo()
