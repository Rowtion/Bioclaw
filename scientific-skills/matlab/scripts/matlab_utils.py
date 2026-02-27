#!/usr/bin/env python3
"""
MATLAB/Octave 集成工具
提供 Python 与 MATLAB/Octave 的交互功能
"""

import subprocess
import os
import tempfile
from pathlib import Path
from typing import Optional, List, Dict, Any, Union


class MatlabEngine:
    """MATLAB/Octave 引擎封装"""
    
    def __init__(self, use_octave: bool = True):
        """
        初始化引擎
        
        Args:
            use_octave: 是否使用 GNU Octave (免费)，False 则使用 MATLAB
        """
        self.use_octave = use_octave
        self.engine_name = "Octave" if use_octave else "MATLAB"
        self._check_engine()
    
    def _check_engine(self):
        """检查引擎是否可用"""
        cmd = "octave" if self.use_octave else "matlab"
        try:
            subprocess.run([cmd, "--version"], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise RuntimeError(f"{self.engine_name} 未安装或不在 PATH 中")
    
    def run_script(self, script_path: str, *args) -> str:
        """
        运行 MATLAB/Octave 脚本
        
        Args:
            script_path: .m 脚本文件路径
            *args: 传递给脚本的参数
        
        Returns:
            脚本输出
        """
        if self.use_octave:
            cmd = ["octave", script_path]
        else:
            arg_str = ",".join(str(a) for a in args)
            cmd = ["matlab", "-nodisplay", "-nosplash", "-r", 
                   f"run('{script_path}'); exit;"]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        return result.stdout + result.stderr
    
    def run_code(self, code: str) -> str:
        """
        直接运行 MATLAB/Octave 代码
        
        Args:
            code: MATLAB/Octave 代码字符串
        
        Returns:
            执行输出
        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.m', delete=False) as f:
            f.write(code)
            temp_path = f.name
        
        try:
            return self.run_script(temp_path)
        finally:
            os.unlink(temp_path)
    
    def matrix_to_matlab(self, matrix: List[List[float]], var_name: str = "A") -> str:
        """
        将 Python 矩阵转换为 MATLAB 矩阵赋值语句
        
        Args:
            matrix: 二维数值列表
            var_name: 变量名
        
        Returns:
            MATLAB 赋值语句
        """
        rows_str = []
        for row in matrix:
            row_str = " ".join(str(x) for x in row)
            rows_str.append(row_str)
        matlab_str = "; ".join(rows_str)
        return f"{var_name} = [{matlab_str}];"
    
    def solve_linear_system(self, A: List[List[float]], b: List[float]) -> Optional[List[float]]:
        """
        使用 MATLAB/Octave 求解线性方程组 Ax = b
        
        Args:
            A: 系数矩阵
            b: 右侧向量
        
        Returns:
            解向量 x
        """
        code = f"""
A = {A};
b = {b}';
x = A \\ b;
disp(x);
"""
        output = self.run_code(code)
        # 解析输出中的数值
        lines = output.strip().split('\n')
        for line in lines:
            line = line.strip()
            if line and not line.startswith('%') and not line.startswith('warning'):
                try:
                    values = [float(x) for x in line.split()]
                    if values:
                        return values
                except ValueError:
                    continue
        return None
    
    def eigendecomposition(self, matrix: List[List[float]]) -> Dict[str, Any]:
        """
        计算矩阵的特征值和特征向量
        
        Args:
            matrix: 方阵
        
        Returns:
            {'eigenvalues': [...], 'eigenvectors': [[...], ...]}
        """
        code = f"""
A = {matrix};
[V, D] = eig(A);
disp('Eigenvalues:');
disp(diag(D)');
disp('Eigenvectors:');
disp(V);
"""
        output = self.run_code(code)
        # 解析输出
        return {"output": output, "matrix": matrix}


class MatlabCodeGenerator:
    """MATLAB/Octave 代码生成器"""
    
    @staticmethod
    def generate_matrix_operations(matrix_size: int = 3) -> str:
        """生成矩阵操作示例代码"""
        return f"""% 矩阵操作示例
% 创建矩阵
A = randn({matrix_size}, {matrix_size});
B = eye({matrix_size});

% 基本运算
disp('矩阵 A:');
disp(A);
disp('转置 A\':');
disp(A');
disp('逆矩阵 inv(A):');
disp(inv(A));
disp('行列式 det(A):');
disp(det(A));
disp('特征值 eig(A):');
disp(eig(A));

% 解线性方程组
x = A \\ B(:,1);
disp('解 x:');
disp(x);
"""
    
    @staticmethod
    def generate_plotting(data_length: int = 100) -> str:
        """生成绘图示例代码"""
        return f"""% 绘图示例
x = linspace(0, 2*pi, {data_length});
y1 = sin(x);
y2 = cos(x);

% 创建图形
figure;
plot(x, y1, 'b-', 'LineWidth', 2);
hold on;
plot(x, y2, 'r--', 'LineWidth', 2);
hold off;

xlabel('x');
ylabel('y');
title('Sine and Cosine Waves');
legend('sin(x)', 'cos(x)');
grid on;

% 保存图形
saveas(gcf, 'plot_output.png');
disp('图形已保存为 plot_output.png');
"""
    
    @staticmethod
    def generate_statistics(data_vector: List[float] = None) -> str:
        """生成统计分析示例代码"""
        if data_vector is None:
            data_str = "randn(100, 1)"
        else:
            data_str = str(data_vector)
        
        return f"""% 统计分析示例
data = {data_str};

% 描述性统计
disp('均值:');
disp(mean(data));
disp('标准差:');
disp(std(data));
disp('中位数:');
disp(median(data));
disp('最小值:');
disp(min(data));
disp('最大值:');
disp(max(data));

% 直方图
figure;
histogram(data, 20);
title('数据分布');
xlabel('值');
ylabel('频率');
saveas(gcf, 'histogram.png');
"""
    
    @staticmethod
    def generate_ode_solver() -> str:
        """生成微分方程求解示例代码"""
        return """% 常微分方程求解示例
% 定义方程: dy/dt = -2*y, y(0) = 1
f = @(t, y) -2*y;

% 求解
[t, y] = ode45(f, [0 5], 1);

% 绘图
figure;
plot(t, y, 'b-', 'LineWidth', 2);
xlabel('时间 t');
ylabel('y(t)');
title('ODE Solution: dy/dt = -2y');
grid on;
saveas(gcf, 'ode_solution.png');

% 解析解对比
y_exact = exp(-2*t);
hold on;
plot(t, y_exact, 'r--', 'LineWidth', 2);
legend('数值解', '解析解');
disp('最大误差:');
disp(max(abs(y - y_exact')));
"""


def convert_python_to_matlab_array(py_array: Union[List, Any]) -> str:
    """
    将 Python 数组/列表转换为 MATLAB 数组字符串
    
    Args:
        py_array: Python 列表或 numpy 数组
    
    Returns:
        MATLAB 数组字符串表示
    """
    try:
        import numpy as np
        if isinstance(py_array, np.ndarray):
            py_array = py_array.tolist()
    except ImportError:
        pass
    
    if not isinstance(py_array, list):
        return str(py_array)
    
    if len(py_array) == 0:
        return "[]"
    
    # 检查是否为二维数组
    if isinstance(py_array[0], list):
        rows = []
        for row in py_array:
            row_str = ", ".join(str(x) for x in row)
            rows.append(row_str)
        return "[\n    " + ";\n    ".join(rows) + "\n]"
    else:
        # 一维数组
        return "[" + ", ".join(str(x) for x in py_array) + "]"


def execute_matlab_expression(expression: str, use_octave: bool = True) -> str:
    """
    执行单个 MATLAB/Octave 表达式
    
    Args:
        expression: MATLAB 表达式
        use_octave: 是否使用 Octave
    
    Returns:
        执行结果输出
    """
    engine = MatlabEngine(use_octave=use_octave)
    return engine.run_code(expression)


if __name__ == "__main__":
    # 示例用法
    print("MATLAB/Octave 工具示例\n")
    
    # 生成示例代码
    generator = MatlabCodeGenerator()
    
    print("=== 矩阵操作示例 ===")
    print(generator.generate_matrix_operations(3))
    
    print("\n=== 绘图示例 ===")
    print(generator.generate_plotting(100))
    
    print("\n=== ODE 求解示例 ===")
    print(generator.generate_ode_solver())
    
    # Python 到 MATLAB 数组转换
    print("\n=== Python 到 MATLAB 数组转换 ===")
    python_matrix = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    matlab_str = convert_python_to_matlab_array(python_matrix)
    print(f"Python: {python_matrix}")
    print(f"MATLAB: A = {matlab_str}")
