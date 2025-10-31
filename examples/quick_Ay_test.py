#!/usr/bin/env python3
"""
氘核-质子Ay快速测试计算
=====================

这是一个简化的190 MeV/u氘核-质子分析能力计算，
用于快速验证程序功能和查看理论预测结果。

作者：Tic-tac团队
日期：2024年10月
"""

import os
import sys
import numpy as np
import subprocess
from pathlib import Path

class QuickAyCalculation:
    """快速Ay计算类"""
    
    def __init__(self, work_dir="output/quick_Ay_test"):
        """初始化计算"""
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # 物理参数
        self.lab_energy_per_nucleon = 190.0  # MeV/u
        self.deuteron_mass = 1875.613  # MeV/c²
        self.proton_mass = 938.272    # MeV/c²
        
        # 计算质心系能量
        self.cm_energy = self._calculate_cm_energy()
        
        # 角度网格（质心系，较少点用于快速测试）
        self.theta_cm = np.linspace(20, 160, 15)  # 15个角度点
        
        print(f"氘核-质子散射 Ay 快速测试计算")
        print(f"实验室系能量: {self.lab_energy_per_nucleon} MeV/u")
        print(f"质心系能量: {self.cm_energy:.2f} MeV")
        print(f"工作目录: {self.work_dir}")
    
    def _calculate_cm_energy(self):
        """计算质心系能量"""
        # 氘核实验室系总动能
        T_lab = 2 * self.lab_energy_per_nucleon  # 氘核总动能
        
        # 计算质心系能量
        s = (self.deuteron_mass + self.proton_mass)**2 + 2 * self.proton_mass * T_lab
        E_cm = np.sqrt(s) - self.deuteron_mass - self.proton_mass
        
        return E_cm
    
    def generate_quick_config(self):
        """生成快速测试配置文件"""
        config_content = f"""# 氘核-质子散射 Ay 快速测试配置
# 190 MeV/u 氘核 + 质子 -> 分析能力计算

# 基本参数（降低精度用于快速测试）
two_J_3N_max=1              # 三体系统最大角动量（快速测试用较小值）
J_2N_max=2                  # 二体系统最大角动量
Np_WP=30                    # p方向波包数量（减少用于快速测试）
Nq_WP=30                    # q方向波包数量

# 网格参数
Nphi=36                     # φ角积分点数
Nx=36                       # x积分点数
Np_per_WP=6                # 每个波包的p点数
Nq_per_WP=6                # 每个波包的q点数

# Chebyshev网格参数
chebyshev_s=200            # 缩放参数
chebyshev_t=1
p_grid_type=chebyshev
q_grid_type=chebyshev

# 并行设置
P123_omp_num_threads=4     # 并行线程数
parallel_run=true

# 物理选项
tensor_force=true          # 包含张量力
isospin_breaking_1S0=true  # 包含同位旋破缺
include_breakup_channels=false  # 不考虑破裂道
solve_faddeev=true
solve_dense=false

# 势能模型
potential_model=LO_internal    # 使用较简单的LO势用于快速测试

# 能量设置
energy_input_file=data/energy_190MeV_quick.txt
output_folder={self.work_dir}/results

# 其他设置
production_run=false       # 非生产运行
midpoint_approx=true       # 使用中点近似加速
"""
        
        config_file = self.work_dir / "input_quick_Ay.txt"
        with open(config_file, 'w') as f:
            f.write(config_content)
        
        print(f"快速测试配置文件生成: {config_file}")
        return config_file
    
    def generate_energy_file(self):
        """生成能量文件"""
        energy_file = Path("data/energy_190MeV_quick.txt")
        energy_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(energy_file, 'w') as f:
            f.write("# 氘核-质子散射能量点（快速测试）\n")
            f.write("# 质心系能量 (MeV)\n")
            f.write(f"{self.cm_energy:.4f}\n")
        
        print(f"能量文件生成: {energy_file}")
        return energy_file
    
    def run_quick_test(self):
        """运行快速测试计算"""
        print("\n=== 开始氘核-质子Ay快速测试 ===")
        
        # 生成必要文件
        config_file = self.generate_quick_config()
        self.generate_energy_file()
        
        # 首先尝试使用重构后的程序
        programs_to_try = ["tic-tac", "CPP/run"]
        program = None
        
        for prog in programs_to_try:
            if Path(prog).exists():
                program = prog
                break
        
        if not program:
            print("错误：未找到可执行程序")
            print("请尝试以下命令之一：")
            print("1. make                    # 编译新的重构版本")
            print("2. cd CPP && make         # 编译原版本")
            return False
        
        print(f"使用程序: {program}")
        
        # 运行计算
        try:
            print("正在运行快速Ay测试...")
            print("预计运行时间：5-15分钟")
            
            if program == "tic-tac":
                cmd = ["./tic-tac", "--input", str(config_file)]
            else:
                # 使用原版本程序
                cmd = ["./CPP/run", "--input", str(config_file)]
            
            print(f"执行命令: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)  # 30分钟超时
            
            if result.returncode == 0:
                print("✓ 快速测试计算成功完成")
                print("标准输出：")
                print(result.stdout[-1000:])  # 显示最后1000个字符
                return True
            else:
                print(f"✗ 计算失败，错误代码：{result.returncode}")
                print("标准输出：", result.stdout[-500:])
                print("错误输出：", result.stderr[-500:])
                return False
                
        except subprocess.TimeoutExpired:
            print("✗ 计算超时（超过30分钟）")
            return False
        except FileNotFoundError:
            print(f"✗ 程序文件不存在: {program}")
            return False
        except Exception as e:
            print(f"✗ 运行过程中出现错误：{e}")
            return False
    
    def analyze_results(self):
        """分析计算结果"""
        print("\n=== 分析快速测试结果 ===")
        
        results_dir = self.work_dir / "results"
        if not results_dir.exists():
            # 也检查其他可能的输出位置
            possible_dirs = [
                Path("CPP/Output"),
                Path("output"),
                self.work_dir
            ]
            
            for pdir in possible_dirs:
                if pdir.exists():
                    results_dir = pdir
                    break
        
        if not results_dir.exists():
            print("错误：未找到结果目录")
            print("可能的原因：")
            print("1. 计算未成功完成")
            print("2. 输出路径配置错误")
            return
        
        print(f"在目录中查找结果: {results_dir}")
        
        # 查找所有可能的结果文件
        result_files = []
        for pattern in ["*.txt", "*.dat", "*.csv", "*.h5"]:
            result_files.extend(list(results_dir.glob(pattern)))
        
        if result_files:
            print(f"找到 {len(result_files)} 个结果文件：")
            for f in result_files[:10]:  # 显示前10个文件
                print(f"  - {f.name}")
            if len(result_files) > 10:
                print(f"  ... 等共{len(result_files)}个文件")
                
            # 尝试找到和显示Ay相关的结果
            self._try_display_ay_results(result_files)
        else:
            print("未找到结果文件")
    
    def _try_display_ay_results(self, result_files):
        """尝试显示Ay结果"""
        
        # 查找可能包含Ay信息的文件
        ay_candidates = []
        for f in result_files:
            fname = f.name.lower()
            if any(keyword in fname for keyword in ['ay', 'analyzing', 'polarization', 'u_pw']):
                ay_candidates.append(f)
        
        if not ay_candidates:
            print("未找到明确的Ay结果文件，尝试显示其他结果...")
            ay_candidates = [f for f in result_files if f.suffix.lower() in ['.txt', '.dat', '.csv']][:3]
        
        for result_file in ay_candidates:
            print(f"\n--- 文件: {result_file.name} ---")
            try:
                with open(result_file, 'r') as f:
                    lines = f.readlines()
                
                # 显示前几行了解文件结构
                print("文件内容预览：")
                for i, line in enumerate(lines[:15]):
                    print(f"{i+1:2d}: {line.rstrip()}")
                
                if len(lines) > 15:
                    print(f"... (文件共{len(lines)}行)")
                
                # 尝试解析数值数据
                self._try_parse_numerical_data(result_file)
                
            except Exception as e:
                print(f"读取文件 {result_file.name} 时出错：{e}")
    
    def _try_parse_numerical_data(self, file_path):
        """尝试解析数值数据"""
        try:
            # 尝试读取数值数据
            data = np.loadtxt(file_path, comments='#')
            
            if data.size == 0:
                return
            
            if data.ndim == 1:
                print(f"一维数据，{len(data)}个数值")
                if len(data) <= 20:
                    print("数值：", data)
                else:
                    print("前10个数值：", data[:10])
                    print("后10个数值：", data[-10:])
            
            elif data.ndim == 2:
                print(f"二维数据，形状：{data.shape}")
                print("前几行数据：")
                print(data[:min(10, data.shape[0])])
                
                # 如果看起来像角度-Ay数据
                if data.shape[1] >= 2:
                    angles = data[:, 0]
                    values = data[:, 1]
                    print(f"\n可能的物理量分析：")
                    print(f"第1列范围：{angles.min():.3f} 到 {angles.max():.3f}")
                    print(f"第2列范围：{values.min():.6f} 到 {values.max():.6f}")
                    
                    if 0 <= angles.min() <= 180 and angles.max() <= 180:
                        print("第1列可能是散射角（度）")
                    
                    if -1 <= values.min() <= 1 and values.max() <= 1:
                        print("第2列可能是分析能力Ay（范围在-1到1之间）")
        
        except Exception as e:
            print(f"数值解析失败：{e}")
    
    def create_simple_plot(self):
        """创建简单的结果图"""
        try:
            import matplotlib.pyplot as plt
            
            results_dir = self.work_dir / "results"
            
            # 查找可能的数据文件
            data_files = list(results_dir.glob("*.txt")) + list(results_dir.glob("*.dat"))
            
            fig, axes = plt.subplots(1, min(2, len(data_files)), figsize=(12, 5))
            if len(data_files) == 1:
                axes = [axes]
            
            plotted = False
            for i, data_file in enumerate(data_files[:2]):
                try:
                    data = np.loadtxt(data_file, comments='#')
                    if data.ndim == 2 and data.shape[1] >= 2:
                        axes[i].plot(data[:, 0], data[:, 1], 'o-')
                        axes[i].set_title(f"{data_file.name}")
                        axes[i].grid(True, alpha=0.3)
                        plotted = True
                except:
                    continue
            
            if plotted:
                plt.tight_layout()
                plot_file = self.work_dir / "quick_results_plot.png"
                plt.savefig(plot_file, dpi=150, bbox_inches='tight')
                plt.show()
                print(f"结果图保存为: {plot_file}")
            else:
                print("无法创建图形：未找到合适的数据")
                
        except ImportError:
            print("matplotlib未安装，跳过绘图")
        except Exception as e:
            print(f"绘图时出错：{e}")

def main():
    """主函数"""
    print("=" * 60)
    print("氘核-质子散射分析能力(Ay)快速测试")
    print("能量: 190 MeV/u")
    print("=" * 60)
    
    # 创建计算实例
    calc = QuickAyCalculation()
    
    # 直接运行快速测试
    success = calc.run_quick_test()
    
    if success:
        calc.analyze_results()
        calc.create_simple_plot()
        print("\n✓ 快速测试完成！")
        print(f"结果保存在: {calc.work_dir}")
    else:
        print("\n✗ 快速测试失败")
        print("可能的解决方案：")
        print("1. 检查程序是否正确编译")
        print("2. 确认所需的库文件是否安装")
        print("3. 查看错误输出信息")

if __name__ == "__main__":
    main()