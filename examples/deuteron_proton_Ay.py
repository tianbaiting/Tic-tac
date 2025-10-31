#!/usr/bin/env python3
"""
氘核-质子散射分析能力(Ay)计算示例
=====================================

计算190 MeV/u极化氘核与质子散射的分析能力(Ay)
Ay反映散射截面随y方向极化的变化：
Ay = (σ↑ - σ↓) / (σ↑ + σ↓)

物理背景：
- 入射粒子：极化氘核，实验室系能量 190 MeV/u
- 靶核：质子（静止）
- 观测量：y方向分析能力 Ay(θ)
- 角度范围：0° - 180°（质心系）

作者：Tic-tac团队
日期：2024年10月
"""

import os
import sys
import numpy as np
import subprocess
from pathlib import Path

class DeuteronProtonAyCalculation:
    """氘核-质子Ay计算类"""
    
    def __init__(self, work_dir="output/deuteron_proton_Ay"):
        """初始化计算"""
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        # 物理参数
        self.lab_energy_per_nucleon = 190.0  # MeV/u
        self.deuteron_mass = 1875.613  # MeV/c²
        self.proton_mass = 938.272    # MeV/c²
        
        # 计算质心系能量
        self.cm_energy = self._calculate_cm_energy()
        
        # 角度网格（质心系）
        self.theta_cm = np.linspace(10, 170, 33)  # 度
        
        print(f"氘核-质子散射 Ay 计算")
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
    
    def generate_config_file(self):
        """生成Tic-tac配置文件"""
        config_content = f"""# 氘核-质子散射 Ay 计算配置
# 190 MeV/u 氘核 + 质子 -> 分析能力计算

# 基本参数
two_J_3N_max=3              # 三体系统最大角动量（需要较高值用于极化计算）
J_2N_max=3                  # 二体系统最大角动量
Np_WP=80                    # p方向波包数量（需要足够精度）
Nq_WP=80                    # q方向波包数量

# 网格参数
Nphi=72                     # φ角积分点数（极化计算需要更多点）
Nx=72                       # x积分点数
Np_per_WP=12               # 每个波包的p点数
Nq_per_WP=12               # 每个波包的q点数

# Chebyshev网格参数
chebyshev_s=300            # 高能需要更大的缩放参数
chebyshev_t=1
p_grid_type=chebyshev
q_grid_type=chebyshev

# 并行设置
P123_omp_num_threads=8     # 使用多线程加速
parallel_run=true

# 物理选项
tensor_force=true          # 必须包含张量力
isospin_breaking_1S0=true  # 包含同位旋破缺
include_breakup_channels=false  # 先不考虑破裂道
solve_faddeev=true
solve_dense=false

# 极化相关设置
calculate_polarization=true    # 计算极化观测量
polarization_type=vector       # 矢量极化（氘核的极化）
analyzing_power=true           # 计算分析能力

# 势能模型（推荐用于中高能）
potential_model=N3LO_Idaho    # N3LO手征势

# 能量设置
energy_input_file=data/energy_190MeV_deuteron.txt
output_folder={self.work_dir}/results

# 角度设置
angle_input_file=data/angles_cm.txt
coordinate_system=center_of_mass

# 数据存储
save_matrix_elements=true
save_phase_shifts=true
save_cross_sections=true
"""
        
        config_file = self.work_dir / "input_Ay.txt"
        with open(config_file, 'w') as f:
            f.write(config_content)
        
        print(f"配置文件生成: {config_file}")
        return config_file
    
    def generate_energy_file(self):
        """生成能量文件"""
        energy_file = Path("data/energy_190MeV_deuteron.txt")
        energy_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(energy_file, 'w') as f:
            f.write("# 氘核-质子散射能量点\n")
            f.write("# 质心系能量 (MeV)\n")
            f.write(f"{self.cm_energy:.4f}\n")
        
        print(f"能量文件生成: {energy_file}")
        return energy_file
    
    def generate_angle_file(self):
        """生成角度文件"""
        angle_file = Path("data/angles_cm.txt")
        angle_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(angle_file, 'w') as f:
            f.write("# 质心系散射角度 (度)\n")
            f.write("# θ_cm 范围: 10° - 170°\n")
            for theta in self.theta_cm:
                f.write(f"{theta:.2f}\n")
        
        print(f"角度文件生成: {angle_file}")
        return angle_file
    
    def run_calculation(self):
        """运行Tic-tac计算"""
        print("\n=== 开始氘核-质子Ay计算 ===")
        
        # 生成必要文件
        config_file = self.generate_config_file()
        self.generate_energy_file()
        self.generate_angle_file()
        
        # 检查tic-tac程序是否存在
        if not Path("tic-tac").exists():
            print("错误：tic-tac程序未找到，请先编译程序")
            print("运行命令：make")
            return False
        
        # 运行计算
        try:
            print("正在运行 Tic-tac 计算...")
            print("注意：这可能需要几个小时的计算时间")
            
            cmd = ["./tic-tac", "--input", str(config_file)]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600*4)  # 4小时超时
            
            if result.returncode == 0:
                print("✓ 计算成功完成")
                return True
            else:
                print(f"✗ 计算失败，错误代码：{result.returncode}")
                print("标准输出：", result.stdout)
                print("错误输出：", result.stderr)
                return False
                
        except subprocess.TimeoutExpired:
            print("✗ 计算超时（超过4小时）")
            return False
        except Exception as e:
            print(f"✗ 运行过程中出现错误：{e}")
            return False
    
    def analyze_results(self):
        """分析计算结果"""
        print("\n=== 分析计算结果 ===")
        
        results_dir = self.work_dir / "results"
        if not results_dir.exists():
            print("错误：结果目录不存在")
            return
        
        # 查找Ay结果文件
        ay_files = list(results_dir.glob("*Ay*.txt"))
        cross_section_files = list(results_dir.glob("*cross_section*.txt"))
        
        if ay_files:
            print(f"找到 {len(ay_files)} 个Ay结果文件：")
            for f in ay_files:
                print(f"  - {f.name}")
                
        if cross_section_files:
            print(f"找到 {len(cross_section_files)} 个截面结果文件：")
            for f in cross_section_files:
                print(f"  - {f.name}")
        
        # 尝试读取和显示部分结果
        self._display_sample_results(ay_files)
    
    def _display_sample_results(self, ay_files):
        """显示示例结果"""
        if not ay_files:
            return
        
        try:
            # 读取第一个Ay文件
            ay_file = ay_files[0]
            print(f"\n显示文件 {ay_file.name} 的部分结果：")
            
            with open(ay_file, 'r') as f:
                lines = f.readlines()
            
            # 显示前10行数据
            print("角度(度)     Ay")
            print("-" * 20)
            for i, line in enumerate(lines[:10]):
                if not line.startswith('#') and line.strip():
                    print(line.strip())
            
            if len(lines) > 10:
                print("...")
                print(f"（总共 {len(lines)} 行数据）")
                
        except Exception as e:
            print(f"读取结果文件时出错：{e}")
    
    def create_analysis_script(self):
        """创建结果分析脚本"""
        script_content = '''#!/usr/bin/env python3
"""
氘核-质子Ay结果分析脚本
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def plot_analyzing_power():
    """绘制分析能力曲线"""
    
    # 读取计算结果
    results_dir = Path("output/deuteron_proton_Ay/results")
    ay_files = list(results_dir.glob("*Ay*.txt"))
    
    if not ay_files:
        print("未找到Ay结果文件")
        return
    
    # 设置绘图
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # 绘制Ay
    for ay_file in ay_files:
        try:
            data = np.loadtxt(ay_file, comments='#')
            if data.ndim == 2 and data.shape[1] >= 2:
                angles = data[:, 0]
                ay_values = data[:, 1]
                
                ax1.plot(angles, ay_values, 'o-', label=ay_file.stem)
                ax1.set_xlabel('散射角 θ_cm (度)')
                ax1.set_ylabel('分析能力 Ay')
                ax1.set_title('氘核-质子散射 Ay (190 MeV/u)')
                ax1.grid(True, alpha=0.3)
                ax1.legend()
        
        except Exception as e:
            print(f"处理文件 {ay_file} 时出错: {e}")
    
    # 绘制微分截面（如果有的话）
    cross_section_files = list(results_dir.glob("*cross_section*.txt"))
    for cs_file in cross_section_files:
        try:
            data = np.loadtxt(cs_file, comments='#')
            if data.ndim == 2 and data.shape[1] >= 2:
                angles = data[:, 0]
                cross_sections = data[:, 1]
                
                ax2.semilogy(angles, cross_sections, 's-', label=cs_file.stem)
                ax2.set_xlabel('散射角 θ_cm (度)')
                ax2.set_ylabel('微分截面 (mb/sr)')
                ax2.set_title('氘核-质子散射微分截面')
                ax2.grid(True, alpha=0.3)
                ax2.legend()
        
        except Exception as e:
            print(f"处理文件 {cs_file} 时出错: {e}")
    
    plt.tight_layout()
    plt.savefig('deuteron_proton_Ay_190MeV.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("图形已保存为: deuteron_proton_Ay_190MeV.png")

def compare_with_experiment():
    """与实验数据比较（如果有的话）"""
    print("实验数据比较功能待实现...")
    print("建议查阅相关实验论文获取190 MeV/u氘核-質子散射的Ay实验数据")

if __name__ == "__main__":
    plot_analyzing_power()
    compare_with_experiment()
'''
        
        script_file = self.work_dir / "analyze_Ay_results.py"
        with open(script_file, 'w') as f:
            f.write(script_content)
        
        # 添加执行权限
        os.chmod(script_file, 0o755)
        print(f"分析脚本生成: {script_file}")
        
    def create_readme(self):
        """创建说明文档"""
        readme_content = f"""# 氘核-质子散射分析能力(Ay)计算

## 物理背景

本计算研究190 MeV/u极化氘核与质子散射的y方向分析能力：

- **入射粒子**: 极化氘核
- **实验室系能量**: 190 MeV/u (每核子)
- **质心系能量**: {self.cm_energy:.2f} MeV
- **靶核**: 质子（静止）
- **观测量**: y方向分析能力 Ay(θ)

## 分析能力的物理意义

分析能力Ay定义为：
```
Ay = (σ↑ - σ↓) / (σ↑ + σ↓)
```

其中：
- σ↑: y方向向上极化的截面
- σ↓: y方向向下极化的截面

Ay反映了散射过程中的自旋-轨道耦合效应。

## 计算参数

- 三体角动量截断: J_max = 3
- 波包数量: Np = Nq = 80
- 势能模型: N3LO Idaho手征势
- 包含张量力和同位旋破缺效应

## 运行方法

1. 确保已编译tic-tac程序：
   ```bash
   make
   ```

2. 运行计算：
   ```python
   python3 examples/deuteron_proton_Ay.py
   ```

3. 分析结果：
   ```python
   python3 output/deuteron_proton_Ay/analyze_Ay_results.py
   ```

## 预期结果

- Ay(θ)在前向角度可能为正值
- 在某些角度处可能出现零点
- 后向角度的行为取决于具体的核力模型

## 实验比较

建议与以下实验数据比较：
- RIKEN的氘核散射实验
- 相关的极化散射测量
- 理论计算的其他结果

## 注意事项

- 这是一个高精度计算，可能需要几个小时
- 结果的精度依赖于波包数量和角动量截断
- 极化计算对数值精度要求较高

生成时间: {self._get_current_time()}
"""
        
        readme_file = self.work_dir / "README.md"
        with open(readme_file, 'w') as f:
            f.write(readme_content)
        
        print(f"说明文档生成: {readme_file}")
    
    def _get_current_time(self):
        """获取当前时间"""
        from datetime import datetime
        return datetime.now().strftime("%Y年%m月%d日 %H:%M")

def main():
    """主函数"""
    print("=" * 60)
    print("氘核-质子散射分析能力(Ay)计算示例")
    print("能量: 190 MeV/u")
    print("=" * 60)
    
    # 创建计算实例
    calc = DeuteronProtonAyCalculation()
    
    # 生成所有必要文件
    calc.create_readme()
    calc.create_analysis_script()
    
    # 询问是否运行计算
    response = input("\n是否立即开始计算？(y/N): ").strip().lower()
    
    if response == 'y':
        success = calc.run_calculation()
        if success:
            calc.analyze_results()
            print("\n计算完成！请查看结果目录和分析脚本。")
        else:
            print("\n计算失败，请检查配置和程序状态。")
    else:
        print("\n文件已准备就绪。")
        print("要开始计算，请运行:")
        print("python3 examples/deuteron_proton_Ay.py")
    
    print(f"\n所有文件位于: {calc.work_dir}")

if __name__ == "__main__":
    main()