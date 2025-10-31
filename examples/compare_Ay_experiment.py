#!/usr/bin/env python3
"""
氘核-质子Ay实验数据对比脚本
============================

这个脚本用于将Tic-tac计算的分析能力结果与已发表的实验数据进行对比。

实验数据来源：
1. RIKEN实验室的极化氘核-质子散射测量
2. TRIUMF的精密散射实验
3. 其他国际实验室的相关测量

作者：Tic-tac团队
日期：2024年10月
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json

class AyDataComparison:
    """Ay数据对比类"""
    
    def __init__(self):
        """初始化"""
        self.experimental_data = {}
        self.theoretical_data = {}
        self.load_experimental_data()
    
    def load_experimental_data(self):
        """加载实验数据（模拟数据，实际使用时需要真实数据）"""
        
        # 模拟的190 MeV/u氘核-质子散射Ay实验数据
        # 实际使用时应该从实验论文或数据库获取真实数据
        
        angles_exp = np.array([
            20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160
        ])  # 质心系角度（度）
        
        # 模拟的Ay值（基于典型的d+p散射行为）
        ay_exp = np.array([
            0.02, 0.04, 0.08, 0.12, 0.15, 0.16, 0.14, 0.10, 0.05, 0.00, 
            -0.05, -0.08, -0.06, -0.02, 0.01
        ])
        
        # 实验误差（模拟）
        ay_err = np.array([
            0.01, 0.01, 0.015, 0.02, 0.02, 0.025, 0.02, 0.015, 0.02, 0.025,
            0.02, 0.025, 0.02, 0.015, 0.015
        ])
        
        self.experimental_data = {
            'angles': angles_exp,
            'ay_values': ay_exp,
            'ay_errors': ay_err,
            'energy': '190 MeV/u',
            'reference': '模拟数据 - 请替换为真实实验数据',
            'note': '基于典型d+p散射Ay行为的模拟数据'
        }
        
        print("已加载实验数据（模拟）")
        print("注意：请将此模拟数据替换为真实的实验数据")
    
    def load_theoretical_results(self, results_dir="output/deuteron_proton_Ay/results"):
        """加载理论计算结果"""
        results_path = Path(results_dir)
        
        if not results_path.exists():
            print(f"错误：结果目录不存在 - {results_path}")
            return False
        
        # 查找Ay结果文件
        ay_files = list(results_path.glob("*Ay*.txt"))
        
        if not ay_files:
            print("未找到Ay计算结果文件")
            return False
        
        # 读取第一个Ay文件
        ay_file = ay_files[0]
        print(f"读取理论结果: {ay_file}")
        
        try:
            data = np.loadtxt(ay_file, comments='#')
            if data.ndim == 2 and data.shape[1] >= 2:
                self.theoretical_data = {
                    'angles': data[:, 0],
                    'ay_values': data[:, 1],
                    'file': ay_file.name
                }
                print(f"成功读取 {len(data)} 个数据点")
                return True
            else:
                print("数据格式不正确")
                return False
                
        except Exception as e:
            print(f"读取文件时出错: {e}")
            return False
    
    def plot_comparison(self, save_plot=True):
        """绘制理论与实验对比图"""
        
        if not self.theoretical_data:
            print("缺少理论计算数据")
            return
        
        # 创建图形
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # 绘制实验数据
        exp_data = self.experimental_data
        ax.errorbar(exp_data['angles'], exp_data['ay_values'], 
                   yerr=exp_data['ay_errors'],
                   fmt='ro', capsize=5, capthick=2, 
                   label=f"实验数据 ({exp_data['reference']})", 
                   markersize=6)
        
        # 绘制理论计算
        theo_data = self.theoretical_data
        ax.plot(theo_data['angles'], theo_data['ay_values'], 
               'b-', linewidth=2, label='Tic-tac 计算')
        ax.plot(theo_data['angles'], theo_data['ay_values'], 
               'bs', markersize=4)
        
        # 设置图形属性
        ax.set_xlabel('散射角 θ_cm (度)', fontsize=12)
        ax.set_ylabel('分析能力 Ay', fontsize=12)
        ax.set_title('氘核-质子散射分析能力对比 (190 MeV/u)', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=11)
        
        # 设置坐标轴范围
        ax.set_xlim(0, 180)
        ax.axhline(y=0, color='k', linestyle='--', alpha=0.5)
        
        # 添加说明文本
        ax.text(0.02, 0.98, f"质心系能量: {exp_data['energy']}", 
               transform=ax.transAxes, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        
        if save_plot:
            filename = 'deuteron_proton_Ay_comparison_190MeV.png'
            plt.savefig(filename, dpi=300, bbox_inches='tight')
            print(f"对比图已保存: {filename}")
        
        plt.show()
    
    def calculate_chi_squared(self):
        """计算卡方值评估理论与实验的符合程度"""
        
        if not self.theoretical_data:
            print("缺少理论计算数据，无法计算卡方")
            return
        
        exp_angles = self.experimental_data['angles']
        exp_ay = self.experimental_data['ay_values']
        exp_err = self.experimental_data['ay_errors']
        
        theo_angles = self.theoretical_data['angles']
        theo_ay = self.theoretical_data['ay_values']
        
        # 插值理论数据到实验角度点
        theo_ay_interp = np.interp(exp_angles, theo_angles, theo_ay)
        
        # 计算卡方
        chi_squared = np.sum(((exp_ay - theo_ay_interp) / exp_err)**2)
        dof = len(exp_angles)  # 自由度
        reduced_chi_squared = chi_squared / dof
        
        print("\n=== 统计分析 ===")
        print(f"数据点数: {len(exp_angles)}")
        print(f"卡方值 χ²: {chi_squared:.2f}")
        print(f"自由度: {dof}")
        print(f"约化卡方 χ²/dof: {reduced_chi_squared:.2f}")
        
        if reduced_chi_squared < 1.5:
            print("✓ 理论与实验符合良好")
        elif reduced_chi_squared < 3.0:
            print("⚠ 理论与实验符合尚可，可能需要改进")
        else:
            print("✗ 理论与实验差异较大，需要检查计算或模型")
        
        return chi_squared, reduced_chi_squared
    
    def generate_report(self, output_file="Ay_comparison_report.txt"):
        """生成对比报告"""
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write("氘核-质子散射分析能力(Ay)对比报告\n")
            f.write("=" * 50 + "\n\n")
            
            f.write("实验数据信息:\n")
            f.write(f"  能量: {self.experimental_data['energy']}\n")
            f.write(f"  数据点数: {len(self.experimental_data['angles'])}\n")
            f.write(f"  角度范围: {self.experimental_data['angles'][0]}° - {self.experimental_data['angles'][-1]}°\n")
            f.write(f"  参考文献: {self.experimental_data['reference']}\n\n")
            
            if self.theoretical_data:
                f.write("理论计算信息:\n")
                f.write(f"  计算文件: {self.theoretical_data['file']}\n")
                f.write(f"  数据点数: {len(self.theoretical_data['angles'])}\n")
                f.write(f"  角度范围: {self.theoretical_data['angles'][0]:.1f}° - {self.theoretical_data['angles'][-1]:.1f}°\n\n")
                
                # 计算统计量
                chi2, reduced_chi2 = self.calculate_chi_squared()
                f.write("统计分析:\n")
                f.write(f"  卡方值: {chi2:.2f}\n")
                f.write(f"  约化卡方: {reduced_chi2:.2f}\n\n")
            
            f.write("详细数据对比:\n")
            f.write("角度(度)  实验Ay   误差    理论Ay   差异\n")
            f.write("-" * 45 + "\n")
            
            exp_angles = self.experimental_data['angles']
            exp_ay = self.experimental_data['ay_values']
            exp_err = self.experimental_data['ay_errors']
            
            if self.theoretical_data:
                theo_angles = self.theoretical_data['angles']
                theo_ay = self.theoretical_data['ay_values']
                theo_ay_interp = np.interp(exp_angles, theo_angles, theo_ay)
                
                for i, angle in enumerate(exp_angles):
                    diff = abs(exp_ay[i] - theo_ay_interp[i])
                    f.write(f"{angle:6.1f}   {exp_ay[i]:7.3f}  {exp_err[i]:6.3f}  {theo_ay_interp[i]:7.3f}  {diff:6.3f}\n")
        
        print(f"对比报告已生成: {output_file}")

def main():
    """主函数"""
    print("氘核-质子散射Ay数据对比工具")
    print("=" * 40)
    
    # 创建对比实例
    comparison = AyDataComparison()
    
    # 尝试加载理论结果
    success = comparison.load_theoretical_results()
    
    if success:
        # 绘制对比图
        comparison.plot_comparison()
        
        # 计算统计量
        comparison.calculate_chi_squared()
        
        # 生成报告
        comparison.generate_report()
        
        print("\n对比分析完成！")
        print("请查看生成的图形和报告文件。")
    else:
        print("\n无法找到理论计算结果。")
        print("请先运行: python3 examples/deuteron_proton_Ay.py")
    
    print("\n重要提醒：")
    print("1. 当前使用的是模拟实验数据")
    print("2. 请将模拟数据替换为真实的实验数据")
    print("3. 真实数据可从相关实验论文或数据库获取")

if __name__ == "__main__":
    main()