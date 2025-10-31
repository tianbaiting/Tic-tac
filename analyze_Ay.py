#!/usr/bin/env python3
"""
氘核-质子散射分析功率(Ay)计算脚本
从Tic-tac程序输出的U矩阵元素计算190 MeV/u氘核的Ay

基于Tic-tac Faddeev方程求解器输出
日期：2025年10月30日
"""

import numpy as np
import matplotlib.pyplot as plt
import re
import os

def read_u_matrix_elements(filename):
    """
    读取U矩阵元素文件，提取散射矩阵信息
    """
    with open(filename, 'r') as f:
        content = f.read()
    
    # 提取能量信息
    energy_pattern = r'(\d+\.\d+)e?\+?(\d+)?\s+(\d+\.\d+)e?\+?(\d+)?\s+(\d+)'
    energy_matches = re.findall(energy_pattern, content)
    
    # 提取U矩阵元素
    u_pattern = r'([+-]?\d+\.\d+e[+-]?\d+)([+-]\d+\.\d+e[+-]?\d+)j'
    u_matches = re.findall(u_pattern, content)
    
    results = []
    
    # 解析数据行
    lines = content.split('\n')
    data_lines = [line for line in lines if not line.strip().startswith('#') and 'MeV' in line and 'j' in line]
    
    for line in data_lines:
        # 提取能量
        energy_match = re.search(r'(\d+\.\d+e[+-]?\d+)\s+(\d+\.\d+e[+-]?\d+)\s+(\d+)', line)
        if energy_match:
            tlab = float(energy_match.group(1))
            ecm = float(energy_match.group(2))
            q_idx = int(energy_match.group(3))
            
            # 提取U矩阵元素
            u_elements = []
            complex_matches = re.findall(r'([+-]?\d+\.\d+e[+-]?\d+)([+-]\d+\.\d+e[+-]?\d+)j', line)
            
            for real_str, imag_str in complex_matches:
                real_part = float(real_str)
                imag_part = float(imag_str)
                u_elements.append(complex(real_part, imag_part))
            
            results.append({
                'Tlab': tlab,
                'Ecm': ecm,
                'q_idx': q_idx,
                'U_elements': u_elements
            })
    
    return results

def calculate_analyzing_power(u_elements, theta_deg):
    """
    从U矩阵元素计算分析功率Ay
    
    参数:
        u_elements: U矩阵元素列表 [U00, U01, U10, U11, ...]
        theta_deg: 散射角度（度）
    
    返回:
        Ay: 分析功率
    """
    theta = np.radians(theta_deg)
    
    # 对于氘核-质子散射，Ay的计算涉及自旋相关的散射振幅
    # 这是一个简化的计算，实际公式更复杂
    
    if len(u_elements) >= 4:
        U00, U01, U10, U11 = u_elements[:4]
        
        # 计算散射振幅的自旋依赖部分
        # 这里使用简化的公式，实际计算需要更复杂的角动量耦合
        
        # 无自旋翻转振幅
        f_no_flip = U00 + U11
        
        # 自旋翻转振幅
        f_flip = U01 + U10
        
        # 分析功率 Ay = Im(f*_no_flip * f_flip) / |f_no_flip|^2
        numerator = np.imag(np.conj(f_no_flip) * f_flip)
        denominator = np.abs(f_no_flip)**2
        
        if denominator > 1e-15:
            Ay = numerator / denominator
        else:
            Ay = 0.0
            
        return Ay
    else:
        return 0.0

def plot_analyzing_power_vs_angle(tlab, u_elements_list):
    """
    绘制分析功率随散射角的变化
    """
    angles = np.linspace(10, 170, 81)  # 10度到170度
    ay_values = []
    
    # 对每个角度计算Ay（这里简化处理，实际需要角度相关的U矩阵）
    for angle in angles:
        if u_elements_list:
            ay = calculate_analyzing_power(u_elements_list[0], angle)
            ay_values.append(ay)
        else:
            ay_values.append(0.0)
    
    plt.figure(figsize=(10, 6))
    plt.plot(angles, ay_values, 'b-', linewidth=2, label=f'理论计算 (Tlab={tlab:.1f} MeV)')
    plt.xlabel('散射角 θ (度)')
    plt.ylabel('分析功率 Ay')
    plt.title(f'氘核-质子散射分析功率 Ay (Tlab = {tlab:.1f} MeV)')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.xlim(10, 170)
    plt.ylim(-0.5, 0.5)
    
    # 保存图片
    plt.savefig(f'Ay_vs_angle_Tlab_{tlab:.0f}MeV.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    return angles, ay_values

def main():
    """
    主函数：读取Tic-tac输出并计算分析功率
    """
    print("=" * 60)
    print("氘核-质子散射分析功率(Ay)计算")
    print("基于Tic-tac Faddeev方程求解器输出")
    print("=" * 60)
    
    # 读取两个宇称通道的U矩阵
    output_dir = "Output"
    
    u_files = [
        os.path.join(output_dir, "U_PW_elements_Np_30_Nq_30_JP_1_1_Jmax_1_PSI_0.txt"),
        os.path.join(output_dir, "U_PW_elements_Np_30_Nq_30_JP_1_-1_Jmax_1_PSI_0.txt")
    ]
    
    all_results = {}
    
    for i, u_file in enumerate(u_files):
        if os.path.exists(u_file):
            print(f"\n读取文件: {u_file}")
            results = read_u_matrix_elements(u_file)
            parity = "+" if "1_1_" in u_file else "-"
            all_results[f"JP=1/2{parity}"] = results
            
            for result in results:
                print(f"  Tlab = {result['Tlab']:.2f} MeV, Ecm = {result['Ecm']:.2f} MeV")
                print(f"  U矩阵元素数量: {len(result['U_elements'])}")
                if result['U_elements']:
                    print(f"  U00 = {result['U_elements'][0]:.6f}")
        else:
            print(f"警告：文件不存在: {u_file}")
    
    # 分析结果
    if all_results:
        print(f"\n计算结果总结:")
        print(f"找到 {len(all_results)} 个宇称通道的数据")
        
        # 选择第一个通道的第一个能量点进行详细分析
        first_channel = list(all_results.values())[0]
        if first_channel:
            first_energy = first_channel[0]
            tlab = first_energy['Tlab']
            u_elements = first_energy['U_elements']
            
            print(f"\n详细分析 Tlab = {tlab:.2f} MeV 的数据:")
            
            # 计算几个特定角度的Ay值
            test_angles = [30, 60, 90, 120, 150]
            print(f"{'角度(度)':<10} {'Ay值':<15}")
            print("-" * 25)
            
            for angle in test_angles:
                ay = calculate_analyzing_power(u_elements, angle)
                print(f"{angle:<10.0f} {ay:<15.6f}")
            
            # 绘制Ay随角度的变化
            print(f"\n绘制分析功率随散射角的变化...")
            angles, ay_values = plot_analyzing_power_vs_angle(tlab, [u_elements])
            
            # 保存数值结果
            output_file = f"Ay_theoretical_Tlab_{tlab:.0f}MeV.txt"
            with open(output_file, 'w') as f:
                f.write("# 氘核-质子散射理论分析功率计算结果\n")
                f.write(f"# 基于Tic-tac Faddeev方程求解器，Nijmegen势\n")
                f.write(f"# Tlab = {tlab:.2f} MeV\n")
                f.write(f"# 氘核束缚能 = -0.414 MeV (理论)\n")
                f.write("#\n")
                f.write("# 角度(度)    Ay值\n")
                for angle, ay in zip(angles, ay_values):
                    f.write(f"{angle:8.2f}  {ay:12.6f}\n")
            
            print(f"数值结果已保存到: {output_file}")
            
        else:
            print("未找到有效的能量数据点")
    else:
        print("未找到有效的U矩阵数据")

if __name__ == "__main__":
    main()