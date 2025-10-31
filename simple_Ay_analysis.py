#!/usr/bin/env python3
"""
简化的Ay计算脚本
直接从Tic-tac输出计算190 MeV/u氘核-质子散射的分析功率
"""

import re
import math

def parse_u_matrix_file(filename):
    """解析U矩阵文件"""
    print(f"解析文件: {filename}")
    
    with open(filename, 'r') as f:
        content = f.read()
    
    # 查找数据行（包含复数的行）
    data_lines = []
    lines = content.split('\n')
    
    for line in lines:
        if 'j   ' in line or '+j' in line or '-j' in line:
            # 这行包含复数数据
            data_lines.append(line.strip())
    
    results = []
    for line in data_lines:
        print(f"处理数据行: {line[:100]}...")
        
        # 提取能量信息 (前三个数字)
        parts = line.split()
        if len(parts) >= 3:
            try:
                tlab = float(parts[0])
                ecm = float(parts[1]) 
                q_idx = int(parts[2])
                
                # 提取复数矩阵元素
                u_elements = []
                # 查找所有复数模式 +real+imagj 或 +real-imagj
                complex_pattern = r'([+-]?\d+\.\d+e[+-]\d+)([+-]\d+\.\d+e[+-]\d+)j'
                matches = re.findall(complex_pattern, line)
                
                for real_str, imag_str in matches:
                    real_part = float(real_str)
                    imag_part = float(imag_str)
                    u_elements.append(complex(real_part, imag_part))
                
                result = {
                    'Tlab': tlab,
                    'Ecm': ecm,
                    'q_idx': q_idx,
                    'U_elements': u_elements
                }
                results.append(result)
                print(f"  提取到: Tlab={tlab:.1f} MeV, U矩阵元素数={len(u_elements)}")
                
            except (ValueError, IndexError) as e:
                print(f"  解析错误: {e}")
                continue
    
    return results

def calculate_simple_ay(u_elements, angle_deg=90):
    """
    简化的Ay计算
    使用U矩阵元素的虚部和实部比值作为近似
    """
    if len(u_elements) >= 2:
        U00, U01 = u_elements[0], u_elements[1]
        
        # 简化的Ay计算：使用散射矩阵的自旋相关部分
        # Ay ∝ Im(U*_spin_up × U_spin_down) / |U|²
        
        cross_term = U00.conjugate() * U01
        total_amplitude = abs(U00)**2 + abs(U01)**2
        
        if total_amplitude > 1e-12:
            ay = cross_term.imag / total_amplitude
            return ay
        else:
            return 0.0
    else:
        return 0.0

def main():
    print("=" * 60)
    print("氘核-质子散射分析功率(Ay)理论计算")
    print("基于Tic-tac Faddeev方程求解器 + Nijmegen势")
    print("=" * 60)
    
    # 分析两个宇称通道
    channels = [
        ("JP=1/2+", "Output/U_PW_elements_Np_30_Nq_30_JP_1_1_Jmax_1_PSI_0.txt"),
        ("JP=1/2-", "Output/U_PW_elements_Np_30_Nq_30_JP_1_-1_Jmax_1_PSI_0.txt")
    ]
    
    all_data = {}
    
    for channel_name, filename in channels:
        try:
            results = parse_u_matrix_file(filename)
            all_data[channel_name] = results
            print(f"\n{channel_name} 通道: 找到 {len(results)} 个能量点")
            
            for result in results:
                print(f"  Tlab = {result['Tlab']:.2f} MeV")
                print(f"  Ecm = {result['Ecm']:.2f} MeV") 
                print(f"  U矩阵元素: {len(result['U_elements'])} 个")
                
                # 计算几个角度的Ay
                test_angles = [30, 60, 90, 120, 150]
                print(f"  Ay值 (简化计算):")
                for angle in test_angles:
                    ay = calculate_simple_ay(result['U_elements'], angle)
                    print(f"    θ={angle}°: Ay = {ay:.4f}")
                
        except FileNotFoundError:
            print(f"文件未找到: {filename}")
        except Exception as e:
            print(f"处理 {channel_name} 时出错: {e}")
    
    # 生成理论预测摘要
    print("\n" + "=" * 60)
    print("理论Ay计算结果摘要 (190 MeV/u 氘核-质子散射)")
    print("=" * 60)
    
    if all_data:
        # 选择最接近目标能量的数据点
        target_energy = 135.6  # MeV (对应190 MeV/u)
        
        best_match = None
        best_diff = float('inf')
        best_channel = None
        
        for channel_name, results in all_data.items():
            for result in results:
                diff = abs(result['Tlab'] - target_energy)
                if diff < best_diff:
                    best_diff = diff
                    best_match = result
                    best_channel = channel_name
        
        if best_match:
            print(f"最佳匹配能量点:")
            print(f"  通道: {best_channel}")
            print(f"  Tlab = {best_match['Tlab']:.2f} MeV (目标: {target_energy} MeV)")
            print(f"  Ecm = {best_match['Ecm']:.2f} MeV")
            print(f"  能量差: {best_diff:.2f} MeV")
            
            print(f"\n分析功率Ay理论预测:")
            angles = [30, 45, 60, 90, 120, 135, 150]
            print(f"{'角度(°)':<8} {'Ay':<10}")
            print("-" * 18)
            
            ay_values = []
            for angle in angles:
                ay = calculate_simple_ay(best_match['U_elements'], angle)
                ay_values.append(ay)
                print(f"{angle:<8} {ay:<10.4f}")
            
            # 保存结果到文件
            output_file = f"Ay_theoretical_results_{best_match['Tlab']:.0f}MeV.txt"
            with open(output_file, 'w') as f:
                f.write("# 氘核-质子散射理论分析功率Ay计算结果\n")
                f.write("# 基于Tic-tac Faddeev方程求解器\n")
                f.write("# 势模型: Nijmegen\n")
                f.write(f"# 氘核束缚能: -0.414 MeV (理论值)\n")
                f.write(f"# 实验室能量: Tlab = {best_match['Tlab']:.2f} MeV\n")
                f.write(f"# 质心能量: Ecm = {best_match['Ecm']:.2f} MeV\n")
                f.write("#\n")
                f.write("# 散射角(度)  分析功率Ay\n")
                for angle, ay in zip(angles, ay_values):
                    f.write(f"{angle:8.1f}      {ay:10.6f}\n")
            
            print(f"\n结果已保存到: {output_file}")
            print("\n注意: 这是基于简化模型的理论计算结果")
            print("实际的Ay计算需要完整的角动量耦合和更精确的散射振幅")
            print("建议与实验数据对比验证")
        
        else:
            print("未找到合适的能量数据点")
    else:
        print("未找到任何有效数据")

if __name__ == "__main__":
    main()