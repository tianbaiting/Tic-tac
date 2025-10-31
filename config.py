#!/usr/bin/env python3
"""
Tic-tac 计算参数配置脚本
用于生成和管理计算参数文件
"""

import os
import sys
from pathlib import Path

class TicTacConfig:
    """Tic-tac计算配置类"""
    
    def __init__(self):
        """初始化默认配置"""
        self.config = {
            # 基本参数
            'two_J_3N_max': 1,          # 三体系统最大总角动量（×2）
            'Np_WP': 50,                # p方向波包数量
            'Nq_WP': 50,                # q方向波包数量
            'J_2N_max': 1,              # 二体系统最大角动量
            
            # 网格参数
            'Nphi': 48,                 # φ角积分点数
            'Nx': 48,                   # x积分点数
            'Np_per_WP': 8,             # 每个波包的p点数
            'Nq_per_WP': 8,             # 每个波包的q点数
            
            # Chebyshev网格参数
            'chebyshev_s': 200,         # Chebyshev网格缩放参数
            'chebyshev_t': 1,           # Chebyshev网格类型参数
            
            # 网格类型
            'p_grid_type': 'chebyshev', # p网格类型
            'q_grid_type': 'chebyshev', # q网格类型
            
            # 网格文件
            'p_grid_filename': 'data/NWP-100-splitnorm.txt',
            'q_grid_filename': 'data/NWP-100-splitnorm.txt',
            
            # 并行设置
            'P123_omp_num_threads': 4,  # OpenMP线程数
            'parallel_run': True,       # 是否并行运行
            
            # 计算选项
            'P123_recovery': False,     # 是否恢复P123计算
            'tensor_force': True,       # 是否包含张量力
            'isospin_breaking_1S0': True, # 是否包含同位旋破缺
            'midpoint_approx': False,   # 是否使用中点近似
            'calculate_and_store_P123': True, # 是否计算并存储P123
            'include_breakup_channels': False, # 是否包含破裂道
            'solve_faddeev': True,      # 是否求解Faddeev方程
            'solve_dense': False,       # 是否使用稠密求解器
            'production_run': True,     # 是否生产运行
            
            # 势能模型
            'potential_model': 'LO_internal', # 势能模型类型
            
            # 参数扫描
            'parameter_walk': False,    # 是否参数扫描
            'parameter_file': '',       # 参数文件
            
            # 能量范围
            'PSI_start': -1,           # 起始PSI
            'PSI_end': -1,             # 结束PSI
            
            # 文件路径
            'energy_input_file': 'data/lab_energies.txt',
            'output_folder': 'output',
            'P123_folder': 'data/permutation_matrices'
        }
    
    def save_config(self, filename='data/input.txt'):
        """保存配置到文件"""
        with open(filename, 'w', encoding='utf-8') as f:
            f.write("# Tic-tac 计算参数配置文件\n")
            f.write("# 自动生成 - 请根据需要修改\n\n")
            
            for key, value in self.config.items():
                if isinstance(value, bool):
                    f.write(f"{key}={str(value).lower()}\n")
                elif isinstance(value, str) and value == '':
                    f.write(f"{key}=\n")
                else:
                    f.write(f"{key}={value}\n")
        
        print(f"配置文件已保存到: {filename}")
    
    def load_config(self, filename):
        """从文件加载配置"""
        if not os.path.exists(filename):
            print(f"配置文件不存在: {filename}")
            return False
        
        with open(filename, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                
                if '=' in line:
                    key, value = line.split('=', 1)
                    key = key.strip()
                    value = value.strip()
                    
                    # 类型转换
                    if value.lower() in ['true', 'false']:
                        self.config[key] = value.lower() == 'true'
                    elif value.isdigit():
                        self.config[key] = int(value)
                    elif value.replace('.', '').isdigit():
                        self.config[key] = float(value)
                    else:
                        self.config[key] = value
        
        print(f"配置已从文件加载: {filename}")
        return True
    
    def set_parameter(self, key, value):
        """设置单个参数"""
        self.config[key] = value
        print(f"参数已设置: {key} = {value}")
    
    def get_parameter(self, key):
        """获取单个参数"""
        return self.config.get(key, None)
    
    def print_config(self):
        """打印当前配置"""
        print("=== Tic-tac 当前配置 ===")
        for key, value in self.config.items():
            print(f"{key:25s} = {value}")
        print("========================")

def main():
    """主函数"""
    config = TicTacConfig()
    
    if len(sys.argv) < 2:
        print("用法:")
        print("  python3 config.py save [filename]    - 保存默认配置")
        print("  python3 config.py load [filename]    - 加载配置")
        print("  python3 config.py show               - 显示当前配置")
        return
    
    command = sys.argv[1]
    
    if command == 'save':
        filename = sys.argv[2] if len(sys.argv) > 2 else 'data/input.txt'
        # 确保目录存在
        os.makedirs(os.path.dirname(filename), exist_ok=True)
        config.save_config(filename)
    
    elif command == 'load':
        filename = sys.argv[2] if len(sys.argv) > 2 else 'data/input.txt'
        config.load_config(filename)
        config.print_config()
    
    elif command == 'show':
        config.print_config()
    
    else:
        print(f"未知命令: {command}")

if __name__ == '__main__':
    main()