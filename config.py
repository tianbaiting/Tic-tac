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
            # 与 set_default_values(...) 对齐的默认值
            'two_J_3N_max': 1,
            'Np_WP': 50,
            'Nq_WP': 50,
            'J_2N_max': 3,
            'Nphi': 48,
            'Nx': 48,
            'Np_per_WP': 8,
            'Nq_per_WP': 8,
            'chebyshev_s': 100,
            'chebyshev_t': 1,
            'channel_idx': -1,
            'parallel_run': False,
            'P123_omp_num_threads': 4,
            'P123_recovery': False,
            'tensor_force': True,
            'isospin_breaking_1S0': True,
            'midpoint_approx': False,
            'calculate_and_store_P123': True,
            'include_breakup_channels': False,
            'solve_faddeev': True,
            'solve_dense': False,
            'production_run': True,
            'potential_model': 'LO_internal',
            'parameter_walk': False,
            'parameter_file': 'none',
            'PSI_start': -1,
            'PSI_end': -1,
            'subfolder': 'Output',
            'p_grid_type': 'chebyshev',
            'p_grid_filename': '',
            'q_grid_type': 'chebyshev',
            'q_grid_filename': '',
            'energy_input_file': 'CPP/Input/lab_energies.txt',
            'output_folder': 'CPP/Output',
            'P123_folder': 'CPP/Output',
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
                    else:
                        try:
                            self.config[key] = int(value)
                        except ValueError:
                            try:
                                self.config[key] = float(value)
                            except ValueError:
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
        directory = os.path.dirname(filename)
        if directory:
            os.makedirs(directory, exist_ok=True)
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
