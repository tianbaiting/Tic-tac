#!/usr/bin/env python3
"""
同步CPP目录到重构后的src目录结构
保持文件内容同步，同时维护重构后的目录组织
"""

import os
import shutil
from pathlib import Path

# 文件映射：CPP目录 -> src目录
FILE_MAPPING = {
    # 主程序
    "CPP/main.cpp": "src/main.cpp",
    
    # 配置相关
    "CPP/run_organizer.cpp": "src/config/run_organizer.cpp",
    "CPP/run_organizer.h": "src/config/run_organizer.h",
    "CPP/set_run_parameters.cpp": "src/config/set_run_parameters.cpp",
    "CPP/set_run_parameters.h": "src/config/set_run_parameters.h",
    
    # IO相关
    "CPP/disk_io_routines.cpp": "src/io/disk_io_routines.cpp",
    "CPP/disk_io_routines.h": "src/io/disk_io_routines.h",
    
    # 工具函数
    "CPP/auxiliary.cpp": "src/utils/auxiliary.cpp",
    "CPP/auxiliary.h": "src/utils/auxiliary.h",
    "CPP/error_management.cpp": "src/utils/error_management.cpp",
    "CPP/error_management.h": "src/utils/error_management.h",
    
    # 核心 - 状态空间
    "CPP/make_pw_symm_states.cpp": "src/core/state_space/make_pw_symm_states.cpp",
    "CPP/make_pw_symm_states.h": "src/core/state_space/make_pw_symm_states.h",
    "CPP/make_wp_states.cpp": "src/core/state_space/make_wp_states.cpp",
    "CPP/make_wp_states.h": "src/core/state_space/make_wp_states.h",
    "CPP/make_swp_states.cpp": "src/core/state_space/make_swp_states.cpp",
    "CPP/make_swp_states.h": "src/core/state_space/make_swp_states.h",
    "CPP/make_permutation_matrix.cpp": "src/core/state_space/make_permutation_matrix.cpp",
    "CPP/make_permutation_matrix.h": "src/core/state_space/make_permutation_matrix.h",
    
    # 核心 - 势矩阵
    "CPP/make_potential_matrix.cpp": "src/core/potential/make_potential_matrix.cpp",
    "CPP/make_potential_matrix.h": "src/core/potential/make_potential_matrix.h",
    
    # 核心 - Resolvent
    "CPP/make_resolvent.cpp": "src/core/resolvent/make_resolvent.cpp",
    "CPP/make_resolvent.h": "src/core/resolvent/make_resolvent.h",
    
    # 核心 - Faddeev求解器
    "CPP/solve_faddeev.cpp": "src/core/faddeev_solver/solve_faddeev.cpp",
    "CPP/solve_faddeev.h": "src/core/faddeev_solver/solve_faddeev.h",
    
    # 头文件
    "CPP/constants.h": "include/constants.h",
    "CPP/type_defs.h": "include/type_defs.h",
}

# General_functions映射
GENERAL_FUNCTIONS_MAPPING = {
    "CPP/General_functions/coupling_coefficients.cpp": "src/utils/coupling_coefficients.cpp",
    "CPP/General_functions/coupling_coefficients.h": "src/utils/coupling_coefficients.h",
    "CPP/General_functions/gauss_legendre.cpp": "src/utils/gauss_legendre.cpp",
    "CPP/General_functions/gauss_legendre.h": "src/utils/gauss_legendre.h",
    "CPP/General_functions/kinetic_conversion.cpp": "src/utils/kinetic_conversion.cpp",
    "CPP/General_functions/kinetic_conversion.h": "src/utils/kinetic_conversion.h",
    "CPP/General_functions/legendre.cpp": "src/utils/legendre.cpp",
    "CPP/General_functions/legendre.h": "src/utils/legendre.h",
    "CPP/General_functions/matrix_routines.cpp": "src/utils/matrix_routines.cpp",
    "CPP/General_functions/matrix_routines.h": "src/utils/matrix_routines.h",
    "CPP/General_functions/templates.h": "src/utils/templates.h",
}

FILE_MAPPING.update(GENERAL_FUNCTIONS_MAPPING)

def fix_includes(content, filepath):
    """修复include路径"""
    lines = content.split('\n')
    fixed_lines = []
    
    for line in lines:
        new_line = line
        
        # 修复General_functions路径
        if 'General_functions/' in line:
            new_line = line.replace('General_functions/', 'utils/')
        
        # 修复相对路径
        if '../' in new_line:
            # 简单处理：移除所有../
            import re
            new_line = re.sub(r'#include\s+"\.\./', '#include "', new_line)
        
        # 对于src目录下的文件，确保不包含include/前缀
        if 'include/' in new_line and '/src/' in filepath:
            new_line = new_line.replace('include/', '')
        
        fixed_lines.append(new_line)
    
    return '\n'.join(fixed_lines)

def sync_files():
    """同步文件从CPP到src"""
    print("开始同步CPP目录到src目录...")
    synced_count = 0
    
    for cpp_path, src_path in FILE_MAPPING.items():
        if os.path.exists(cpp_path):
            # 创建目标目录
            os.makedirs(os.path.dirname(src_path), exist_ok=True)
            
            # 读取源文件
            with open(cpp_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            # 修复include路径
            fixed_content = fix_includes(content, src_path)
            
            # 写入目标文件
            with open(src_path, 'w', encoding='utf-8') as f:
                f.write(fixed_content)
            
            print(f"✓ {cpp_path} -> {src_path}")
            synced_count += 1
        else:
            print(f"✗ 源文件不存在: {cpp_path}")
    
    print(f"\n同步完成！共处理 {synced_count} 个文件")

if __name__ == "__main__":
    sync_files()
